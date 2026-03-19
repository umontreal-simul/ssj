package ift6561examples;

import umontreal.ssj.simevents.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.stat.*;
import java.io.*;
import java.util.*;

/*
 * Simulates a bank with two types of customers and two types of employees
 * (tellers and counselors) as in homework 3 in 2020 and 2023.
 * See the homework statement and 2023 solution for the specifications.  
 */
public class BankTwoTypes {
   static final double HOUR = 3600.0; // Time is in seconds.
   static final double MINUTE = 60.0;

   // Model data
   // All times are converted into seconds when reading them from the input file. 
   // After that, all the times in the code are in seconds.
   int numPeriods;     // Number of periods in the day.
   double lengthPeriod;// The length of each period (read in minutes).
   double[] lambda;    // Type A arrival rate lambda_j for each period j (read per hour).
   int[] numTellers;   // Number of tellers (type C) for each period.
   int numSlots;       // Number of appointment time slots.
   double slotLength;  // Length of appointment time slots (read in minutes).
   int numCounselors;  // Number of counselors (type D) for the day.
   double r;           // Probability of having an appointment, for any time slot.
   double p;           // Probability that customer does not show up for the appointment.
   double s;           // Threshold for counselors to become tellers (read in minutes).

   // Global variables
   double nExpectedA;  // Expected number of customers of type A per day.
   double nExpectedB;  // Expected number of customers of type B per day.
   int nArrivalsA;     // Number of arrivals of type A today.
   int nArrivalsB;     // Number of arrivals of type B today.
   double arrRateA;    // Current arrival rate for type A.
   int nTellers;       // Number of tellers (type C) in current period.
   int nTellersBusy;   // Number of tellers (type C) that are occupied.
   int currentSlot;    // The current appointment time slot (first one is 0).
   double nextSlotTime; // Start time of next appointment slot. 
   boolean counselorBusy[]; // For each counselor, indicates if he is occupied. 
   CustomerB custBWaiting[]; // The customer waiting, for each counselor, null if none. 
   boolean appoint[][]; // Indicators of appointments, for each counselor and slot.
   int numMatchesAD;    // Number of type A customers served by a counselor today.
   
   Event nextArrivalA = new ArrivalA(); // The next Arrival event for type A.
                        // We sometimes have to change its occurrence time.
   LinkedList<CustomerA> waitListA = new LinkedList<CustomerA>();

   RandomStream streamArrA = new MRG31k3p();    // For arrivals of type A.
   RandomStream streamAppoint = new MRG31k3p(); // For appointments and showing up.  
   LognormalGen genServA, genServB; // For service times, types A and B.
   NormalGen genDelayB;  // To generate arrival delays R for customers type B.
   RandomStreamManager myStreams = new RandomStreamManager(); // List of all random streams.
   
   Tally statArrivalsA = new Tally("Number of arrivals per day, type A");  // To check.
   Tally statArrivalsB = new Tally("Number of arrivals per day, type B");
   Tally statAD  = new Tally("Number of times a counselor (type D) serves a type A");
   Tally statWaitsA = new Tally("For average waiting times for customers type A");
   Tally statWaitsB = new Tally("For average waiting times for customers type B");
   Tally statWaitsInDayA = new Tally("Waiting times within a day for type A");   
   Tally statWaitsInDayB = new Tally("Waiting times within a day for type B");
   // Tally allWaitsA = new TallyStore("All waiting times, for all days, for type A");
   // Tally allWaitsB = new TallyStore("All Waiting times, for all days, for type B");

   
   // The constructor: reads data in file and creates generators after their param. are read.
   // Also put all the streams in a list, for easier management.
   public BankTwoTypes(String fileName) throws IOException {
      readData(fileName);
      // Add all the streams to my list of streams.
      myStreams.add(streamArrA);
      myStreams.add(streamAppoint);
      myStreams.add(genServA.getStream());
      myStreams.add(genServB.getStream());
      myStreams.add(genDelayB.getStream());
   }

   // Reads data, constructs some arrays, and creates generators after params are read.
   // To simplify the code, all times are converted into seconds right after reading them. 
   public void readData(String dataFile) throws IOException {
      Locale loc = Locale.getDefault();
      Locale.setDefault(Locale.US); // to read reals as 8.3 instead of 8,3
      BufferedReader input = new BufferedReader(new FileReader(dataFile));
      Scanner scan = new Scanner(input);
      numPeriods = scan.nextInt();   
      lengthPeriod = scan.nextDouble() * MINUTE;  scan.nextLine();
      lambda = new double[numPeriods];
      numTellers = new int[numPeriods];
      nExpectedA = 0.0;
      for (int j = 0; j < numPeriods; j++) {
         numTellers[j] = scan.nextInt();
         lambda[j] = scan.nextDouble() / HOUR;
         nExpectedA += lambda[j] * lengthPeriod;
         scan.nextLine();
      }
      numCounselors = scan.nextInt();   scan.nextLine();
      numSlots = scan.nextInt();
      slotLength = scan.nextDouble() * MINUTE;   scan.nextLine();
      appoint = new boolean[numCounselors][numSlots+1]; // We put one extra empty slot.
      counselorBusy = new boolean[numCounselors];       // Indicates if counselor is busy, for each.
      custBWaiting = new CustomerB[numCounselors];  
      r = scan.nextDouble();   p = scan.nextDouble();  scan.nextLine();
      s = scan.nextDouble() * MINUTE;   scan.nextLine();
      nExpectedB = numCounselors * numSlots * r * (1.0 - p);
            
      double mu, sigma; // Mean and stdev for normal and lognormal times.
      mu = scan.nextDouble();    sigma = scan.nextDouble();   scan.nextLine();
      genServA = new LognormalGen (new MRG31k3p(), new LognormalDistFromMoments (mu, sigma * sigma));
      // The service time parameters for type B are read in minutes.
      mu = scan.nextDouble() * MINUTE;    sigma = scan.nextDouble() * MINUTE;   scan.nextLine();
      genServB = new LognormalGen (new MRG31k3p(), new LognormalDistFromMoments (mu, sigma * sigma));
      mu = scan.nextDouble();    sigma = scan.nextDouble();   scan.nextLine();
      genDelayB = new NormalGen (new MRG31k3p(), mu, sigma);
      scan.close();
      Locale.setDefault(loc);
   }

   // A customer of type A, created by an ArrivalA event.
   class CustomerA {
      double arrivalTime, serviceTime;

      // Constructor: the customer is created when it arrives.
      public CustomerA() {
         serviceTime = genServA.nextDouble(); // Generate service time.
         if (nTellersBusy < nTellers) { // Start service immediately.
            nTellersBusy++;
            statWaitsInDayA.add(0.0);  // allWaitsA.add(0.0);;
            new ServiceCompletionA().schedule(serviceTime);  
            return;
         }
         // Check if can be served by a counselor.   
         for (int c = 0; c < numCounselors; c++)
            if (tryMatchAD (c, this))  { 
               statWaitsInDayA.add(0.0);  // allWaitsA.add(0.0);
               return;
            }
         // Must join the queue.   
         arrivalTime = Sim.time();
         waitListA.addLast(this);
      }

      // This customer was waiting and now starts its service with a teller.
      // We calculate the wait time and schedule the end of service. 
      public void endWait() {  
         double wait = Sim.time() - arrivalTime;
         statWaitsInDayA.add(wait);   // allWaitsA.add(wait);
         new ServiceCompletionA().schedule(serviceTime);
      }
   }

   // A customer of type B. 
   class CustomerB {
      double arrivalTime, serviceTime;
      int c, slot;  

      // Constructor: this creates the customer which is arriving now.
      // We need to know the counselor and time slot for this customer.
      public CustomerB(int counselor, int slot) {
         // System.out.println("Cust B arrives, current time = " + Sim.time());
         c = counselor;   this.slot = slot;
         serviceTime = genServB.nextDouble(); // Generate service time.
         if (!counselorBusy[c]) { // Start service immediately, no wait.
            // System.out.println("Cust B arrives, no wait, current time = " + Sim.time());
            counselorBusy[c] = true;
            statWaitsInDayB.add(0.0);  // allWaitsB.add(0.0);
            new ServiceCompletionB(c, slot).schedule(serviceTime);
         } else { 
            // Place this customer in wait for counselor.
            // System.out.println("Cust B arrives and must wait, current time = " + Sim.time());
            custBWaiting[c] = this;
            arrivalTime = Sim.time();
         }
      }

      public void endWait() {
         // Same as for type A.
         double wait = Sim.time() - arrivalTime;
         statWaitsInDayB.add(wait);   // allWaitsB.add(wait);
         new ServiceCompletionB(c, slot).schedule(serviceTime);
      }
   }

   // Event: A customer of type A arrives.
   class ArrivalA extends Event {
      public void actions() {
         nextArrivalA.schedule(ExponentialDist.inverseF(arrRateA, streamArrA.nextDouble()));
         nArrivalsA++;
         new CustomerA(); // Creates the customer of type A that just arrived.
      }
   }

   // Event: A customer of type B arrives.
   class ArrivalB extends Event {
      int c, slot;   // We need to know which counselor and time slot.

      public ArrivalB (int counselor, int slot) { 
         c = counselor;  this.slot = slot;
      }
      
      public void actions() {  
         nArrivalsB++;
         new CustomerB (c, slot); // New customer of type B just arrived for c.
      }
   }

   // Event: A customer of type A ends its service with a teller.
   class ServiceCompletionA extends Event {
      public void actions() {
         nTellersBusy--;
         checkQueueA();
      }
   }
 
   // Event: A customer of type A ends its service with a counselor.
   class ServiceCompletionAD extends Event {
      int c; // We need to know which counselor.

      public ServiceCompletionAD (int counselor) { c = counselor; }
      
      public void actions() { 
         if (custBWaiting[c] != null) {
            custBWaiting[c].endWait();
            custBWaiting[c] = null;
         } else {
            counselorBusy[c] = false;
            // Counselor c now free but could have an appointment either now
            // or in next slot.  We much check tryMatchAD for this c.         
            if (waitListA.size() > 0)
               if (tryMatchAD (c, waitListA.getFirst()))
                  (waitListA.removeFirst()).endWait();
         }
      }
   }

   // Event: A customer of type B ends its appointment.
   // We need to know which counselor and what slot,
   // because it can be the customer from the previous slot!
   class ServiceCompletionB extends Event {
      
      int c, slot; // We need to know which counselor and time slot.

      public ServiceCompletionB (int counselor, int slot) {
         c = counselor;  this.slot = slot; 
      }

      public void actions() {
         appoint[c][slot] = false;  // Done.
         if (custBWaiting[c] != null) {
            custBWaiting[c].endWait();
            custBWaiting[c] = null;
         } else {
            counselorBusy[c] = false;
            if (waitListA.size() > 0)
               if (tryMatchAD (c, waitListA.getFirst())) 
                  (waitListA.removeFirst()).endWait();
         }
      }
   }
  
   // Start serving new customers of type A if tellers are free and queue not empty.
   public void checkQueueA() {
      while ((waitListA.size() > 0) && (nTellersBusy < nTellers)) {
         (waitListA.removeFirst()).endWait();
         nTellersBusy++;
      }
   }
  
   // Checks if counselor c is available to serve this customer.
   // If yes, match them, start service, and return true.
   public boolean tryMatchAD (int c, CustomerA cust) {
      // Check first if c is currently free, then check for how long. 
      if ((!appoint[c][currentSlot]) & (!counselorBusy[c])) {
         if ((nextSlotTime - Sim.time() > s) | (!appoint[c][currentSlot+1])) {
            // System.out.println("A counselor starts serving a cust A, current time = " + Sim.time());
            // System.out.println("  nextSlotTime = " + nextSlotTime + ",  s = " + s);
            numMatchesAD++;
            counselorBusy[c] = true;
            new ServiceCompletionAD(c).schedule (cust.serviceTime);             
            return true;  // Found a match!
         }
      }
      return false;
   }

   // Event: A new time period begins.
   class NextPeriod extends Event {
      int j; // Number of the new period.

      public NextPeriod(int period) {
         j = period;
      }

      // This is done when the new period starts. 
      public void actions() {
         if (j < numPeriods) {
            arrRateA = lambda[j];
            nTellers = numTellers[j];
            if (j == 0) {
               nextArrivalA.schedule(ExponentialDist.inverseF(arrRateA, streamArrA.nextDouble()));
            } else {
               checkQueueA();
               nextArrivalA.reschedule((nextArrivalA.time() - Sim.time()) * lambda[j - 1] / lambda[j]);
            }
            new NextPeriod(j + 1).schedule(lengthPeriod);
         } else
            nextArrivalA.cancel(); // End of the day.
      }
   }
   
   // Event: start a new appointment time slot.
   class NextSlotTime extends Event {

      public void actions() {
         currentSlot++;    
         nextSlotTime += slotLength;
         for (int c = 0; c < numCounselors; c++) 
            if (appoint[c][currentSlot])
               if (streamAppoint.nextDouble() < 1-p) {
                  // The customer shows up, we schedule its arrival event.
                  // Here we assume that it cannot start service before the start of its slot!
                  new ArrivalB(c, currentSlot).schedule (Math.max (0.0, genDelayB.nextDouble()));
               }

         if (currentSlot+1 < numSlots)
            new NextSlotTime().schedule(slotLength);
      }
   }
   
   // Determines the presence of appointments for all counselors in all slots.
   public void generateAllAppointments() {
      for (int c = 0; c < numCounselors; c++) {
         counselorBusy[c] = false; // Counselor c is free and no customer is waiting for him.
         custBWaiting[c] = null;   // No customer B is waiting for counselor c.
         for (int i = 0; i < numSlots; i++) 
            appoint[c][i] = (streamAppoint.nextDouble() < r);
         appoint[c][numSlots] = false;  // Extra slot needed for the test in TryMatchAD.
      }
   }
   
   public void simulateOneDay() {
      Sim.init();
      statWaitsInDayA.init();
      statWaitsInDayB.init();
      nArrivalsA = 0;
      nArrivalsB = 0;
      numMatchesAD = 0;
      nTellersBusy = 0;
      currentSlot = -1;  // Will be put to 0 when first slot begins at time 0.
      nextSlotTime = 0.0;
    
      generateAllAppointments();  
      new NextPeriod(0).schedule(0.0);  // Start at time 0.
      new NextSlotTime().schedule(0.0);
      Sim.start();
      // Here the simulation is running......... we exit when simulation stops.

      // Day is over. Add observations for the day to statistical collectors.
      statArrivalsA.add((double) nArrivalsA);
      statArrivalsB.add((double) nArrivalsB);
      statAD.add((double) numMatchesAD);
      statWaitsA.add(statWaitsInDayA.sum() / nExpectedA);
      statWaitsB.add(statWaitsInDayB.sum() / nExpectedB);   
   }

   static public void main(String[] args) throws IOException {
      
      int n = 10000;  // Number of simulation runs.
      BankTwoTypes cc = new BankTwoTypes(
         args.length == 1 ? args[0] : "src/main/docs/examples/ift6561examples/BankTwoTypes.dat");
      for (int i = 0; i < n; i++) {
         // Reset each stream to a new substream.
         // cc.myStreams.resetNextSubstream();
         cc.simulateOneDay();
      }  
      System.out.println("Expected number of type A customers per day = " + cc.nExpectedA + "\n");
      System.out.println("Expected number of type B customers per day = " + cc.nExpectedB + "\n");
      System.out.println (cc.statArrivalsA.report());
      System.out.println (cc.statArrivalsB.report());
      System.out.println (cc.statAD.reportAndCIStudent(0.95, 4));
      System.out.println (cc.statWaitsA.reportAndCIStudent(0.95, 4));
      System.out.println (cc.statWaitsB.reportAndCIStudent(0.95, 4));
   }
}
