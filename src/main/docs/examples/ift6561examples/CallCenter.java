package ift6561examples;
import umontreal.ssj.simevents.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.stat.Tally;
import java.io.*;
import java.util.StringTokenizer;
import java.util.LinkedList;

public class CallCenter {

   static final double HOUR = 3600.0;  // Time is in seconds. 

   // Data
   // Arrival rates are per hour, service and patience times are in seconds.
   double openingTime;    // Opening time of the center (in hours).
   int numPeriods;        // Number of working periods (hours) in the day.
   int[] numAgents;       // Number of agents for each period.   
   double[] lambda;       // Base arrival rate lambda_j for each j.
   double alpha0;         // Parameter of gamma distribution for B.
   double p;              // Probability that patience time is 0.
   double nu;             // Parameter of exponential for patience time.
   double alpha, gamma;   // Parameters of gamma service time distribution. 
   double s;              // Want stats on waiting times smaller than s.

   // Variables
   double busyness;       // Current value of B.
   double arrRate = 0.0;  // Current arrival rate.
   int nAgents;           // Number of agents in current period.
   int nBusy;             // Number of agents occupied;
   int nArrivals;         // Number of arrivals today;
   int nAbandon;          // Number of abandonments during the day.
   int nGoodQoS;          // Number of waiting times less than s today.
   double nCallsExpected; // Expected number of calls per day.

   Event nextArrival = new Arrival();           // The next Arrival event.
   LinkedList<Call> waitList = new LinkedList<Call> ();

   RandomStream streamB        = new MRG32k3a(); // For B.
   RandomStream streamArr      = new MRG32k3a(); // For arrivals.
   RandomStream streamPatience = new MRG32k3a(); // For patience times.
   GammaGen genServ;   // For service times; created in readData().
 
   Tally statArrivals = new Tally ("Number of arrivals per day");
   Tally statWaits    = new Tally ("Average waiting time per customer");
   Tally statWaitsDay = new Tally ("Waiting times within a day");
   Tally statGoodQoS  = new Tally ("Proportion of waiting times < s");
   Tally statAbandon  = new Tally ("Proportion of calls lost");

   public CallCenter (String fileName) throws IOException {
      readData (fileName);
      // genServ can be created only after its parameters are read.
      // The acceptanc/rejection method is much faster than inversion.
      genServ = new GammaAcceptanceRejectionGen (new MRG32k3a(),
                    new GammaDist (alpha, gamma));
   }

   // Reads data and construct arrays.
   public void readData (String fileName) throws IOException {
      BufferedReader input = new BufferedReader (new FileReader (fileName));
      StringTokenizer line = new StringTokenizer (input.readLine());
      openingTime = Double.parseDouble (line.nextToken());
      line = new StringTokenizer (input.readLine());
      numPeriods  = Integer.parseInt (line.nextToken());

      numAgents = new int[numPeriods];
      lambda = new double[numPeriods];
      nCallsExpected = 0.0;
      for (int j=0; j < numPeriods; j++) {
         line = new StringTokenizer (input.readLine());
         numAgents[j] = Integer.parseInt (line.nextToken());
         lambda[j]    = Double.parseDouble (line.nextToken());
         nCallsExpected += lambda[j];
      }
      line = new StringTokenizer (input.readLine());
      alpha0 = Double.parseDouble (line.nextToken());
      line = new StringTokenizer (input.readLine());
      p = Double.parseDouble (line.nextToken());
      line = new StringTokenizer (input.readLine());
      nu = Double.parseDouble (line.nextToken());
      line = new StringTokenizer (input.readLine());
      alpha = Double.parseDouble (line.nextToken());
      line = new StringTokenizer (input.readLine());
      gamma = Double.parseDouble (line.nextToken());
      line = new StringTokenizer (input.readLine());
      s = Double.parseDouble (line.nextToken());
      input.close();
   }

   // A phone call.
   class Call { 

      double arrivalTime, serviceTime, patienceTime;

      public Call() {
         serviceTime = genServ.nextDouble(); // Generate service time.
         if (nBusy < nAgents) {           // Start service immediately.
            nBusy++;
            nGoodQoS++;
            statWaitsDay.add (0.0);
            new CallCompletion().schedule (serviceTime);
         } else {                         // Join the queue.
            patienceTime = generPatience();
            arrivalTime = Sim.time();
            waitList.addLast (this);
         }
      }
      public void endWait() {
         double wait = Sim.time() - arrivalTime;
         if (patienceTime < wait) { // Caller has abandoned.
            nAbandon++;
            wait = patienceTime;    // Effective waiting time.
         }
         else {
            nBusy++;
            new CallCompletion().schedule (serviceTime);
         }
         if (wait < s) nGoodQoS++;
         statWaitsDay.add (wait);
      }
   } 

   // Event: A new period begins.
   class NextPeriod extends Event {
      int j;     // Number of the new period.
      public NextPeriod (int period) { j = period; }
      public void actions() {
         if (j < numPeriods) {
            nAgents = numAgents[j];
            arrRate = busyness * lambda[j] / HOUR;
            if (j == 0) 
               nextArrival.schedule 
                  (ExponentialDist.inverseF (arrRate, streamArr.nextDouble()));
            else {
               checkQueue();
               nextArrival.reschedule ((nextArrival.time() - Sim.time()) 
                                       * lambda[j-1] / lambda[j]);
            }
            new NextPeriod(j+1).schedule (1.0 * HOUR);
         }
         else
            nextArrival.cancel();  // End of the day.
      }
   }

   // Event: A call arrives.
   class Arrival extends Event {
      public void actions() {
         nextArrival.schedule 
            (ExponentialDist.inverseF (arrRate, streamArr.nextDouble()));
         nArrivals++;
         new Call();               // Call just arrived.
      }
   }

   // Event: A call is completed.
   class CallCompletion extends Event {
      public void actions() { nBusy--;   checkQueue(); }
   }

   // Start answering new calls if agents are free and queue not empty.
   public void checkQueue() {
      while ((waitList.size() > 0) && (nBusy < nAgents))
         ((Call)waitList.removeFirst()).endWait();
   }

   // Generates the patience time for a call.
   public double generPatience() {
      double u = streamPatience.nextDouble();
      if (u <= p) 
         return 0.0;
      else 
         return ExponentialDist.inverseF (nu, (1.0-u) / (1.0-p));
   }

   public void simulateOneDay (double busyness) { 
      Sim.init();        statWaitsDay.init();
      nArrivals = 0;     nAbandon = 0;     
      nGoodQoS = 0;      nBusy = 0;
      this.busyness = busyness;

      new NextPeriod(0).schedule (openingTime * HOUR);
      Sim.start();
      // Here the simulation is running...

      statArrivals.add ((double)nArrivals);
      statAbandon.add ((double)nAbandon / nCallsExpected);
      statGoodQoS.add ((double)nGoodQoS / nCallsExpected);
      statWaits.add (statWaitsDay.sum() / nCallsExpected);
   }

   public void simulateOneDay () {
      simulateOneDay (GammaDist.inverseF (alpha0, alpha0, 8, 
                                       streamB.nextDouble()));
   }

   static public void main (String[] args) throws IOException { 
      CallCenter cc = new CallCenter ("CallCenter.dat"); 
      for (int i=1; i <= 1000; i++)  cc.simulateOneDay();
      System.out.println ("\nNum. calls expected = " + cc.nCallsExpected +"\n");
      System.out.println (
         cc.statArrivals.reportAndCIStudent (0.9) +
         cc.statWaits.reportAndCIStudent (0.9) +
         cc.statGoodQoS.reportAndCIStudent (0.9) +
         cc.statAbandon.reportAndCIStudent (0.9)); 
   }
}
