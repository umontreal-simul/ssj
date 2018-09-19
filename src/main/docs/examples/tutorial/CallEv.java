import umontreal.ssj.simevents.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.stat.Tally;
import java.io.*;
import java.util.*;

public class CallEv {

   static final double HOUR = 3600.0;  // Time is in seconds.

   // Data
   // Arrival rates are per hour, service and patience times are in seconds.
   int numDays;           // Number of days to simulate.
   double openingTime;    // Opening time of the center (in hours).
   int numPeriods;        // Number of working periods (hours) in the day.
   int[] numAgents;       // Number of agents for each period.
   double[] lambda;       // Base arrival rate lambda_j for each j.
   double alpha0;         // Parameter of gamma distribution for W.
   double p;              // Probability that patience time is 0.
   double nu;             // Parameter of exponential for patience time.
   double alpha, beta;    // Parameters of gamma service time distribution.
   double s;              // Want stats on waiting times smaller than s.

   // Variables
   double busyness;       // Current value of W.
   double arrRate = 0.0;  // Current arrival rate.
   int nAgents;           // Number of agents in current period.
   int nBusy;             // Number of agents occupied;
   int nArrivals;         // Number of arrivals today;
   int nAbandon;          // Number of abandonments during the day.
   int nGoodQoS;          // Number of waiting times less than s today.
   double nCallsExpected; // Expected number of calls per day.

   Event nextArrival = new Arrival();           // The next Arrival event.

   RandomStream streamW        = new MRG32k3a(); // For W.
   RandomStream streamArr      = new MRG32k3a(); // For arrivals.
   RandomStream streamPatience = new MRG32k3a(); // For patience times.
   GammaGen genServ;   // For service times; created in readData().

   LinkedListStat<Call> waitList = new LinkedListStat<Call> ("Waiting calls");

   Tally statArrivals = new Tally ("Number of arrivals per day");
   Tally statWaits    = new Tally ("Average waiting time per customer");
   Tally statWaitsDay = new Tally ("Waiting times within a day");
   Tally statGoodQoS  = new Tally ("Proportion of waiting times < s");
   Tally statAbandon  = new Tally ("Proportion of calls lost");

   class Call { double arrivTime, servTime, patienceTime; }

   static public void main (String[] args) throws IOException {
      new CallEv(args.length == 1 ? args[0] : "CallEv.dat");
   }

   public CallEv(String fileName) throws IOException {
      readData(fileName);
      for (int i=1; i <= numDays; i++)  simulOneDay();
      System.out.println ("\n Num. calls expected = " + nCallsExpected + "\n");
      statArrivals.setConfidenceIntervalStudent();
      statWaits.setConfidenceIntervalStudent();
      statGoodQoS.setConfidenceIntervalStudent();
      statAbandon.setConfidenceIntervalStudent();
      System.out.println (statArrivals.report (0.9, 3));
      System.out.println (statWaits.report (0.9, 3));
      System.out.println (statGoodQoS.report (0.9, 3));
      System.out.println (statAbandon.report (0.9, 3));
   }

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
            new NextPeriod(j + 1).schedule (1.0 * HOUR);
         } else
            nextArrival.cancel();  // End of the day.
      }
   }

   class Arrival extends Event {
      public void actions() {
         nextArrival.schedule
            (ExponentialDist.inverseF (arrRate, streamArr.nextDouble()));
         nArrivals++;
         Call call = new Call();               // Call just arrived.
         call.servTime = genServ.nextDouble(); // Generate service time.
         if (nBusy < nAgents) {          // Start service immediately.
            nBusy++;
            nGoodQoS++;
            statWaitsDay.add (0.0);
            new CallCompletion().schedule (call.servTime);
         } else {                        // Join the queue.
            call.patienceTime = generPatience();
            call.arrivTime = Sim.time();
            waitList.addLast (call);
         }
      }
   }

   class CallCompletion extends Event {
      public void actions() { nBusy--;   checkQueue(); }
   }

   public void checkQueue() {
      // Start answering new calls if agents are free and queue not empty.
      while ((waitList.size() > 0) && (nBusy < nAgents)) {
         Call call = waitList.removeFirst();
         double wait = Sim.time() - call.arrivTime;
         if (call.patienceTime < wait) { // Caller has abandoned.
            nAbandon++;
            wait = call.patienceTime;    // Effective waiting time.
         } else {
            nBusy++;
            new CallCompletion().schedule (call.servTime);
         }
         if (wait < s) nGoodQoS++;
         statWaitsDay.add (wait);
      }
   }

   public double generPatience() {
      // Generates the patience time for a call.
      double u = streamPatience.nextDouble();
      if (u <= p)
         return 0.0;
      else
         return ExponentialDist.inverseF (nu, (1.0 - u) / (1.0 - p));
   }

   public void readData(String fileName) throws IOException {
      // Reads data and construct arrays.
      Locale loc = Locale.getDefault();
      Locale.setDefault(Locale.US); // to read reals as 8.3 instead of 8,3
      BufferedReader input = new BufferedReader (new FileReader (fileName));
      Scanner scan = new Scanner(input);
      numDays = scan.nextInt();
      scan.nextLine();
      openingTime = scan.nextDouble();
      scan.nextLine();
      numPeriods = scan.nextInt();
      scan.nextLine();
      numAgents = new int[numPeriods];
      lambda = new double[numPeriods];
      nCallsExpected = 0.0;
      for (int j = 0; j < numPeriods; j++) {
         numAgents[j] = scan.nextInt();
         lambda[j] = scan.nextDouble();
         nCallsExpected += lambda[j];
         scan.nextLine();
      }
      alpha0 = scan.nextDouble();
      scan.nextLine();
      p = scan.nextDouble();
      scan.nextLine();
      nu = scan.nextDouble();
      scan.nextLine();
      alpha = scan.nextDouble();
      scan.nextLine();
      beta = scan.nextDouble();
      scan.nextLine();
      s = scan.nextDouble();
      scan.close();
      Locale.setDefault(loc);

      // genServ can be created only after its parameters are known.
      genServ = new GammaAcceptanceRejectionGen ( // Faster than inversion
                   new MRG32k3a(), alpha, beta);
   }

   public void simulOneDay() {
      Sim.init();        statWaitsDay.init();
      nArrivals = 0;     nAbandon = 0;
      nGoodQoS = 0;      nBusy = 0;
      busyness = GammaDist.inverseF (alpha0, alpha0, 8, streamW.nextDouble());
      new NextPeriod(0).schedule (openingTime * HOUR);
      Sim.start();
      // Here the simulation is running...
      statArrivals.add ((double)nArrivals);
      statAbandon.add ((double)nAbandon / nCallsExpected);
      statGoodQoS.add ((double)nGoodQoS / nCallsExpected);
      statWaits.add (statWaitsDay.sum() / nCallsExpected);
   }
}
