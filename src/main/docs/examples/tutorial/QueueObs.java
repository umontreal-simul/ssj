package tutorial;
import java.util.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;

public class QueueObs {

   Tally waitingTimes = new Tally ("Waiting times");
   Tally averageWaits = new Tally ("Average wait");
   RandomVariateGen genArr;  // For interarrival times.
   RandomVariateGen genServ; // For service times.
   int cust;    // Number of the current customer.

   public QueueObs (double lambda, double mu, int step) {
      genArr = new ExponentialGen (new MRG32k3a(), lambda);
      genServ = new ExponentialGen (new MRG32k3a(), mu);
      waitingTimes.setBroadcasting (true);
      waitingTimes.addObservationListener (new ObservationTrace (step));
      waitingTimes.addObservationListener (new LargeWaitsCollector (2.0));
   }

   public double simulate (int numCust) {
      waitingTimes.init();
      double Wi = 0.0;
      waitingTimes.add (Wi);
      for (cust = 2; cust <= numCust; cust++) {
         Wi += genServ.nextDouble() - genArr.nextDouble();
         if (Wi < 0.0) Wi = 0.0;
         waitingTimes.add (Wi);
      }
      return waitingTimes.average();
   }

   public void simulateRuns (int n, int numCust) {
      averageWaits.init();
      for (int i=0; i<n; i++)
	     averageWaits.add (simulate (numCust));
   }

   // A listener that observes each waiting time and prints every `step`th one.
   public class ObservationTrace implements ObservationListener {
      private int step;

      public ObservationTrace (int step) { this.step = step; }

      public void newObservation (StatProbe probe, double x) {
         if (cust % step == 0)
            System.out.println ("Customer " + cust + " waited " 
                   + x + " time units.");
      }
   }

   // A listener that observes waiting times and collects those larger than threshold.
   public class LargeWaitsCollector implements ObservationListener {
      double threshold;
      ArrayList<Double> largeWaits = new ArrayList<Double>();

      public LargeWaitsCollector (double threshold) {
         this.threshold = threshold;
      }

      public void newObservation (StatProbe probe, double x) {
         if (x > threshold) largeWaits.add (x);
      }

	  // Maybe print the list largeWaits. 
      public String formatLargeWaits () {
	     return "not yet implemented...";
      }
   }

   public static void main (String[] args) { 
      QueueObs queue = new QueueObs (1.0, 2.0, 5);
      queue.simulateRuns (2, 100);
      System.out.println ("\n\n" + queue.averageWaits.report());
   }
}
