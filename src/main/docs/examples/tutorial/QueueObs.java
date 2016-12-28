package tutorial;
import java.util.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.simevents.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;

public class QueueObs {

   Tally waitingTimes = new Tally ("Waiting times");
   Tally averageWaits = new Tally ("Average wait");
   RandomVariateGen genArr;
   RandomVariateGen genServ;
   int cust;    // Number of the current customer.

   public QueueObs (double lambda, double mu, int step) {
      genArr = new ExponentialGen (new MRG32k3a(), lambda);
      genServ = new ExponentialGen (new MRG32k3a(), mu);
      waitingTimes.setBroadcasting (true);
      waitingTimes.addObservationListener (new ObservationTrace (step));
      waitingTimes.addObservationListener (new LargeWaitsCollector (2.0));
   }

   public double simulateOneRun (int numCust) {
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
	  averageWaits.add (simulateOneRun (numCust));
   }

   public class ObservationTrace implements ObservationListener {
      private int step;

      public ObservationTrace (int step) { this.step = step; }

      public void newObservation (StatProbe probe, double x) {
         if (cust % step == 0)
            System.out.println ("Customer " + cust + " waited " 
                   + x + " time units.");
      }
   }

   public class LargeWaitsCollector implements ObservationListener {
      double threshold;
      ArrayList<Double> largeWaits = new ArrayList<Double>();

      public LargeWaitsCollector (double threshold) {
         this.threshold = threshold;
      }

      public void newObservation (StatProbe probe, double x) {
         if (x > threshold) largeWaits.add (x);
      }

      public String formatLargeWaits () {
	  // Should print the list largeWaits. 
	  return "not yet implemented...";
      }
   }

   public static void main (String[] args) { 
      QueueObs queue = new QueueObs (1.0, 2.0, 5);
      queue.simulateRuns (2, 100);
      System.out.println ("\n\n" + queue.averageWaits.report());
   }
}
