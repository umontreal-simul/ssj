package tutorial;
import umontreal.ssj.stat.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.ExponentialDist;
import umontreal.ssj.util.Chrono;

public class QueueLindley {

   RandomStream streamArr  = new MRG32k3a();
   RandomStream streamServ = new MRG32k3a();
   Tally averageWaits = new Tally ("Average waits");
 
   public double simulate (int numCust, double lambda, double mu) {
      double Wi = 0.0;
      double sumWi = 0.0;
      for (int i = 2; i <= numCust; i++) {
         Wi += ExponentialDist.inverseF (mu, streamServ.nextDouble()) -
               ExponentialDist.inverseF (lambda, streamArr.nextDouble());
         if (Wi < 0.0) Wi = 0.0;
         sumWi += Wi;
      }
      return sumWi / numCust;
   }

   public void simulateRuns (int n, int numCust, double lambda, double mu) {
      averageWaits.init();
      for (int i=0; i<n; i++)
	  averageWaits.add (simulate (numCust, lambda, mu));
   }

   public static void main (String[] args) { 
      // Chrono timer = new Chrono();
      QueueLindley queue = new QueueLindley();
      queue.simulateRuns (100, 10000, 1.0, 2.0);
      System.out.println (queue.averageWaits.report());
      // System.out.println ("Total CPU time: " + timer.format());
   }
}
