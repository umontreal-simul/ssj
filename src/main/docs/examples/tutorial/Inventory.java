package tutorial;
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.*;

public class Inventory {

   double lambda;  // Mean demand size for a day.
   double c;       // Sale price per item.
   double h;       // Inventory cost per item per day.
   double K;       // Fixed ordering cost per order.
   double k;       // Marginal ordering cost per item.
   double p;       // Probability that an order arrives.

   RandomVariateGenInt genDemand;
   RandomStream streamDemand = new MRG32k3a();
   RandomStream streamOrder  = new MRG32k3a();
   Tally statProfit          = new Tally ("stats on profit");

   public Inventory (double lambda, double c, double h,
                     double K, double k, double p) {
      this.lambda = lambda;
      this.c = c;  this.h = h;  this.K = K;  this.k = k;  this.p = p;
      genDemand = new PoissonGen (streamDemand, new PoissonDist (lambda));
   }

   // Simulates the system for m days, with the (s,S) policy,
   // and returns the average profit per day.
   public double simulate (int m, int s, int S) {
      int Xj = S;         // Stock in the morning.
      int Yj;             // Stock in the evening.
      double profit = 0.0;    // Cumulated profit.
      for (int j = 0; j < m; j++) {
         Yj = Xj - genDemand.nextInt(); // Subtract demand for the day.
         if (Yj < 0) Yj = 0;            // Lost demand.
         profit += c * (Xj - Yj) - h * Yj;
         if ((Yj < s) && (streamOrder.nextDouble() < p)) {
            // We have a successful order.
            profit -= K + k * (S - Yj);
            Xj = S;
         } else
            Xj = Yj;
      }
      return profit / m;
   }

   // Performs n independent simulation runs of the system for m days,
   // with the (s,S) policy, and returns a report with a 90% confidence
   // interval on the expected average profit per day.
   public void simulateRuns (int n, int m, int s, int S) {
      for (int i = 0; i < n; i++)
         statProfit.add (simulate (m, s, S));
   }

   public static void main (String[] args) {
      Inventory system = new Inventory (100.0, 2.0, 0.1, 10.0, 1.0, 0.95);
	  Chrono timer = new Chrono();
      system.simulateRuns (500, 2000, 80, 200);
      system.statProfit.setConfidenceIntervalStudent();
      System.out.println (system.statProfit.report (0.9, 3));
      System.out.println ("Total CPU time: " + timer.format() + "\n");
   }
}
