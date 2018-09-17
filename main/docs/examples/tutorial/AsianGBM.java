package tutorial;

import java.io.IOException;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.*;

public class AsianGBM {
   double strike;    // Strike price.
   int d;            // Number of observation times.
   double discount;  // Discount factor, exp(-r * zeta[d]).
   double[] muDelta; // muDelta[j] = (zeta[j+1] - zeta[j]) * (r - sigma^2/2).
   double[] sigmaSqrtDelta; // sqrt(zeta[j+1] - zeta[j]) * sigma.
   double[] logS;    // Log of the GBM process: logS[t] = log (S[t]).

   // Array zeta[0..s] must contain zeta[0]=0.0, plus the d observation times.
   // This constructor precomputes several quantities to speedup the simulation.
   public AsianGBM (double r, double sigma, double strike,
                    double s0, int d, double[] zeta) {
      this.strike = strike;
      this.d = d;
      discount = Math.exp (-r * zeta[d]);
      double mu = r - 0.5 * sigma * sigma;
      muDelta = new double[d];
      sigmaSqrtDelta = new double[d];
      logS = new double[d+1];
      double delta;
      for (int j = 0; j < d; j++) {
         delta = zeta[j+1] - zeta[j];
         muDelta[j] = mu * delta;
         sigmaSqrtDelta[j] = sigma * Math.sqrt (delta);
      }
      logS[0] = Math.log (s0);
   }

   // Generates the log of the process S.
   public void generatePath (RandomStream stream) {
       for (int j = 0; j < d; j++)
          logS[j+1] = logS[j] + muDelta[j] + sigmaSqrtDelta[j]
                   * NormalDist.inverseF01 (stream.nextDouble());
   }

   // Computes and returns the discounted option payoff.
   public double getPayoff () {
       double average = 0.0;  // Average of the GBM process.
       for (int j = 1; j <= d; j++) average += Math.exp (logS[j]);
       average /= d;
       if (average > strike) return discount * (average - strike);
       else return 0.0;
   }

   // Performs n simulation runs using stream and collects statistics in statValue.
   public void simulateRuns (int n, RandomStream stream, Tally statValue) {
      statValue.init();
      for (int i=0; i<n; i++) {
         generatePath (stream);
         statValue.add (getPayoff ());
         stream.resetNextSubstream();
      }
   }

   public static void main (String[] args)  throws IOException {
      int d = 12;
      double[] zeta = new double[d+1];   zeta[0] = 0.0;
      for (int j=1; j<=d; j++)
         zeta[j] = (double)j / (double)d;
      AsianGBM process = new AsianGBM (0.05, 0.5, 100.0, 100.0, d, zeta);
      Tally statValue = new Tally ("Stats on value of Asian option");

      Chrono timer = new Chrono();
      int n = 1000000;
      process.simulateRuns (n, new MRG32k3a(), statValue);
      statValue.setConfidenceIntervalStudent();
      System.out.println (statValue.report (0.95, 3));
      System.out.println ("Total CPU time:      " + timer.format() + "\n");
   }
}
