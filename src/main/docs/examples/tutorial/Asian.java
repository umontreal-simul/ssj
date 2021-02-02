package tutorial;

import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.*;

public class Asian {
   double strike;    // Strike price.
   int s;            // Number of observation times.
   double discount;  // Discount factor exp(-r * zeta[t]).
   double[] muDelta; // Differences * (r - sigma^2/2).
   double[] sigmaSqrtDelta; // Square roots of differences * sigma.
   double[] logS;    // Log of the GBM process: logS[t] = log (S[t]).

   // Array zeta[0..s] must contain zeta[0]=0.0, plus the s observation times.
   public Asian (double r, double sigma, double strike,
                 double s0, int s, double[] zeta) {
      this.strike = strike;
      this.s = s;
      discount = Math.exp (-r * zeta[s]);
      double mu = r - 0.5 * sigma * sigma;
      muDelta = new double[s];
      sigmaSqrtDelta = new double[s];
      logS = new double[s+1];
      double delta;
      for (int j = 0; j < s; j++) {
         delta = zeta[j+1] - zeta[j];
         muDelta[j] = mu * delta;
         sigmaSqrtDelta[j] = sigma * Math.sqrt (delta);
      }
      logS[0] = Math.log (s0);
   }

   // Generates the process S.
   public void generatePath (RandomStream stream) {
       for (int j = 0; j < s; j++)
          logS[j+1] = logS[j] + muDelta[j] + sigmaSqrtDelta[j]
                   * NormalDist.inverseF01 (stream.nextDouble());
   }

   // Computes and returns the discounted option payoff.
   public double getPayoff () {
       double average = 0.0;  // Average of the GBM process.
       for (int j = 1; j <= s; j++) average += Math.exp (logS[j]);
       average /= s;
       if (average > strike) return discount * (average - strike);
       else return 0.0;
   }

   // Performs n indep. runs using stream and collects statistics in statValue.
   public void simulateRuns (int n, RandomStream stream, Tally statValue) {
      statValue.init();
      for (int i=0; i<n; i++) {
         generatePath (stream);
         statValue.add (getPayoff ());
         stream.resetNextSubstream();
      }
   }

   public static void main (String[] args) {
      int s = 12;
      double[] zeta = new double[s+1];   zeta[0] = 0.0;
      for (int j=1; j<=s; j++)
         zeta[j] = (double)j / (double)s;
      Asian process = new Asian (0.05, 0.5, 100.0, 100.0, s, zeta);
      Tally statValue = new Tally ("Stats on value of Asian option");

      Chrono timer = new Chrono();
      int n = 100000;
      process.simulateRuns (n, new MRG32k3a(), statValue);
      statValue.setConfidenceIntervalStudent();
      System.out.println (statValue.report (0.95, 3));
      System.out.println ("Total CPU time:      " + timer.format() + "\n");
   }
}
