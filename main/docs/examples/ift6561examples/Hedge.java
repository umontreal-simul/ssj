package ift6561examples;
import umontreal.ssj.charts.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.stat.*;
import umontreal.ssj.util.*;

public class Hedge {
   static double r;
   static double sigma;
   static double mu;
   double strike;    // Strike price.
   int s;            // Number of observation times.
   double[] zeta;
   double discount;  // Discount factor exp(-r * zeta[t]).
   double[] muDelta; // Differences * (r - sigma^2/2).
   double[] sigmaSqrtDelta; // Square roots of differences * sigma.
   double[] sigmaSqrtTime;  // Square roots of time left * sigma.
   double[] logS;    // Log of the GBM process: logS[t] = log (S[t]).

   // Array zeta[0..s+1] must contain zeta[0]=0.0, plus the s observation times.
   public Hedge (double r, double sigma, double strike, 
                 double s0, int s, double[] zeta) {
      this.r = r;
      this.sigma = sigma;
      this.strike = strike;
      this.s = s;
      this.zeta = zeta;
      discount = Math.exp (-r * zeta[s]);
      mu = r - 0.5 * sigma * sigma;
      muDelta = new double[s];
      sigmaSqrtDelta = new double[s];
      sigmaSqrtTime = new double[s];
      logS = new double[s+1];
      double delta;
      for (int j = 0; j < s; j++) {
         delta = zeta[j+1] - zeta[j];
         muDelta[j] = mu * delta;
         sigmaSqrtDelta[j] = sigma * Math.sqrt (delta);
         sigmaSqrtTime[j] = sigma * Math.sqrt (zeta[s] - zeta[j]);
      }   
      logS[0] = Math.log (s0);
   }

   // Generates the process S. 
   public void generatePath (RandomStream stream) {
       for (int j = 0; j < s; j++)
          logS[j+1] = logS[j] + muDelta[j] + sigmaSqrtDelta[j] 
                   * NormalDist.inverseF01 (stream.nextDouble());
   }

   // Computes and returns the option payoff.
   public double getPayoff () {
      return Math.max (0, (Math.exp (logS[s]) - strike));
   }

   // Returns expected payoff for given remaining time and state s0, 
   // computed by Black-Scholes formula.
   public static double getExpectedPayoff (double r, double sigma, double strike,
                                   double s0, double time) {
      if (time <= 0)
         return Math.max (0, s0 - strike);
      double sigmaSqrtTime = Math.sqrt (time) * sigma;
      double z0 = (Math.log (strike / s0) - mu*time) / sigmaSqrtTime;
      double phi1 = NormalDist.cdf01 (-z0 + sigmaSqrtTime);
      double phi2 = NormalDist.cdf01 (-z0);
      double discount = Math.exp (-r * time);
      return s0 * phi1 - strike * discount * phi2;
   }

   // Computes and returns the derivative of the expected payoff with
   // respect to s0, for given remaining time (the delta of the option).
   public static double getDeltaExpectedPayoff (double r, double sigma, double strike,
                                   double s0, double time) {
      if (time <= 0) {
         if (s0 < strike)
            return 0;
         return 1;
      }
      double sigmaSqrtTime = Math.sqrt (time) * sigma;
      double mdz0 = 1 / (s0 * sigmaSqrtTime);
      double z0 = (Math.log (strike / s0) - mu*time) / sigmaSqrtTime;
      double x = -z0 + sigmaSqrtTime;
      double phi1 = NormalDist.cdf01 (x);
      double dphi1 = NormalDist.density (0, 1, x) * mdz0;
      double dphi2 = NormalDist.density (0, 1, -z0) * mdz0;
      double discount = Math.exp (-r * time);
      return phi1 + dphi1 * s0 - strike * discount * dphi2;
   }

   // Performs delta hedging, and returns the final value of the portfolio.
   public double performDeltaHedging () {
      double sj = Math.exp (logS[0]);
      double v = getExpectedPayoff (r, sigma, strike, sj, zeta[s]);
      double wj = getDeltaExpectedPayoff (r, sigma, strike, sj, zeta[s]);
      double yj = v - wj*sj;
      for (int j = 0; j < s; j++) {
         double sj1 = Math.exp (logS[j + 1]);
         double wj1 = getDeltaExpectedPayoff (r, sigma, strike, sj1, zeta[s] - zeta[j]);
         double yj1 = yj * Math.exp (r*(zeta[j + 1] - zeta[j])) - sj1 * (wj1 - wj);
         wj = wj1;
         yj = yj1;
         sj = sj1;
      }
      return wj * sj + yj;   // Portfolio value at time T = zeta[s].
   }

   // Performs n indep. runs using stream and collects statistics in statValue.
   public void simulateRuns (int n, RandomStream stream) {
      TallyStore statPayoff = new TallyStore ("Payoff", n);
      TallyStore statPortfolio = new TallyStore ("Final value of portfolio", n);
      TallyStore statLoss = new TallyStore ("Net loss (portfolio value - payoff)", n);
      for (int i=0; i<n; i++) {
         generatePath (stream);
         double payoff = getPayoff();
         double value = performDeltaHedging ();
         statPayoff.add (payoff);
         statPortfolio.add (value);
         statLoss.add (value - payoff);
         //  stream.resetNextSubstream();
      }
      // statPayoff.getDoubleArrayList().trimToSize();   //  Do we really need this here ???
      // statPortfolio.getDoubleArrayList().trimToSize();
      System.out.println (statPayoff.reportAndCIStudent (0.95));
      System.out.println (statPortfolio.reportAndCIStudent (0.95));
      System.out.println (statLoss.reportAndCIStudent (0.95));
      System.out.printf ("Expected payoff: %f%n", getExpectedPayoff (r, sigma, strike,
                                                        Math.exp (logS[0]), zeta[s]));

      double[][] data = new double[][] { statPayoff.getArray(), statPortfolio.getArray() };
      new ScatterChart ("Hedging of an option at " + (zeta.length - 1) + " times",
                        "g(s(T))", "p(T)", data).view (800, 600);
      new HistogramChart ("Distribution of the loss, for " + (zeta.length - 1) + " hedging points", 
                          "g(s(T)) - p(T)", "Count", statLoss.getDoubleArrayList()).view (800, 600);
   }

   public static void main (String[] args) { 
      RandomStream stream = new MRG32k3a();
      double r = 0.03;
      double sigma = 0.2;
      double strike = 100;
      double s0 = 100;
      int n = 500;

      for (int s : new int[] { 12, 52, 250, 1000 }) {
         System.out.printf ("\n\n Simulating with %d evaluation times%n", s);
         double[] zeta = new double[s+1];   zeta[0] = 0.0;
         for (int j = 1; j <= s; j++) 
            zeta[j] = (double)j / (double)s;
         Hedge process = new Hedge (r, sigma, strike, s0, s, zeta);
 
         process.simulateRuns (n, stream);
      }
   }
}
