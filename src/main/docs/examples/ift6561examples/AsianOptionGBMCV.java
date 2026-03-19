package ift6561examples;

import umontreal.ssj.stochprocess.*;
import umontreal.ssj.mcqmctools.*;
import umontreal.ssj.probdist.NormalDist;

/**
 * This class represents an Asian Option based on a geometric Brownian motion,
 * for which we use the payoff under a geometric average (instead of the
 * arithmetic average) as a control variate. It is a subclass of `AsianOption`.
 */

public class AsianOptionGBMCV extends AsianOption implements MonteCarloModelCV {

   double expectedGeo;  // Expected value of first CV (payoff under geometric average).
   double expectedSumS; // Expected value of second CV (sum of values of S(t_j)).

   /**
    * Array <TT>obsTimes[0..d+1]</TT> must contain <TT>obsTimes[0] = 0</TT>, plus
    * the <SPAN CLASS="MATH"><I>d</I></SPAN> observation times.
    * 
    */
   public AsianOptionGBMCV(double r, int d, double[] obsTimes, double strike) {
      super(r, d, obsTimes, strike);
   }

   public AsianOptionGBMCV(GeometricBrownianMotion sp, double r, int d, double[] obsTimes, double strike) {
      super(r, d, obsTimes, strike);
      setProcess(sp); // Must make sure we compute expectedGeo and expectedSumS.
   }

   public AsianOptionGBMCV(double r, int d, double T1, double T, double strike) {
      super(r, d, T1, T, strike);
   }

   public AsianOptionGBMCV(GeometricBrownianMotion sp, double r, int d, double T1, double T, double strike) {
      super(r, d, T1, T, strike);
      setProcess(sp); // Must make sure we compute expectedGeo and expectedSumS.
   }

   /**
    * Sets (or resets) the process to `sp', which much be a GBM in this class, and
    * resets the observation times of `sp` to those of this `AsianOption` object.
    * Also computes the expected values of the CVs.  
    */
   public void setProcess(GeometricBrownianMotion sp) {
      // Reset the process to sp. Assumes that obsTimes have been set.
      priceProcess = sp;
      sp.setObservationTimes(obsTimes, d);
      computeExpectedGeo();
      computeExpectedSumS();
   }

   /**
    * Computes the expected value of geometric CV.  This is done only once.
    */
   public double computeExpectedGeo() {
      // First, get the underlying BM and its parameters.
      BrownianMotion bm = ((GeometricBrownianMotion) priceProcess).getBrownianMotion();
      double s0 = priceProcess.getX0(); // Initial value of the GBM.
      double mu = bm.getMu();
      double sigma = bm.getSigma();
      double my = 0;
      double s2y = 0;
      for (int j = 1; j <= d; j++) {
         my += obsTimes[j];
         s2y += (obsTimes[j] - obsTimes[j - 1]) * (d - j + 1) * (d - j + 1);
      }
      my = Math.log(s0) + my * mu / d;
      s2y *= sigma * sigma / (d * d);
      double dd = (-Math.log(strike) + my) / Math.sqrt(s2y);
      return expectedGeo = discount
            * (Math.exp(my + 0.5 * s2y) * NormalDist.cdf01(dd + Math.sqrt(s2y)) 
               - strike * NormalDist.cdf01(dd));
   }

   /**
    * Computes the expected value of second CV, the sum of S(t_j). Done only once.
    */
   public double computeExpectedSumS() {
      double s0 = priceProcess.getX0(); // Initial value of the GBM.
      double sum = 0.0;
      for (int j = 1; j <= d; j++)
         sum += Math.exp (r * obsTimes[j]);
      return expectedSumS = sum * s0;
   }
   
   /**
    * Returns the expected value of geometric CV (must have been computed before).
    * 
    */
   public double getExpectedGeo() {
      return expectedGeo;
   }
   
   /**
    * Returns the expected value of second CV (the sum, must have been computed before).
    * 
    */
   public double getExpectedSumS() {
      return expectedSumS;
   }

   /**
    * Computes and returns discounted payoff for geometric average. Assumes that
    * the path has been generated.
    * 
    */
   public double getPayoffGeo() {
      double[] pathBM = ((GeometricBrownianMotion) priceProcess).getBrownianMotion().getPath();
      double average = 0.0; // Average over BM sample path.
      for (int j = 1; j <= d; j++)
         average += pathBM[j];
      average /= d;
      average = path[0] * Math.exp(average);
      if (average > strike)
         return discount * (average - strike);
      else
         return 0.0;
   }

   /**
    * Computes and returns the value of the second CV, the sum of S(t_j). Assumes that
    * the path has been generated.
    * 
    */
   public double getSumS() {
      double sum = 0.0; 
      for (int j = 1; j <= d; j++)
          sum += path[j];
      return sum;
   }

   // Recovers the realizations of the control variates from the last run.
   // Returns the geometric average minus its expectation. Assumes that sp is a GBM.
   // This is useful if we want to use `MonteCarloExperiment.simulateRunsCV()`. 
   public double[] getValuesCV() {
      double[] cv = new double[2];
      cv[0] = (getPayoffGeo() - expectedGeo);
      cv[1] = (getSumS() - expectedSumS);
      return cv;
   }

   // Returns the number of control variates.
   public int getNumberCV() {
      return 2;
   }

   // Returns the realization of the first control variate, centered.
   public double getCV1() {
      return (getPayoffGeo() - expectedGeo);
   }
   
   // Returns the realization of the second control variate, centered.
   public double getCV2() {
      return (getSumS() - expectedSumS);
   }
   
   public String toString() {
      return "Asian option model with " + d + " observation times,"
            + "under a GBM process,\n with payoff for geometric average as a control variate.";
   }
}
