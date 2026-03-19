package ift6561examples;
// package umontreal.ssj.finance;

import umontreal.ssj.stat.Tally;
import umontreal.ssj.stochprocess.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.mcqmctools.*;

/**
 * This is the same as `AsianOption`, but with a slightly modified payoff function.
 * Here the payoff is based on the average price at the observation times, 
 * except that the prices at time 0 and at the final time count for 1/2
 * compared with the other observations. It is also assumed that the observation times
 * are equally spaced, stating at 0.  This payoff is to approximate a continuous
 * average, expressed as an integral of the price over the time interval [0,T].
 * 
 */
public class AsianCallOptionPLin extends AsianOption {

   /**
    * Constructor. The d+1 observation times are equally spaced, from 0 to T.
    */
   public AsianCallOptionPLin (double r, int d, double T, double strike) {
      super (r, d, 0, T, strike);
   }

   /**
    * Computes and returns discounted payoff. Assumes path has been generated.
    */
   public double getPerformance() {
      double average = 0.5 * path[0]; // Average over sample path.
      for (int j = 1; j < d; j++)
         average += path[j];
      average += 0.5 * path[d];
      average /= d;
      if (average > strike)
         return discount * (average - strike);
      else
         return 0.0;
   }

   public String toString() {
      return "Asian option model with " + d + " observation times, piecewise linear payoff";
   }

}
