/*
 * Class:        BernoulliDist
 * Description:  Bernoulli distribution
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        August 2010
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package umontreal.ssj.probdist;

/**
 * Extends the class  @ref DiscreteDistributionInt for the *Bernoulli*
 * distribution @cite sLAW00a&thinsp; with parameter @f$p@f$, where
 * @f$0\le p\le1@f$. Its mass function is given by
 * @anchor REF_probdist_BernoulliDist_eq_fmass_bernoulli
 * @f[
 *   f(x) = \begin{cases}
 *    1 - p, \qquad
 *    & 
 *    \text{if $x = 0$;} 
 *    \\ 
 *   p, 
 *    & 
 *    \text{if $x = 1$;} \tag{fmass-bernoulli} 
 *    \\ 
 *   0, 
 *    & 
 *    \text{otherwise.} 
 *   \end{cases}
 * @f]
 * Its distribution function is
 * @anchor REF_probdist_BernoulliDist_eq_cdf_bernoulli
 * @f[
 *   F(x) = \begin{cases}
 *    0, 
 *    & 
 *    \text{if $x < 0$;} 
 *    \\ 
 *   1 - p, \qquad
 *    & 
 *    \text{if $0 \le x < 1$;} \tag{cdf-bernoulli} 
 *    \\ 
 *   1, 
 *    & 
 *    \text{if $x \ge1$.} 
 *   \end{cases}
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_discrete
 */
public class BernoulliDist extends DiscreteDistributionInt {
   private double p;
   private double q;

   /**
    * Creates a Bernoulli distribution object.
    */
   public BernoulliDist (double p) {
      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in [0,1]");
      this.p = p;
      q = 1 - p;
      supportA = 0;
      supportB = 1;
   }

   public double prob (int x) {
      if (1 == x) return p;
      if (0 == x) return q;
      return 0.0;
   }

   public double cdf (int x) {
      if (x < 0) return 0.0;
      if (x < 1) return q;
      return 1.0;
   }

   public double barF (int x) {
      if (x > 1) return 0.0;
      if (x > 0) return p;
      return 1.0;
   }

   public int inverseFInt (double u) {
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u not in [0,1]");
      if (u > q) return 1;
      return 0;
   }

   public double getMean() {
      return BernoulliDist.getMean (p);
   }

   public double getVariance() {
      return BernoulliDist.getVariance (p);
   }

   public double getStandardDeviation() {
      return BernoulliDist.getStandardDeviation (p);
   }

/**
 * Returns the Bernoulli probability @f$f(x)@f$ with parameter @f$p@f$ (see
 * eq. ( {@link REF_probdist_BernoulliDist_eq_fmass_bernoulli
 * fmass-bernoulli} )).
 */
public static double prob (double p, int x) {
      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in [0,1]");
      if (1 == x) return p;
      if (0 == x) return 1.0 - p;
      return 0.0;
   }

   /**
    * Returns the Bernoulli distribution function @f$F(x)@f$ with
    * parameter @f$p@f$ (see eq. (
    * {@link REF_probdist_BernoulliDist_eq_cdf_bernoulli
    * cdf-bernoulli} )).
    */
   public static double cdf (double p, int x) {
      if (p < 0.0 | p > 1.0)
         throw new IllegalArgumentException ("p not in [0,1]");
      if (x < 0) return 0.0;
      if (x < 1) return 1.0 - p;
      return 1.0;
   }

   /**
    * Returns the complementary Bernoulli distribution function
    * @f$\bar{F}(x) = P[X \ge x]@f$ with parameter @f$p@f$.
    */
   public static double barF (double p, int x) {
      if (p < 0.0 | p > 1.0)
         throw new IllegalArgumentException ("p not in [0,1]");
      if (x > 1) return 0.0;
      if (x > 0) return p;
      return 1.0;
   }

   /**
    * Returns the inverse of the Bernoulli distribution function with
    * parameter @f$p@f$ at @f$u@f$.
    */
   public static int inverseF (double p, double u) {
      if (p < 0.0 | p > 1.0)
         throw new IllegalArgumentException ("p not in [0,1]");
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u not in [0,1]");
      if (u > 1.0 - p) return 1;
      return 0;
    }

   /**
    * Estimates the parameters @f$p@f$ of the Bernoulli distribution using
    * the maximum likelihood method, from the @f$m@f$ observations
    * @f$x[i]@f$, @f$i = 0, 1,…, m-1@f$. The estimate is returned in a
    * one-element array: [@f$p@f$].
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param m            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameter [@f$\hat{p}@f$]
    */
   public static double[] getMLE (int[] x, int m) {
      if (m < 2)
         throw new UnsupportedOperationException(" m < 2");

      // compute the empirical mean
      double sum = 0.0;
      for (int i = 0; i < m; i++)
         sum += x[i];
      sum /= (double) m;

      double param[] = new double[1];
      param[0] = sum;
      return param;
   }

   /**
    * Creates a new instance of a Bernoulli distribution with parameter
    * @f$p@f$ estimated using the maximum likelihood method, from the
    * @f$m@f$ observations @f$x[i]@f$, @f$i = 0, 1, …, m-1@f$.
    *  @param x            the list of observations to use to estimate the
    *                      parameters
    *  @param m            the number of observations to use to estimate
    *                      the parameters
    */
   public static BernoulliDist getInstanceFromMLE (int[] x, int m) {
      double param[] = new double[1];
      param = getMLE (x, m);
      return new BernoulliDist (param[0]);
   }

   /**
    * Returns the mean @f$E[X] = p@f$ of the Bernoulli distribution with
    * parameter @f$p@f$.
    *  @return the mean of the Bernoulli distribution @f$E[X] = np@f$
    */
   public static double getMean (double p) {
      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in [0, 1]");

      return (p);
   }

   /**
    * Computes the variance @f$\mbox{Var}[X] = p(1 - p)@f$ of the
    * Bernoulli distribution with parameter @f$p@f$.
    *  @return the variance of the Bernoulli distribution
    */
   public static double getVariance (double p) {
      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in [0, 1]");
      return (p * (1.0 - p));
   }

   /**
    * Computes the standard deviation of the Bernoulli distribution with
    * parameter @f$p@f$.
    *  @return the standard deviation of the Bernoulli distribution
    */
   public static double getStandardDeviation (double p) {
      return Math.sqrt (BernoulliDist.getVariance (p));
   }

   /**
    * Returns the parameter @f$p@f$ of this object.
    */
   public double getP() {
      return p;
   }

   /**
    * Returns an array that contains the parameter @f$p@f$ of the current
    * distribution: [@f$p@f$].
    */
   public double[] getParams () {
      double[] retour = {p};
      return retour;
   }

   /**
    * Resets the parameter to this new value.
    */
   public void setParams (double p) {
      this.p = p;
      q = 1 - p;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : p = " + p;
   }

}