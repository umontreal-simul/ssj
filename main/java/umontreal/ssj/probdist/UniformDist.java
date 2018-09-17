/*
 * Class:        UniformDist
 * Description:  uniform distribution over the reals
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author
 * @since
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
 * Extends the class  @ref ContinuousDistribution for the *uniform*
 * distribution @cite tJOH95b&thinsp; (page 276) over the interval
 * @f$[a,b]@f$. Its density is
 * @anchor REF_probdist_UniformDist_eq_funiform
 * @f[
 *   f(x) = 1/(b-a) \qquad\mbox{ for } a\le x\le b \tag{funiform}
 * @f]
 * and 0 elsewhere. The distribution function is
 * @anchor REF_probdist_UniformDist_eq_cdfuniform
 * @f[
 *   F(x) = (x-a)/(b-a) \qquad\mbox{ for } a\le x\le b \tag{cdfuniform}
 * @f]
 * and its inverse is
 * @anchor REF_probdist_UniformDist_eq_cdinvfuniform
 * @f[
 *   F^{-1}(u) = a + (b - a)u \qquad\mbox{for }0 \le u \le1. \tag{cdinvfuniform}
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class UniformDist extends ContinuousDistribution {
   private double a;
   private double b;

   /**
    * Constructs a uniform distribution over the interval @f$(a,b) =
    * (0,1)@f$.
    */
   public UniformDist() {
      setParams (0.0, 1.0);
   }

   /**
    * Constructs a uniform distribution over the interval @f$(a,b)@f$.
    */
   public UniformDist (double a, double b) {
      setParams (a, b);
   }


   public double density (double x) {
      return density (a, b, x);
   }

   public double cdf (double x) {
      return cdf (a, b, x);
   }

   public double barF (double x) {
      return barF (a, b, x);
   }

   public double inverseF (double u) {
      return inverseF (a, b, u);
   }

   public double getMean() {
      return UniformDist.getMean (a, b);
   }

   public double getVariance() {
      return UniformDist.getVariance (a, b);
   }

   public double getStandardDeviation() {
      return UniformDist.getStandardDeviation (a, b);
   }

/**
 * Computes the uniform density function @f$f(x)@f$ in (
 * {@link REF_probdist_UniformDist_eq_funiform funiform} ).
 */
public static double density (double a, double b, double x) {
      if (b <= a)
         throw new IllegalArgumentException ("b <= a");
      if (x <= a || x >= b)
         return 0.0;
      return 1.0 / (b - a);
   }

   /**
    * Computes the uniform distribution function as in (
    * {@link REF_probdist_UniformDist_eq_cdfuniform
    * cdfuniform} ).
    */
   public static double cdf (double a, double b, double x) {
      if (b <= a)
        throw new IllegalArgumentException ("b <= a");
      if (x <= a)
         return 0.0;
      if (x >= b)
         return 1.0;
      return (x - a)/(b - a);
   }

   /**
    * Computes the uniform complementary distribution function
    * @f$\bar{F}(x)@f$.
    */
   public static double barF (double a, double b, double x) {
      if (b <= a)
        throw new IllegalArgumentException ("b <= a");
      if (x <= a)
         return 1.0;
      if (x >= b)
         return 0.0;
      return (b - x)/(b - a);
   }

   /**
    * Computes the inverse of the uniform distribution function (
    * {@link REF_probdist_UniformDist_eq_cdinvfuniform
    * cdinvfuniform} ).
    */
   public static double inverseF (double a, double b, double u) {
       if (b <= a)
           throw new IllegalArgumentException ("b <= a");

       if (u > 1.0 || u < 0.0)
           throw new IllegalArgumentException ("u not in [0, 1]");

       if (u <= 0.0)
           return a;
       if (u >= 1.0)
           return b;
       return a + (b - a)*u;
   }

   /**
    * Estimates the parameter @f$(a, b)@f$ of the uniform distribution
    * using the maximum likelihood method, from the @f$n@f$ observations
    * @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$. The estimates are returned in a
    * two-element array, in regular order: [@f$a@f$, @f$b@f$].  The
    * maximum likelihood estimators are the values @f$(\hat{a}@f$,
    * @f$\hat{b})@f$ that satisfy the equations
    * @f{align*}{
    *    \hat{a} 
    *    & 
    *    = 
    *    \min_i \{x_i\}
    *    \\ 
    *   \hat{b} 
    *    & 
    *    = 
    *    \max_i \{x_i\}.
    * @f}
    * See @cite sLAW00a&thinsp; (page 300).
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{a}@f$, @f$\hat{b}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double parameters[] = new double[2];
      parameters[0] = Double.POSITIVE_INFINITY;
      parameters[1] = Double.NEGATIVE_INFINITY;
      for (int i = 0; i < n; i++) {
         if (x[i] < parameters[0])
            parameters[0] = x[i];
         if (x[i] > parameters[1])
            parameters[1] = x[i];
      }

      return parameters;
   }

   /**
    * Creates a new instance of a uniform distribution with parameters
    * @f$a@f$ and @f$b@f$ estimated using the maximum likelihood method
    * based on the @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1, …,
    * n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static UniformDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new UniformDist (parameters[0], parameters[1]);
   }

   /**
    * Computes and returns the mean @f$E[X] = (a + b)/2@f$ of the uniform
    * distribution with parameters @f$a@f$ and @f$b@f$.
    *  @return the mean of the uniform distribution @f$E[X] = (a + b) /
    * 2@f$
    */
   public static double getMean (double a, double b) {
      if (b <= a)
         throw new IllegalArgumentException ("b <= a");

      return ((a + b) / 2);
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] = (b - a)^2/12@f$
    * of the uniform distribution with parameters @f$a@f$ and @f$b@f$.
    *  @return the variance of the uniform distribution @f$\mbox{Var}[X] =
    * (b - a)^2 / 12@f$
    */
   public static double getVariance (double a, double b) {
      if (b <= a)
         throw new IllegalArgumentException ("b <= a");

      return ((b - a) * (b - a) / 12);
   }

   /**
    * Computes and returns the standard deviation of the uniform
    * distribution with parameters @f$a@f$ and @f$b@f$.
    *  @return the standard deviation of the uniform distribution
    */
   public static double getStandardDeviation (double a, double b) {
      return Math.sqrt (UniformDist.getVariance (a, b));
   }

   /**
    * Returns the parameter @f$a@f$.
    */
   public double getA() {
      return a;
   }

   /**
    * Returns the parameter @f$b@f$.
    */
   public double getB() {
      return b;
   }

   /**
    * Sets the parameters @f$a@f$ and @f$b@f$ for this object.
    */
   public void setParams (double a, double b) {
      if (b <= a)
         throw new IllegalArgumentException ("b <= a");
      this.a = a;
      this.b = b;
      supportA = a;
      supportB = b;
   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$a@f$,
    * @f$b@f$].
    */
   public double[] getParams () {
      double[] retour = {a, b};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : a = " + a + ", b = " + b;
   }

}