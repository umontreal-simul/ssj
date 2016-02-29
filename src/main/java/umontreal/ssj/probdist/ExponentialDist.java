/*
 * Class:        ExponentialDist
 * Description:  exponential distribution
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
import  umontreal.ssj.util.Num;

/**
 * Extends the class  @ref ContinuousDistribution for the *exponential*
 * distribution @cite tJOH95a&thinsp; (page 494) with mean @f$1/\lambda@f$
 * where @f$\lambda> 0@f$. Its density is
 * @anchor REF_probdist_ExponentialDist_eq_fexpon
 * @f[
 *   f(x) = \lambda e^{-\lambda x} \qquad\mbox{for }x\ge0, \tag{fexpon}
 * @f]
 * its distribution function is
 * @anchor REF_probdist_ExponentialDist_eq_Fexpon
 * @f[
 *   F(x) = 1 - e^{-\lambda x},\qquad\mbox{for }x \ge0, \tag{Fexpon}
 * @f]
 * and its inverse distribution function is
 * @f[
 *   F^{-1}(u) = -\ln(1-u)/\lambda, \qquad\mbox{for } 0 < u < 1.
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class ExponentialDist extends ContinuousDistribution {
   private double lambda;

   /**
    * Constructs an `ExponentialDist` object with parameter @f$\lambda@f$
    * = 1.
    */
   public ExponentialDist() {
      setLambda (1.0);
   }

   /**
    * Constructs an `ExponentialDist` object with parameter @f$\lambda@f$
    * = `lambda`.
    */
   public ExponentialDist (double lambda) {
      setLambda (lambda);
  }


   public double density (double x) {
      return density (lambda, x);
   }

   public double cdf (double x) {
      return cdf (lambda, x);
   }

   public double barF (double x) {
      return barF (lambda, x);
   }

   public double inverseF (double u) {
      return inverseF (lambda, u);
   }

   public double getMean() {
      return ExponentialDist.getMean (lambda);
   }

   public double getVariance() {
      return ExponentialDist.getVariance (lambda);
   }

   public double getStandardDeviation() {
      return ExponentialDist.getStandardDeviation (lambda);
   }

/**
 * Computes the density function.
 */
public static double density (double lambda, double x) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      return x < 0 ? 0 : lambda*Math.exp (-lambda*x);
   }

   /**
    * Computes the distribution function.
    */
   public static double cdf (double lambda, double x) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (x <= 0.0)
         return 0.0;
      double y = lambda * x;
      if (y >= XBIG)
         return 1.0;
      return -Math.expm1 (-y);
   }

   /**
    * Computes the complementary distribution function.
    */
   public static double barF (double lambda, double x) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (x <= 0.0)
         return 1.0;
      if (lambda*x >= XBIGM)
         return 0.0;
         return Math.exp (-lambda*x);
   }

   /**
    * Computes the inverse distribution function.
    */
   public static double inverseF (double lambda, double u) {
        if (lambda <= 0)
           throw new IllegalArgumentException ("lambda <= 0");
        if (u < 0.0 || u > 1.0)
            throw new IllegalArgumentException ("u not in [0,1]");
        if (u >= 1.0)
            return Double.POSITIVE_INFINITY;
        if (u <= 0.0)
            return 0.0;
        return -Math.log1p (-u)/lambda;
   }

   /**
    * Estimates the parameter @f$\lambda@f$ of the exponential
    * distribution using the maximum likelihood method, from the @f$n@f$
    * observations @f$x[i]@f$, @f$i = 0, 1,…, n-1@f$. The estimate is
    * returned in a one-element array, as element 0.  The equation of the
    * maximum likelihood is defined as @f$\hat{\lambda} = 1/\bar{x}_n@f$,
    * where @f$\bar{x}_n@f$ is the average of @f$x[0],…,x[n-1]@f$ (see
    * @cite tJOH95a&thinsp; (page 506)).
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameter [@f$\hat{\lambda}@f$]
    */
   public static double[] getMLE (double[] x, int n)
   {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double parameters[];
      double sum = 0.0;
      parameters = new double[1];
      for (int i = 0; i < n; i++)
         sum+= x[i];
      parameters[0] = (double) n / sum;

      return parameters;
   }

   /**
    * Creates a new instance of an exponential distribution with parameter
    * @f$\lambda@f$ estimated using the maximum likelihood method based
    * on the @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static ExponentialDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new ExponentialDist (parameters[0]);
   }

   /**
    * Computes and returns the mean, @f$E[X] = 1/\lambda@f$, of the
    * exponential distribution with parameter @f$\lambda@f$.
    *  @return the mean of the exponential distribution @f$E[X] = 1 /
    * \lambda@f$
    */
   public static double getMean (double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return (1 / lambda);
   }

   /**
    * Computes and returns the variance, @f$\mbox{Var}[X] =
    * 1/\lambda^2@f$, of the exponential distribution with parameter
    * @f$\lambda@f$.
    *  @return the variance of the Exponential distribution
    * @f$\mbox{Var}[X] = 1 / \lambda^2@f$
    */
   public static double getVariance (double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return (1 / (lambda * lambda));
   }

   /**
    * Computes and returns the standard deviation of the exponential
    * distribution with parameter @f$\lambda@f$.
    *  @return the standard deviation of the exponential distribution
    */
   public static double getStandardDeviation (double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return (1 / lambda);
   }

   /**
    * Returns the value of @f$\lambda@f$ for this object.
    */
   public double getLambda() {
      return lambda;
   }

   /**
    * Sets the value of @f$\lambda@f$ for this object.
    */
   public void setLambda (double lambda) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      this.lambda = lambda;
      supportA = 0.0;
   }

   /**
    * Return a table containing the parameters of the current
    * distribution.
    */
   public double[] getParams () {
      double[] retour = {lambda};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : lambda = " + lambda;
   }

}