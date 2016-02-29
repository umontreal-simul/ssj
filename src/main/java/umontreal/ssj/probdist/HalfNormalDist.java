/*
 * Class:        HalfNormalDist
 * Description:  half-normal distribution
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
package  umontreal.ssj.probdist;
import umontreal.ssj.util.*;
import umontreal.ssj.functions.MathFunction;
import optimization.*;

/**
 * Extends the class  @ref ContinuousDistribution for the *half-normal*
 * distribution with parameters @f$\mu@f$ and @f$\sigma> 0@f$. Its density
 * is
 * @anchor REF_probdist_HalfNormalDist_eq_fHalfNormal
 * @f{align}{
 *    f(x) 
 *    & 
 *   =
 *    \frac{1}{\sigma}\sqrt{\frac{2}{\pi}}\; e^{-(x-\mu)^2/2\sigma^2}, \qquad\mbox{for } x \ge\mu. 
 *    \\ 
 *    \tag{fHalfNormal} f(x) 
 *    & 
 *   =
 *    0, \qquad\mbox{for } x < \mu. \nonumber
 * @f}
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class HalfNormalDist extends ContinuousDistribution {
   protected double mu;
   protected double sigma;
   protected double C1;

   /**
    * Constructs a `HalfNormalDist` object with parameters @f$\mu=@f$
    * `mu` and @f$\sigma=@f$ `sigma`.
    */
   public HalfNormalDist (double mu, double sigma) {
      setParams (mu, sigma);
   }


   public double density (double x) {
      final double z = (x-mu)/sigma;
      if (z < 0.0)
         return 0.0;
      return C1 * Math.exp(-z*z/2.0);
   }

   public double cdf (double x) {
      return cdf (mu, sigma, x);
   }

   public double barF (double x) {
      return barF (mu, sigma, x);
   }

   public double inverseF (double u) {
      return inverseF (mu, sigma, u);
   }

   public double getMean() {
      return HalfNormalDist.getMean (mu, sigma);
   }

   public double getVariance() {
      return HalfNormalDist.getVariance (mu, sigma);
   }

   public double getStandardDeviation() {
      return HalfNormalDist.getStandardDeviation (mu, sigma);
   }

/**
 * Computes the density function of the *half-normal* distribution.
 *  @param mu           the parameter mu
 *  @param sigma        the parameter sigma
 *  @param x            the value at which the density is evaluated
 *  @return returns the density function
 */
public static double density (double mu, double sigma, double x) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      final double Z = (x-mu)/sigma;
      if (Z < 0.0) return 0.0;
      return Math.sqrt(2.0/Math.PI) / sigma * Math.exp(-Z*Z/2.0);
   }

   /**
    * Computes the distribution function.
    *  @param mu           the parameter mu
    *  @param sigma        the parameter sigma
    *  @param x            the value at which the distribution is
    *                      evaluated
    *  @return returns the cdf function
    */
   public static double cdf (double mu, double sigma, double x) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      final double Z = (x-mu)/sigma;
      if (Z <= 0.0) return 0.0;
      return Num.erf(Z/Num.RAC2);
   }

   /**
    * Computes the complementary distribution function.
    *  @param mu           the parameter mu
    *  @param sigma        the parameter sigma
    *  @param x            the value at which the complementary
    *                      distribution is evaluated
    *  @return returns the complementary distribution function
    */
   public static double barF (double mu, double sigma, double x) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      final double Z = (x-mu)/sigma;
      if (Z <= 0.0) return 1.0;
      return Num.erfc(Z/Num.RAC2);
   }

   /**
    * Computes the inverse of the distribution function.
    *  @param mu           the parameter mu
    *  @param sigma        the parameter sigma
    *  @param u            the value at which the inverse distribution is
    *                      evaluated
    *  @return returns the inverse distribution function
    */
   public static double inverseF (double mu, double sigma, double u) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      if (u > 1.0 || u < 0.0)
         throw new IllegalArgumentException ("u not in [0,1]");
      if (u <= 0.0) return mu;
      if (u >= 1.0)
         return Double.POSITIVE_INFINITY;

      final double Z = Num.RAC2 * Num.erfInv(u);
      return mu + sigma * Z;
   }

   /**
    * Estimates the parameters @f$\mu@f$ and @f$\sigma@f$ of the
    * half-normal distribution using the maximum likelihood method from
    * the @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$. The
    * estimates are returned in a two-element array: [@f$\mu@f$,
    * @f$\sigma@f$].  The maximum likelihood estimators are the values
    * @f$\hat{\mu}@f$ and @f$\hat{\sigma}@f$ that satisfy the equation
    * @f{align*}{
    *    \hat{\mu}= \min_j \{x_j\}, 
    *    \\ 
    *    \hat{\sigma}= \sqrt{\frac{1}{n}\Sigma_j(x_j-\hat{\mu})^2}.
    * @f}
    * @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\mu@f$, @f$\sigma@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double mu = Double.MAX_VALUE;
      for (int i = 0 ; i < n ; ++i)
         if (x[i] < mu)
            mu = x[i];

      double sigma = 0.0;
      for (int i = 0 ; i < n ; ++i)
         sigma += (x[i]-mu)*(x[i]-mu);

      double[] parametres = new double [2];
      parametres[0] = mu;
      parametres[1] = Math.sqrt(sigma/n);
      return parametres;
   }

   /**
    * Estimates the parameter @f$\sigma@f$ of the half-normal
    * distribution using the maximum likelihood method from the @f$n@f$
    * observations @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$ and the parameter
    * @f$\mu@f$ = `mu`. The estimate is returned in a one-element array:
    * [@f$\sigma@f$].  The maximum likelihood estimator is the value
    * @f$\hat{\sigma}@f$ that satisfies the equation
    * @f{align*}{
    *    \hat{\sigma}= \sqrt{\frac{1}{n}\Sigma_j(x_j-\mu)^2}.
    * @f}
    * @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameter
    *  @param mu           the parameter mu
    *  @return returns the parameter [@f$\sigma@f$]
    */
   public static double[] getMLE (double[] x, int n, double mu) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double sigma = 0.0;
      for (int i = 0 ; i < n ; ++i)
         sigma += (x[i]-mu)*(x[i]-mu);

      double[] parametres = new double [1];
      parametres[0] = Math.sqrt(sigma/n);
      return parametres;
   }

   /**
    * Computes and returns the mean @f$ E[X] = \mu+ \sigma\sqrt{2 /
    * \pi}. @f$
    *  @param mu           the parameter mu
    *  @param sigma        the parameter sigma
    *  @return returns the mean
    */
   public static double getMean (double mu, double sigma) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      return mu + sigma*Math.sqrt(2.0/Math.PI);
   }

   /**
    * Computes and returns the variance @f$ \mbox{Var}[X] =
    * \left(1-2/\pi\right)\sigma^2. @f$
    *  @param mu           the parameter mu
    *  @param sigma        the parameter sigma
    *  @return returns the variance
    */
   public static double getVariance (double mu, double sigma) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      return (1.0 - 2.0/Math.PI)*sigma*sigma;
   }

   /**
    * Computes the standard deviation of the half-normal distribution with
    * parameters @f$\mu@f$ and @f$\sigma@f$.
    *  @param mu           the parameter mu
    *  @param sigma        the parameter sigma
    *  @return returns the standard deviation
    */
   public static double getStandardDeviation (double mu, double sigma) {
      return Math.sqrt (HalfNormalDist.getVariance (mu, sigma));
   }

   /**
    * Returns the parameter @f$\mu@f$ of this object.
    *  @return returns the parameter mu
    */
   public double getMu() {
      return mu;
   }

   /**
    * Returns the parameter @f$\sigma@f$ of this object.
    *  @return returns the parameter sigma
    */
   public double getSigma() {
      return sigma;
   }

   /**
    * Sets the parameters @f$\mu@f$ and @f$\sigma@f$.
    *  @param mu           the parameter mu
    *  @param sigma        the parameter sigma
    */
   public void setParams (double mu, double sigma) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      this.mu = mu;
      this.sigma = sigma;
      C1 = Math.sqrt(2.0/Math.PI) / sigma;
    }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$\mu@f$,
    * @f$\sigma@f$].
    *  @return returns the parameters [@f$\mu@f$, @f$\sigma@f$]
    */
   public double[] getParams () {
      double[] retour = {mu, sigma};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    *  @return returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : mu = " + mu + ", sigma = " + sigma;
   }

}