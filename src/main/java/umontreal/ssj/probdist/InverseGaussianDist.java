/*
 * Class:        InverseGaussianDist
 * Description:  inverse Gaussian distribution
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
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.util.*;
import umontreal.ssj.functions.MathFunction;

/**
 * Extends the class  @ref ContinuousDistribution for the *inverse Gaussian*
 * distribution with location parameter @f$\mu> 0@f$ and scale parameter
 * @f$\lambda> 0@f$. Its density is
 * @anchor REF_probdist_InverseGaussianDist_eq_fInverseGaussian
 * @f[
 *   f(x) = \sqrt{\frac{\lambda}{2\pi x^3}}\; e^{{-\lambda(x - \mu)^2}/{(2\mu^2x)}}, \qquad\mbox{for } x > 0. \tag{fInverseGaussian}
 * @f]
 * The distribution function is given by
 * @anchor REF_probdist_InverseGaussianDist_eq_FInverseGaussian
 * @f[
 *   F(x) = \Phi\left(\sqrt{\frac{\lambda}{x}}\left(\frac{x}{\mu} - 1\right)\right) + e^{({2\lambda}/{\mu})}\Phi\left(-\sqrt{\frac{\lambda}{x}}\left(\frac{x}{\mu} + 1\right)\right), \tag{FInverseGaussian}
 * @f]
 * where @f$\Phi@f$ is the standard normal distribution function.
 *
 * The non-static versions of the methods `cdf`, `barF`, and `inverseF` call
 * the static version of the same name.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class InverseGaussianDist extends ContinuousDistribution {
   protected double mu;
   protected double lambda;

   private static class Function implements MathFunction {
      protected double mu;
      protected double lambda;
      protected double u;

      public Function (double mu, double lambda, double u) {
         this.mu = mu;
         this.lambda = lambda;
         this.u = u;
      }

      public double evaluate (double x) {
         return u - cdf(mu, lambda, x);
      }
   }

   /**
    * Constructs the *inverse Gaussian* distribution with parameters
    * @f$\mu@f$ and @f$\lambda@f$.
    */
   public InverseGaussianDist (double mu, double lambda) {
      setParams (mu, lambda);
   }


   public double density (double x) {
      return density (mu, lambda, x);
   }

   public double cdf (double x) {
      return cdf (mu, lambda, x);
   }

   public double barF (double x) {
      return barF (mu, lambda, x);
   }

   public double inverseF (double u) {
      return inverseF (mu, lambda, u);
   }

   public double getMean() {
      return getMean (mu, lambda);
   }

   public double getVariance() {
      return getVariance (mu, lambda);
   }

   public double getStandardDeviation() {
      return getStandardDeviation (mu, lambda);
   }

/**
 * Computes the density function (
 * {@link REF_probdist_InverseGaussianDist_eq_fInverseGaussian
 * fInverseGaussian} ) for the inverse gaussian distribution with parameters
 * @f$\mu@f$ and @f$\lambda@f$, evaluated at @f$x@f$.
 */
public static double density (double mu, double lambda, double x) {
      if (mu <= 0.0)
         throw new IllegalArgumentException ("mu <= 0");
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (x <= 0.0)
         return 0.0;

      double sqrtX = Math.sqrt (x);

      return (Math.sqrt (lambda / (2 * Math.PI)) / (sqrtX * sqrtX * sqrtX) *
              Math.exp (-lambda * (x - 2 * mu + (mu * mu / x)) / (2 * mu * mu)));
   }

   /**
    * Computes the distribution function (
    * {@link REF_probdist_InverseGaussianDist_eq_FInverseGaussian
    * FInverseGaussian} ) of the inverse gaussian distribution with
    * parameters @f$\mu@f$ and @f$\lambda@f$, evaluated at @f$x@f$.
    */
   public static double cdf (double mu, double lambda, double x) {
      if (mu <= 0.0)
         throw new IllegalArgumentException ("mu <= 0");
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (x <= 0.0)
         return 0.0;
      double temp = Math.sqrt (lambda / x);
      double z = temp * (x / mu - 1.0);
      double w = temp * (x / mu + 1.0);

      // C'est bien un + dans    exp (2 * lambda / mu)
      return (NormalDist.cdf01 (z) +
              Math.exp (2 * lambda / mu) * NormalDist.cdf01 (-w));
   }

   /**
    * Computes the complementary distribution function of the inverse
    * gaussian distribution with parameters @f$\mu@f$ and @f$\lambda@f$,
    * evaluated at @f$x@f$.
    */
   public static double barF (double mu, double lambda, double x) {
      return 1.0 - cdf (mu, lambda, x);
   }

   /**
    * Computes the inverse of the inverse gaussian distribution with
    * parameters @f$\mu@f$ and @f$\lambda@f$.
    */
   public static double inverseF (double mu, double lambda, double u) {
      if (mu <= 0.0)
         throw new IllegalArgumentException ("mu <= 0");
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u must be in [0,1]");
      if (u == 1.0)
         return Double.POSITIVE_INFINITY;
      if (u == 0.0)
         return 0.0;

      Function f = new Function (mu, lambda, u);

      // Find interval containing root = x*
      double sig = getStandardDeviation(mu, lambda);
      double x0 = 0.0;
      double x = mu;
      double v = cdf(mu, lambda, x);
      while (v < u) {
         x0 = x;
         x += 3.0*sig;
         v = cdf(mu, lambda, x);
      }

      return RootFinder.brentDekker (x0, x, f, 1e-12);
   }

   /**
    * Estimates the parameters @f$(\mu, \lambda)@f$ of the inverse
    * gaussian distribution using the maximum likelihood method, from the
    * @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1,…, n-1@f$. The
    * estimates are returned in a two-element array, in regular order:
    * [@f$\mu@f$, @f$\lambda@f$].  The maximum likelihood estimators are
    * the values @f$(\hat{\mu}, \hat{\lambda})@f$ that satisfy the
    * equations:
    * @f{align*}{
    *    \hat{\mu} 
    *    & 
    *    = 
    *    \bar{x}_n 
    *    \\ 
    *    \frac{1}{\hat{\lambda}} 
    *    & 
    *    = 
    *    \frac{1}{n} \sum_{i=1}^n \left(\frac{1}{x_i} - \frac{1}{\hat{\mu}}\right),
    * @f}
    * where @f$\bar{x}_n@f$ is the average of @f$x[0],…,x[n-1]@f$,
    * @cite tJOH95a&thinsp; (page 271).
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{\mu}@f$,
    * @f$\hat{\lambda}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double parameters[];
      parameters = new double[2];
      double sum = 0;
      for (int i = 0; i < n; i++) {
         sum += x[i];   
      }
      parameters[0] = sum / (double) n;

      sum = 0;
      for (int i = 0; i < n; i++) {
         sum += ((1.0 / (double) x[i]) - (1.0 / parameters[0]));
      }
      parameters[1] = (double) n / (double) sum;

      return parameters;
   }

   /**
    * Creates a new instance of an inverse gaussian distribution with
    * parameters @f$\mu@f$ and @f$\lambda@f$ estimated using the maximum
    * likelihood method based on the @f$n@f$ observations @f$x[i]@f$, @f$i
    * = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static InverseGaussianDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new InverseGaussianDist (parameters[0], parameters[1]);
   }

   /**
    * Returns the mean @f$E[X] = \mu@f$ of the inverse gaussian
    * distribution with parameters @f$\mu@f$ and @f$\lambda@f$.
    *  @return the mean of the inverse gaussian distribution @f$E[X] =
    * \mu@f$
    */
   public static double getMean (double mu, double lambda) {
      if (mu <= 0.0)
         throw new IllegalArgumentException ("mu <= 0");
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return mu;      
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] =
    * \mu^3/\lambda@f$ of the inverse gaussian distribution with
    * parameters @f$\mu@f$ and @f$\lambda@f$.
    *  @return the variance of the inverse gaussian distribution
    * @f$\mbox{Var}[X] = \mu^3 / \lambda@f$
    */
   public static double getVariance (double mu, double lambda) {
      if (mu <= 0.0)
         throw new IllegalArgumentException ("mu <= 0");
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return (mu * mu * mu / lambda);
   }

   /**
    * Computes and returns the standard deviation of the inverse gaussian
    * distribution with parameters @f$\mu@f$ and @f$\lambda@f$.
    *  @return the standard deviation of the inverse gaussian distribution
    */
   public static double getStandardDeviation (double mu, double lambda) {
      return Math.sqrt (getVariance (mu, lambda));
   }

   /**
    * Returns the parameter @f$\lambda@f$ of this object.
    */
   public double getLambda() {
      return lambda;
   }

   /**
    * Returns the parameter @f$\mu@f$ of this object.
    */
   public double getMu() {
      return mu;
   }

   /**
    * Sets the parameters @f$\mu@f$ and @f$\lambda@f$ of this object.
    */
   public void setParams (double mu, double lambda) {
      if (mu <= 0.0)
         throw new IllegalArgumentException ("mu <= 0");
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      this.mu = mu;
      this.lambda = lambda;
      supportA = 0.0;
   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$\mu@f$,
    * @f$\lambda@f$].
    */
   public double[] getParams () {
      double[] retour = {mu, lambda};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : mu = " + mu + ", lambda = " + lambda;
   }

}