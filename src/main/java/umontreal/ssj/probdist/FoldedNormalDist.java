/*
 * Class:        FoldedNormalDist
 * Description:  folded normal distribution
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
 * Extends the class  @ref ContinuousDistribution for the *folded normal*
 * distribution with parameters @f$\mu\ge0@f$ and @f$\sigma> 0@f$. The
 * density is
 * @anchor REF_probdist_FoldedNormalDist_eq_fFoldedNormal
 * @f[
 *   f(x) = \phi\left(\frac{x-\mu}{\sigma}\right) + \phi\left(\frac{-x-\mu}{\sigma}\right) \qquad\mbox{for } x \ge0, \tag{fFoldedNormal}
 * @f]
 * @f[
 *   f(x) = 0, \qquad\mbox{ for } x < 0,
 * @f]
 * where @f$ \phi@f$ denotes the density function of a standard normal
 * distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class FoldedNormalDist extends ContinuousDistribution {
   protected double mu;
   protected double sigma;
   private static final double RACPI = 1.7724538509055160273; // Sqrt[PI]

   private static class FunctionInverse implements MathFunction {
        private double u, mu, sigma;

        public FunctionInverse (double mu, double sigma, double u) {
            this.u = u;
            this.mu = mu;
            this.sigma = sigma;
        }

        public double evaluate (double x) {
            return u - cdf(mu, sigma, x);
        }
    }

   /**
    * Constructs a `FoldedNormalDist` object with parameters @f$\mu=@f$
    * `mu` and @f$\sigma=@f$ `sigma`.
    */
   public FoldedNormalDist (double mu, double sigma) {
      setParams (mu, sigma);
   }


   public double density (double x) {
      return density (mu, sigma, x);
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
      return FoldedNormalDist.getMean (mu, sigma);
   }

   public double getVariance() {
      return FoldedNormalDist.getVariance (mu, sigma);
   }

   public double getStandardDeviation() {
      return FoldedNormalDist.getStandardDeviation (mu, sigma);
   }

/**
 * Computes the density function of the *folded normal* distribution.
 *  @param mu           the parameter mu
 *  @param sigma        the parameter sigma
 *  @param x            the value at which the density is evaluated
 *  @return returns the density function
 */
public static double density (double mu, double sigma, double x) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      if (mu < 0.0)
         throw new IllegalArgumentException ("mu < 0");
      if (x < 0.0) return 0.0;
      return NormalDist.density(mu,sigma,x) + NormalDist.density(mu,sigma,-x);
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
      if (mu < 0.0)
         throw new IllegalArgumentException ("mu < 0");
      if (x <= 0.0) return 0.0;
      return NormalDist.cdf01((x-mu)/sigma) - NormalDist.cdf01((-x-mu)/sigma);
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
      if (mu < 0.0)
         throw new IllegalArgumentException ("mu < 0");
      if (x <= 0.0) return 1.0;
      return NormalDist.barF01((x-mu)/sigma) - NormalDist.barF01((-x-mu)/sigma);
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
      if (mu < 0.0)
         throw new IllegalArgumentException ("mu < 0");
      if (u > 1.0 || u < 0.0)
         throw new IllegalArgumentException ("u not in [0,1]");
      if (u <= 0.0) return 0.0;
      if (u >= 1.0)
         return Double.POSITIVE_INFINITY;

      MathFunction f = new FunctionInverse (mu, sigma, u);
      return RootFinder.brentDekker (0.0, mu + 10.0*sigma, f, 1.0e-14);
   }

   /**
    * Computes and returns the mean
    * @f[
    *   E[X] = \sigma\sqrt{\frac{2}{\pi}}\; e^{-\mu^2/(2\sigma^2)} + \mu \mbox{erf}\left(\frac{\mu}{\sigma\sqrt{2}}\right),
    * @f]
    * where erf@f$(z)@f$ is the error function.
    *  @param mu           the parameter mu
    *  @param sigma        the parameter sigma
    *  @return returns the mean
    */
   public static double getMean (double mu, double sigma) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      if (mu < 0.0)
         throw new IllegalArgumentException ("mu < 0");

      return sigma * Num.RAC2 / RACPI * Math.exp(-mu*mu/(2.0*sigma*sigma))
             + mu * Num.erf(mu/(sigma*Num.RAC2));
   }

   /**
    * Computes and returns the variance
    * @f[
    *   \mbox{Var}[X] = \mu^2 + \sigma^2 - E[X]^2.
    * @f]
    * @param mu           the parameter mu
    *  @param sigma        the parameter sigma
    *  @return returns the variance
    */
   public static double getVariance (double mu, double sigma) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      if (mu < 0.0)
         throw new IllegalArgumentException ("mu < 0");
      double mean = sigma * Num.RAC2 / RACPI * Math.exp(-mu*mu/(2.0*sigma*sigma))
                    + mu * Num.erf(mu/(sigma*Num.RAC2));
      return mu*mu + sigma*sigma - mean*mean;
   }

   /**
    * Computes the standard deviation of the folded normal distribution
    * with parameters @f$\mu@f$ and @f$\sigma@f$.
    *  @param mu           the parameter mu
    *  @param sigma        the parameter sigma
    *  @return returns the standard deviation
    */
   public static double getStandardDeviation (double mu, double sigma) {
      return Math.sqrt (FoldedNormalDist.getVariance (mu, sigma));
   }

   /**
    * NOT IMPLEMENTED. Les formules pour le MLE sont données dans
    * @cite tLEO61a&thinsp;.
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{\mu}@f$,
    * @f$\hat{\sigma}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      throw new UnsupportedOperationException("getMLE is not implemented ");
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
    * Sets the parameters @f$\mu@f$ and @f$\sigma@f$ for this object.
    *  @param mu           the parameter mu
    *  @param sigma        the parameter sigma
    */
   public void setParams (double mu, double sigma) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      if (mu < 0.0)
         throw new IllegalArgumentException ("mu < 0");
      this.mu = mu;
      this.sigma = sigma;
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