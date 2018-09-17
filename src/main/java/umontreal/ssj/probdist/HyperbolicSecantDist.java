/*
 * Class:        HyperbolicSecantDist
 * Description:  hyperbolic secant distribution
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
import optimization.*;

/**
 * Extends the class  @ref ContinuousDistribution for the *hyperbolic secant*
 * distribution with location parameter @f$\mu@f$ and scale parameter
 * @f$\sigma> 0@f$. Its density is
 * @anchor REF_probdist_HyperbolicSecantDist_eq_fHyperbolicSecant
 * @f[
 *   f(x) = \frac{1}{2 \sigma} \mbox{ sech}\left(\frac{\pi}{2} \frac{(x - \mu)}{\sigma}\right) \tag{fHyperbolicSecant}
 * @f]
 * The distribution function is given by
 * @anchor REF_probdist_HyperbolicSecantDist_eq_FHyperbolicSecant
 * @f[
 *   F(x) = \frac{2}{\pi} \tan^{-1}\left[\exp{\left(\frac{\pi}{2} \frac{(x - \mu)}{\sigma}\right)}\right] \tag{FHyperbolicSecant}
 * @f]
 * The non-static versions of the methods `cdf`, `barF`, and `inverseF` call
 * the static version of the same name.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class HyperbolicSecantDist extends ContinuousDistribution {
   protected double mu;
   protected double sigma;
   private static final double ZLIMB = 500.0;
   private static final double ZLIMS = 50.0;

   private static class Optim implements Uncmin_methods
   {
      private int n;
      private double[] xi;

      public Optim (double[] x, int n)
      {
         this.n = n;
         this.xi = new double[n];
         System.arraycopy (x, 0, this.xi, 0, n);
      }

      public double f_to_minimize (double[] p)
      {
         double sum = 0.0;

         if (p[2] <= 0.0)
            return 1e200;

         for (int i = 0; i < n; i++)
            sum -= Math.log (density (p[1], p[2], xi[i]));

         return sum;
      }

      public void gradient (double[] x, double[] g)
      {
      }

      public void hessian (double[] x, double[][] h)
      {
      }
   }

   /**
    * Constructs a hyperbolic secant distribution with parameters
    * @f$\mu@f$ and @f$\sigma@f$.
    */
   public HyperbolicSecantDist (double mu, double sigma) {
      setParams (mu, sigma);
   }


   public double density (double x) {
      return HyperbolicSecantDist.density (mu, sigma, x);
   }

   public double cdf (double x) {
      return HyperbolicSecantDist.cdf (mu, sigma, x);
   }

   public double barF (double x) {
      return HyperbolicSecantDist.barF (mu, sigma, x);
   }

   public double inverseF (double u) {
      return HyperbolicSecantDist.inverseF (mu, sigma, u);
   }

   public double getMean() {
      return HyperbolicSecantDist.getMean (mu, sigma);
   }

   public double getVariance() {
      return HyperbolicSecantDist.getVariance (mu, sigma);
   }

   public double getStandardDeviation() {
      return HyperbolicSecantDist.getStandardDeviation (mu, sigma);
   }

/**
 * Computes the density function (
 * {@link REF_probdist_HyperbolicSecantDist_eq_fHyperbolicSecant
 * fHyperbolicSecant} ) for a hyperbolic secant distribution with parameters
 * @f$\mu@f$ and @f$\sigma@f$.
 */
public static double density (double mu, double sigma, double x) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      double y = (x - mu) / sigma;
      if (Math.abs(y) >= ZLIMB)
         return 0.0;
      else
         return (1.0 / (Math.cosh (Math.PI * y / 2.0) * 2.0 * sigma));
   }

   /**
    * Computes the distribution function of the hyperbolic secant
    * distribution with parameters @f$\mu@f$ and @f$\sigma@f$.
    */
   public static double cdf (double mu, double sigma, double x) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      double y = (x - mu) / sigma;
      if (y >= ZLIMS)
         return 1.0;
      else if (y <= -ZLIMB)
      	 return 0.0;
      else
         return (2.0 * Math.atan (Math.exp (Math.PI * y / 2.0))) / Math.PI;
   }

   /**
    * Computes the complementary distribution function of the hyperbolic
    * secant distribution with parameters @f$\mu@f$ and @f$\sigma@f$.
    */
   public static double barF (double mu, double sigma, double x) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");

      double y = (x - mu) / sigma;
      if (y >= ZLIMB)
         return 0.0;
      else if (y <= -ZLIMS)
      	 return 1.0;
      else
         return 2.0 / Math.PI * Math.atan (Math.exp (-Math.PI * y / 2.0));
   }

   /**
    * Computes the inverse of the hyperbolic secant distribution with
    * parameters @f$\mu@f$ and @f$\sigma@f$.
    */
   public static double inverseF (double mu, double sigma, double u) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u not in [0,1]");

      if (u >= 1.0)
         return Double.POSITIVE_INFINITY;
      else if (u <= 0.0)
         return Double.NEGATIVE_INFINITY;
      else
         return (mu + (2.0 * sigma / Math.PI * Math.log (Math.tan (Math.PI / 2.0 * u))));
   }

   /**
    * Estimates the parameters @f$(\mu, \sigma)@f$ of the hyperbolic
    * secant distribution using the maximum likelihood method, from the
    * @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1,…, n-1@f$. The
    * estimates are returned in a two-element array, in regular order:
    * [@f$\mu@f$, @f$\sigma@f$].  The estimate of the parameters is
    * given by maximizing numerically the log-likelihood function, using
    * the Uncmin package @cite iSCHa, @cite iVERa&thinsp;.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{\mu}@f$,
    * @f$\hat{\sigma}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      double sum;

      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      Optim system = new Optim (x, n);

      double[] parameters = new double[2];
      double[] xpls = new double[3];
      double[] param = new double[3];
      double[] fpls = new double[3];
      double[] gpls = new double[3];
      int[] itrcmd = new int[2];
      double[][] a = new double[3][3];
      double[] udiag = new double[3];

      sum = 0.0;
      for (int i = 0; i < n; i++)
         sum += x[i];
      param[1] = sum / (double) n;

      sum = 0.0;
      for (int i = 0; i < n; i++)
         sum += (x[i] - param[1]) * (x[i] - param[1]);
      param[2] = Math.sqrt (sum / (double) n);

      Uncmin_f77.optif0_f77 (2, param, system, xpls, fpls, gpls, itrcmd, a, udiag);

      for (int i = 0; i < 2; i++)
         parameters[i] = xpls[i+1];

      return parameters;
   }

   /**
    * Creates a new instance of a hyperbolic secant distribution with
    * parameters @f$\mu@f$ and @f$\sigma@f$ estimated using the maximum
    * likelihood method based on the @f$n@f$ observations @f$x[i]@f$, @f$i
    * = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static HyperbolicSecantDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new HyperbolicSecantDist (parameters[0], parameters[1]);
   }

   /**
    * Computes and returns the mean @f$E[X] = \mu@f$ of the hyperbolic
    * secant distribution with parameters @f$\mu@f$ and @f$\sigma@f$.
    *  @return the mean of the hyperbolic secant distribution @f$E[X] =
    * \mu@f$
    */
   public static double getMean (double mu, double sigma) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      return mu;
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] = \sigma^2@f$ of
    * the hyperbolic secant distribution with parameters @f$\mu@f$ and
    * @f$\sigma@f$.
    *  @return the variance of the hyperbolic secant distribution
    * @f$\mbox{Var}[X] = \sigma^2@f$
    */
   public static double getVariance (double mu, double sigma) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");

      return (sigma * sigma);
   }

   /**
    * Computes and returns the standard deviation of the hyperbolic secant
    * distribution with parameters @f$\mu@f$ and @f$\sigma@f$.
    *  @return the standard deviation of the hyperbolic secant
    * distribution
    */
   public static double getStandardDeviation (double mu, double sigma) {
      return Math.sqrt (HyperbolicSecantDist.getVariance (mu, sigma));
   }

   /**
    * Returns the parameter @f$\mu@f$ of this object.
    */
   public double getMu() {
      return mu;
   }

   /**
    * Returns the parameter @f$\sigma@f$ of this object.
    */
   public double getSigma() {
      return sigma;
   }

   /**
    * Sets the parameters @f$\mu@f$ and @f$\sigma@f$ of this object.
    */
   public void setParams (double mu, double sigma) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");

      this.mu = mu;
      this.sigma = sigma;
   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$\mu@f$,
    * @f$\sigma@f$].
    */
   public double[] getParams () {
      double[] retour = {mu, sigma};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : mu = " + mu + ", sigma = " + sigma;
   }

}