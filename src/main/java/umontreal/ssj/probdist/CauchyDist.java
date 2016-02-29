/*
 * Class:        CauchyDist
 * Description:  Cauchy distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        March 2009
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
import umontreal.ssj.util.Misc;
import optimization.*;

/**
 * Extends the class  @ref ContinuousDistribution for the *Cauchy*
 * distribution @cite tJOH95a&thinsp; (page 299) with location parameter
 * @f$\alpha@f$ and scale parameter @f$\beta> 0@f$. The density function is
 * given by
 * @anchor REF_probdist_CauchyDist_eq_fcuachy
 * @f[
 *   f (x) = \frac{\beta}{\pi[(x-\alpha)^2 + \beta^2]}, \qquad\qquad\mbox{for } -\infty< x < \infty. \tag{fcuachy}
 * @f]
 * The distribution function is
 * @f[
 *   F (x) = \frac{1}{2} + \frac{\arctan((x - \alpha)/\beta)}{\pi}, \qquad\qquad\mbox{for } -\infty< x < \infty,
 * @f]
 * and its inverse is
 * @f[
 *   F^{-1} (u) = \alpha+ \beta\tan(\pi(u - 1/2)). \qquad\mbox{for } 0 < u < 1.
 * @f]
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class CauchyDist extends ContinuousDistribution {
   private double alpha;
   private double beta;

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

         if (p[2] <= 0.0)               // barrier at 0
            return 1.0e200;

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
    * Constructs a `CauchyDist` object with parameters @f$\alpha=0@f$ and
    * @f$\beta=1@f$.
    */
   public CauchyDist() {
      setParams (0.0, 1.0);
   }

   /**
    * Constructs a `CauchyDist` object with parameters @f$\alpha=@f$
    * `alpha` and @f$\beta=@f$ `beta`.
    */
   public CauchyDist (double alpha, double beta) {
      setParams (alpha, beta);
   }


   public double density (double x) {
      return density (alpha, beta, x);
   }

   public double cdf (double x) {
      return cdf (alpha, beta, x);
   }

   public double barF (double x) {
      return barF (alpha, beta, x);
   }

   public double inverseF (double u){
      return inverseF (alpha, beta, u);
   }

   public double getMean() {
      return CauchyDist.getMean (alpha, beta);
   }

   public double getVariance() {
      return CauchyDist.getVariance (alpha, beta);
   }

   public double getStandardDeviation() {
      return CauchyDist.getStandardDeviation (alpha, beta);
   }

/**
 * Computes the density function.
 */
public static double density (double alpha, double beta, double x) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      double t = (x - alpha)/beta;
      return 1.0/(beta * Math.PI*(1 + t*t));
   }

   /**
    * Computes the distribution function.
    */
   public static double cdf (double alpha, double beta, double x) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      double z = (x - alpha)/beta;
      if (z < -0.5)
         return Math.atan(-1./z)/Math.PI;
      return Math.atan(z)/Math.PI + 0.5;
   }

   /**
    * Computes the complementary distribution.
    */
   public static double barF (double alpha, double beta, double x) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      double z = (x - alpha)/beta;
      if (z > 0.5)
         return Math.atan(1./z)/Math.PI;
      return 0.5 - Math.atan(z)/Math.PI;
   }

   /**
    * Computes the inverse of the distribution.
    */
   public static double inverseF (double alpha, double beta, double u) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
     if (u < 0.0 || u > 1.0)
        throw new IllegalArgumentException ("u must be in [0,1]");
     if (u <= 0.0)
        return Double.NEGATIVE_INFINITY;
     if (u >= 1.0)
        return Double.POSITIVE_INFINITY;
     if (u < 0.5)
        return alpha - 1.0/Math.tan (Math.PI*u) * beta;
     return alpha + Math.tan (Math.PI*(u - 0.5)) * beta;
   }

   /**
    * Estimates the parameters @f$(\alpha,\beta)@f$ of the Cauchy
    * distribution using the maximum likelihood method, from the @f$n@f$
    * observations @f$x[i]@f$, @f$i = 0, 1,…, n-1@f$. The estimates are
    * returned in a two-element array, in regular order: [@f$\alpha@f$,
    * @f$\beta@f$].  The estimates of the parameters are given by
    * maximizing numerically the log-likelihood function, using the Uncmin
    * package @cite iSCHa, @cite iVERa&thinsp;.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{\alpha}@f$,
    * @f$\hat{\beta}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      double sum = 0.0;

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

      param[1] = EmpiricalDist.getMedian (x, n);

      int m = Math.round ((float) n / 4.0f);
      double q3 = Misc.quickSelect (x, n, 3 * m);
      double q1 = Misc.quickSelect (x, n, m);
      param[2] = (q3 - q1) / 2.0;

      Uncmin_f77.optif0_f77 (2, param, system, xpls, fpls, gpls, itrcmd, a, udiag);

      for (int i = 0; i < 2; i++)
         parameters[i] = xpls[i+1];

      return parameters;
   }

   /**
    * Creates a new instance of a Cauchy distribution with parameters
    * @f$\alpha@f$ and @f$\beta@f$ estimated using the maximum
    * likelihood method based on the @f$n@f$ observations @f$x[i]@f$, @f$i
    * = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static CauchyDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new CauchyDist (parameters[0], parameters[1]);
   }

   /**
    * Throws an exception since the mean does not exist.
    *  @exception UnsupportedOperationException the mean of the Cauchy
    * distribution is undefined.
    */
   public static double getMean (double alpha, double beta) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");

      throw new UnsupportedOperationException("Undefined mean");
   }

   /**
    * Returns @f$\infty@f$ since the variance does not exist.
    *  @return @f$\infty@f$.
    */
   public static double getVariance (double alpha, double beta) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");

      return Double.POSITIVE_INFINITY;
   }

   /**
    * Returns @f$\infty@f$ since the standard deviation does not exist.
    *  @return @f$\infty@f$
    */
   public static double getStandardDeviation (double alpha, double beta) {
      return Double.POSITIVE_INFINITY;
   }

   /**
    * Returns the value of @f$\alpha@f$ for this object.
    */
   public double getAlpha() {
      return alpha;
   }

   /**
    * Returns the value of @f$\beta@f$ for this object.
    */
   public double getBeta() {
      return beta;
   }

   /**
    * Sets the value of the parameters @f$\alpha@f$ and @f$\beta@f$ for
    * this object.
    */
   public void setParams (double alpha, double beta) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      this.alpha = alpha;
      this.beta = beta;
   }

   /**
    * Return a table containing parameters of the current distribution.
    * This table is put in regular order: [@f$\alpha@f$, @f$\beta@f$].
    */
   public double[] getParams () {
      double[] retour = {alpha, beta};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : alpha = " + alpha + ", beta = " + beta;
   }

}