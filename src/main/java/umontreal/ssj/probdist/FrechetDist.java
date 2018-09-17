/*
 * Class:        FrechetDist
 * Description:  Fréchet distribution
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
import umontreal.ssj.util.*;
import optimization.*;
import umontreal.ssj.functions.MathFunction;

/**
 * Extends the class  @ref ContinuousDistribution for the <em>Fréchet</em>
 * distribution @cite tJOH95b&thinsp; (page 3), with location parameter
 * @f$\delta@f$, scale parameter @f$\beta> 0@f$, and shape parameter
 * @f$\alpha> 0@f$, where we use the notation @f$z = (x-\delta)/\beta@f$.
 * It has density
 * @f[
 *   f (x) = \frac{\alpha e^{-z^{-\alpha}}}{\beta z^{\alpha+1}}, \qquad\mbox{for } x > \delta
 * @f]
 * and distribution function
 * @f[
 *   F(x) = e^{-z^{-\alpha}}, \qquad\mbox{for } x > \delta.
 * @f]
 * Both the density and the distribution are 0 for @f$x \le\delta@f$.
 *
 * The mean is given by
 * @f[
 *   E[X] = \delta+ \beta\Gamma\!\left(1 - \frac{1}{\alpha}\right),
 * @f]
 * where @f$\Gamma(x)@f$ is the gamma function. The variance is
 * @f[
 *   \mbox{Var}[X] = \beta^2 \left[\Gamma\!\left(1 - \frac{2}{\alpha}\right) - \Gamma^2\!\left(1 - \frac{1}{\alpha}\right)\right].
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class FrechetDist extends ContinuousDistribution {
   private double delta;
   private double beta;
   private double alpha;


   private static class Optim implements Lmder_fcn {
      protected double[] x;
      protected int n;

      public Optim (double[] x, int n) {
         this.n = n;
         this.x = x;
      }

      public void fcn (int m, int n, double[] par, double[] fvec, double[][] fjac, int iflag[])
      {
         if (par[1] <= 0.0 || par[2] <= 0.0) {
            final double BIG = 1.0e100;
            fvec[1] = BIG;
            fvec[2] = BIG;
            fvec[3] = BIG;
            return;
         }

         double sum1, sum2, sumb, sum4, sum5;
         double z, w, v;
         double alpha = par[1];
         double beta = par[2];
         double mu = par[3];

         if (iflag[1] == 1) {
            sum1 = sum2 = sumb = sum4 = sum5 = 0;
            for (int i = 0; i < n; i++) {
               z = (x[i] - mu) / beta;
               sum1 += 1.0 / z;
               v = Math.pow(z, -alpha);
               sum2 += v / z;
               sumb += v;
               w = Math.log(z);
               sum4 += w;
               sum5 += v * w;
            }

            fvec[2] = sumb - n;   // eq. for beta
            fvec[3] = (alpha + 1) * sum1 - alpha * sum2;   // eq. for mu
            fvec[1] = n / alpha + sum5 - sum4;   // eq. for alpha

         } else if (iflag[1] == 2) {
            throw new IllegalArgumentException ("iflag = 2");
            // The 3 X 3 Jacobian must be calculated and put in fjac
         }
      }
   }


   private static class Function implements MathFunction {
      private int n;
      private double[] x;
      private double delta;
      public double sumxi;
      public double dif;

      public Function (double[] y, int n, double delta) {
         this.n = n;
         this.x = y;
         this.delta = delta;
         double xmin = Double.MAX_VALUE;
         for (int i = 0; i < n; i++) {
            if ((y[i] < xmin) && (y[i] > delta))
               xmin = y[i];
         }
         dif = xmin - delta;
      }

      public double evaluate (double alpha) {
         if (alpha <= 0.0) return 1.0e100;
         double v, w;
         double sum1 = 0, sum2 = 0, sum3 = 0;
         for (int i = 0; i < n; i++) {
            if (x[i] <= delta)
               continue;
            v = Math.log(x[i] - delta);
            w = Math.pow(dif / (x[i] - delta), alpha);
            sum1 += v;
            sum2 += w;
            sum3 += v * w;
         }

         sum1 /= n;
         sumxi = sum2 / n;
         return 1 / alpha + sum3 / sum2 - sum1;
      }
   }

   /**
    * Constructor for the standard <em>Fréchet</em> distribution with
    * parameters @f$\beta@f$ = 1 and @f$\delta@f$ = 0.
    */
   public FrechetDist (double alpha) {
      setParams (alpha, 1.0, 0.0);
   }

   /**
    * Constructs a `FrechetDist` object with parameters @f$\alpha@f$ =
    * `alpha`, @f$\beta@f$ = `beta` and @f$\delta@f$ = `delta`.
    */
   public FrechetDist (double alpha, double beta, double delta) {
      setParams (alpha, beta, delta);
   }


   public double density (double x) {
      return density (alpha, beta, delta, x);
   }

   public double cdf (double x) {
      return cdf (alpha, beta, delta, x);
   }

   public double barF (double x) {
      return barF (alpha, beta, delta, x);
   }

   public double inverseF (double u) {
      return inverseF (alpha, beta, delta, u);
   }

   public double getMean() {
      return getMean (alpha, beta, delta);
   }

   public double getVariance() {
      return getVariance (alpha, beta, delta);
   }

   public double getStandardDeviation() {
      return getStandardDeviation (alpha, beta, delta);
   }

/**
 * Computes and returns the density function.
 */
public static double density (double alpha, double beta, double delta,
                                 double x) {
      if (beta <= 0)
         throw new IllegalArgumentException ("beta <= 0");
      if (alpha <= 0)
         throw new IllegalArgumentException ("alpha <= 0");
      final double z = (x - delta)/beta;
      if (z <= 0.0)
         return 0.0;
      double t = Math.pow (z, -alpha);
      return  alpha * t * Math.exp (-t) / (z * beta);
   }

   /**
    * Computes and returns the distribution function.
    */
   public static double cdf (double alpha, double beta, double delta,
                             double x) {
      if (beta <= 0)
         throw new IllegalArgumentException ("beta <= 0");
      if (alpha <= 0)
         throw new IllegalArgumentException ("alpha <= 0");
      final double z = (x - delta)/beta;
      if (z <= 0.0)
         return 0.0;
      double t = Math.pow (z, -alpha);
      return  Math.exp (-t);
   }

   /**
    * Computes and returns the complementary distribution function @f$1 -
    * F(x)@f$.
    */
   public static double barF (double alpha, double beta, double delta,
                              double x) {
      if (beta <= 0)
         throw new IllegalArgumentException ("beta <= 0");
      if (alpha <= 0)
         throw new IllegalArgumentException ("alpha <= 0");
      final double z = (x - delta)/beta;
      if (z <= 0.0)
         return 1.0;
      double t = Math.pow (z, -alpha);
      return  -Math.expm1 (-t);
   }

   /**
    * Computes and returns the inverse distribution function.
    */
   public static double inverseF (double alpha, double beta, double delta,
                                  double u) {
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u not in [0, 1]");
      if (beta <= 0)
         throw new IllegalArgumentException ("beta <= 0");
      if (alpha <= 0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (u >= 1.0)
         return Double.POSITIVE_INFINITY;
      if (u <= 0.0)
         return delta;
      double t = Math.pow (-Math.log (u), 1.0/alpha);
      if (t <= Double.MIN_NORMAL)
         return Double.MAX_VALUE;
      return delta + beta / t;
   }

   /**
    * Given @f$\delta=@f$ `delta`, estimates the parameters @f$(\alpha,
    * \beta)@f$ of the <em>Fréchet</em> distribution using the maximum
    * likelihood method with the @f$n@f$ observations @f$x[i]@f$, @f$i =
    * 0, 1,…, n-1@f$. The estimates are returned in a two-element array,
    * in regular order: [@f$\alpha@f$, @f$\beta@f$].  The maximum
    * likelihood estimators are the values @f$(\hat{\alpha},
    * \hat{\beta})@f$ that satisfy the equations:
    * @f{align*}{
    *    \hat{\beta}
    *    & 
    *    = 
    *    \left(\frac{1}{n} \sum_{i=0}^{n-1} (x_i - \delta)^{-\hat{\alpha}}\right)^{\!\!-1/\hat{\alpha}} 
    *    \\ 
    *    \frac{1}{n} \sum_{i=0}^{n-1} \ln(x_i - \delta) 
    *    & 
    *    = 
    *    \frac{1}{\hat{\alpha}} + \frac{\sum_{i=0}^{n-1} (x_i - \delta)^{-\hat{\alpha}}\ln(x_i - \delta)}{\sum_{i=0}^{n-1} (x_i - \delta)^{-\hat{\alpha}}}.
    * @f}
    * @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @param delta        location parameter
    *  @return returns the parameters [@f$\hat{\alpha}@f$,
    * @f$\hat{\beta}@f$]
    */
   public static double[] getMLE (double[] x, int n, double delta) {
      if (n <= 1)
         throw new IllegalArgumentException ("n <= 1");

      Function func = new Function (x, n, delta);
      double a = 1e-4;
      double b = 1.0e12;
      double alpha = RootFinder.brentDekker (a, b, func, 1e-12);
      double par[] = new double[2];
      par[0] = alpha;
      par[1] = func.dif * Math.pow (func.sumxi, -1.0/alpha);
      return par;
   }

   /**
    * Given @f$\delta=@f$ `delta`, creates a new instance of a
    * <em>Fréchet</em> distribution with parameters @f$\alpha@f$ and
    * @f$\beta@f$ estimated using the maximum likelihood method based on
    * the @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @param delta        location parameter
    */
   public static FrechetDist getInstanceFromMLE (double[] x, int n,
                                                 double delta) {
      double par[] = getMLE (x, n, delta);
      return new FrechetDist (par[0], par[1], delta);
   }

   /**
    * Returns the mean of the <em>Fréchet</em> distribution with
    * parameters @f$\alpha@f$, @f$\beta@f$ and @f$\delta@f$.
    *  @return the mean
    */
   public static double getMean (double alpha, double beta, double delta) {
      if (beta <= 0)
         throw new IllegalArgumentException ("beta <= 0");
      if (alpha <= 1)
         throw new IllegalArgumentException ("alpha <= 1");
      double t = Num.lnGamma(1.0 - 1.0/alpha);
      return delta + beta * Math.exp(t);
   }

   /**
    * Returns the variance of the <em>Fréchet</em> distribution with
    * parameters @f$\alpha@f$, @f$\beta@f$ and @f$\delta@f$.
    *  @return the variance
    */
   public static double getVariance (double alpha, double beta,
                                     double delta) {
      if (beta <= 0)
         throw new IllegalArgumentException ("beta <= 0");
      if (alpha <= 2)
         throw new IllegalArgumentException ("alpha <= 2");
      double t = Num.lnGamma(1.0 - 1.0/alpha);
      double mu = Math.exp(t);
      double v = Math.exp(Num.lnGamma(1.0 - 2.0/alpha));
      return beta * beta * (v - mu * mu);
   }

   /**
    * Returns the standard deviation of the <em>Fréchet</em> distribution
    * with parameters @f$\alpha@f$, @f$\beta@f$ and @f$\delta@f$.
    *  @return the standard deviation
    */
   public static double getStandardDeviation (double alpha, double beta,
                                              double delta) {
      return  Math.sqrt(getVariance (alpha, beta, delta));
   }

   /**
    * Returns the parameter @f$\alpha@f$ of this object.
    */
   public double getAlpha() {
      return alpha;
   }

   /**
    * Returns the parameter @f$\beta@f$ of this object.
    */
   public double getBeta() {
      return beta;
   }

   /**
    * Returns the parameter @f$\delta@f$ of this object.
    */
   public double getDelta() {
      return delta;
   }

   /**
    * Sets the parameters @f$\alpha@f$, @f$\beta@f$ and @f$\delta@f$ of
    * this object.
    */
   public void setParams (double alpha, double beta, double delta) {
      if (beta <= 0)
         throw new IllegalArgumentException ("beta <= 0");
      if (alpha <= 0)
         throw new IllegalArgumentException ("alpha <= 0");
      this.delta  = delta;
      this.beta = beta;
      this.alpha = alpha;
   }

   /**
    * Return an array containing the parameters of the current object in
    * regular order: [@f$\alpha@f$, @f$\beta@f$, @f$\delta@f$].
    */
   public double[] getParams() {
      double[] retour = {alpha, beta, delta};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : alpha = " + alpha + ", beta = " + beta + ", delta = " + delta;
   }

}