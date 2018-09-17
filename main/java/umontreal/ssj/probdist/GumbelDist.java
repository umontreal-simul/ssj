/*
 * Class:        GumbelDist
 * Description:  Gumbel distribution
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
import umontreal.ssj.util.Num;
import umontreal.ssj.util.RootFinder;
import umontreal.ssj.functions.MathFunction;

/**
 * Extends the class  @ref ContinuousDistribution for the *Gumbel*
 * distribution @cite tJOH95b&thinsp; (page 2), with location parameter
 * @f$\delta@f$ and scale parameter @f$\beta\neq0@f$. Using the notation
 * @f$z = (x-\delta)/\beta@f$, it has density
 * @anchor REF_probdist_GumbelDist_eq_densgumbel
 * @f[
 *   f (x) = \frac{e^{-z} e^{-e^{-z}}}{|\beta|}, \qquad\mbox{for } -\infty< x < \infty
 *    \tag{densgumbel}
 * @f]
 * and distribution function
 * @f[
 *   F(x) = \left\{ \begin{array}{ll}
 *    e^{-e^{-z}}, \qquad
 *    & 
 *    \mbox{for } \beta> 0 
 *    \\ 
 *   1 - e^{-e^{-z}}, \qquad
 *    & 
 *    \mbox{for } \beta< 0. 
 *   \end{array} \right.
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class GumbelDist extends ContinuousDistribution {
   private double delta;
   private double beta;

   private static class FunctionPlus implements MathFunction {
      // when beta > 0
      protected int n;
      protected double mean;
      protected double[] x;
      private double minx;   // min of all {x[i]}

      public FunctionPlus (double[] y, int n, double mean, double minx) {
         this.n = n;
         this.mean = mean;
         this.x = y;
         this.minx = minx;
      }

      public double evaluate (double lam) {
         if (lam <= 0.0) return 1.0e100;
         double tem;
         double sum2 = 0.0;
         double sum1 = 0.0;

         for (int i = 0; i < n; i++) {
            tem = Math.exp (-(x[i] - minx)* lam);
            sum1 += tem;
            sum2 += x[i] * tem;
         }

         return (mean - 1/lam) * sum1 - sum2;
      }
   }


   private static class FunctionMinus implements MathFunction {
      // when beta < 0
      protected int n;
      protected double mean;
      protected double[] x;
      protected double xmax;

      public FunctionMinus (double[] y, int n, double mean, double xmax) {
         this.n = n;
         this.mean = mean;
         this.x = y;
         this.xmax = xmax;
      }

      public double evaluate (double lam) {
         if (lam >= 0.0) return 1.0e100;
         double tem;
         double sum2 = 0.0;
         double sum1 = 0.0;

         for (int i = 0; i < n; i++) {
            tem = Math.exp ((xmax - x[i]) * lam);
            sum1 += tem;
            sum2 += x[i] * tem;
         }

         return (mean - 1/lam) * sum1 - sum2;
      }
   }

   /**
    * Constructor for the standard Gumbel distribution with parameters
    * @f$\beta@f$ = 1 and @f$\delta@f$ = 0.
    */
   public GumbelDist() {
      setParams (1.0, 0.0);
   }

   /**
    * Constructs a `GumbelDist` object with parameters @f$\beta@f$ =
    * `beta` and @f$\delta@f$ = `delta`.
    */
   public GumbelDist (double beta, double delta) {
      setParams (beta, delta);
   }


   public double density (double x) {
      return density (beta, delta, x);
   }

   public double cdf (double x) {
      return cdf (beta, delta, x);
   }

   public double barF (double x) {
      return barF (beta, delta, x);
   }

   public double inverseF (double u) {
      return inverseF (beta, delta, u);
   }

   public double getMean() {
      return GumbelDist.getMean (beta, delta);
   }

   public double getVariance() {
      return GumbelDist.getVariance (beta, delta);
   }

   public double getStandardDeviation() {
      return GumbelDist.getStandardDeviation (beta, delta);
   }

/**
 * Computes and returns the density function.
 */
public static double density (double beta, double delta, double x) {
      if (beta == 0.)
         throw new IllegalArgumentException ("beta = 0");
      final double z = (x - delta)/beta;
      if (z <= -10.0)
         return 0.0;
      double t = Math.exp (-z);
      return  t * Math.exp (-t)/Math.abs(beta);
   }

   /**
    * Computes and returns the distribution function.
    */
   public static double cdf (double beta, double delta, double x) {
      if (beta == 0.)
         throw new IllegalArgumentException ("beta = 0");
      final double z = (x - delta)/beta;
      if (beta > 0.) {
         if (z <= -7.0)
            return 0.0;
         return Math.exp (-Math.exp (-z));
      } else {   // beta < 0
          if (z <= -7.0)
            return 1.0;
         return -Math.expm1 (-Math.exp (-z));
     }
   }

   /**
    * Computes and returns the complementary distribution function @f$1 -
    * F(x)@f$.
    */
   public static double barF (double beta, double delta, double x) {
      if (beta == 0.)
         throw new IllegalArgumentException ("beta = 0");
      final double z = (x - delta)/beta;
      if (beta > 0.) {
         if (z <= -7.0)
            return 1.0;
         return -Math.expm1 (-Math.exp (-z));
      } else {   // beta < 0
          if (z <= -7.0)
            return 0.0;
         return Math.exp (-Math.exp (-z));
      }
   }

   /**
    * Computes and returns the inverse distribution function.
    */
   public static double inverseF (double beta, double delta, double u) {
       if (u < 0.0 || u > 1.0)
          throw new IllegalArgumentException ("u not in [0, 1]");
      if (beta == 0.)
         throw new IllegalArgumentException ("beta = 0");
       if (u >= 1.0)
           return Double.POSITIVE_INFINITY;
       if (u <= 0.0)
           return Double.NEGATIVE_INFINITY;
       if (beta > 0.)
          return delta - Math.log (-Math.log (u))*beta;
       else
          return delta - Math.log (-Math.log1p(-u))*beta;
   }

   /**
    * Estimates the parameters @f$(\beta,\delta)@f$ of the Gumbel
    * distribution, <em>assuming that @f$\beta> 0@f$</em>, and using the
    * maximum likelihood method with the @f$n@f$ observations @f$x[i]@f$,
    * @f$i = 0, 1,…, n-1@f$. The estimates are returned in a two-element
    * array, in regular order: [@f$\beta@f$, @f$\delta@f$].  The maximum
    * likelihood estimators are the values @f$(\hat{\beta},
    * \hat{\delta})@f$ that satisfy the equations:
    * @f{align*}{
    *    \hat{\beta}
    *    & 
    *    = 
    *    \bar{x}_n - \frac{\sum_{i=1}^n x_i  e^{- x_i/\hat{\beta}}}{\sum_{i=1}^n e^{- x_i / \hat{\beta}}}
    *    \\ 
    *    \hat{\delta} 
    *    & 
    *    = 
    *    -{\hat{\beta}} \ln\left( \frac{1}{n} \sum_{i=1}^n e^{-x_i/\hat{\beta}} \right),
    * @f}
    * where @f$\bar{x}_n@f$ is the average of @f$x[0],…,x[n-1]@f$.
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{\delta}@f$,
    * @f$\hat{\beta}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      if (n <= 1)
         throw new IllegalArgumentException ("n <= 1");
      int i;
      double par[] = new double[2];

      double xmin = Double.MAX_VALUE;
      double sum = 0;
      for (i = 0; i < n; i++) {
         sum += x[i];
         if (x[i] < xmin)
            xmin = x[i];
      }
      double mean = sum / (double) n;

      sum = 0;
      for (i = 0; i < n; i++)
         sum += (x[i] - mean) * (x[i] - mean);
      double variance = sum / (n - 1.0);

      FunctionPlus func = new FunctionPlus (x, n, mean, xmin);

      double lam = 1.0 / (0.7797*Math.sqrt (variance));
      final double EPS = 0.02;
      double a = (1.0 - EPS)*lam - 5.0;
      if (a <= 0)
         a = 1e-15;
      double b = (1.0 + EPS)*lam + 5.0;
      lam = RootFinder.brentDekker (a, b, func, 1e-8);
      par[0] = 1.0 / lam;

      sum = 0;
      for (i = 0; i < n; i++)
           sum += Math.exp (-(x[i] - xmin) * lam);
      par[1] = xmin - Math.log (sum/n) / lam;
      return par;
   }

   /**
    * Similar to  #getMLE, but <em>for the case @f$\beta< 0@f$</em>.
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{\delta}@f$,
    * @f$\hat{\beta}@f$]
    */
   public static double[] getMLEmin (double[] x, int n) {
      if (n <= 1)
         throw new IllegalArgumentException ("n <= 1");

      int i;
      double par[] = new double[2];
      double xmax = -Double.MAX_VALUE;
      double sum = 0.0;
      for (i = 0; i < n; i++) {
         sum += x[i];
         if (x[i] > xmax)
            xmax = x[i];
      }
      double mean = sum / (double) n;

      sum = 0.0;
      for (i = 0; i < n; i++)
         sum += (x[i] - mean) * (x[i] - mean);
      double variance = sum / (n - 1.0);

      FunctionMinus func = new FunctionMinus (x, n, mean, xmax);

      double lam = -1.0 / (0.7797*Math.sqrt (variance));
      final double EPS = 0.02;
      double a = (1.0 + EPS)*lam - 2.0;
      double b = (1.0 - EPS)*lam + 2.0;
      if (b >= 0)
         b = -1e-15;
      lam = RootFinder.brentDekker (a, b, func, 1e-12);
      par[0] = 1.0 / lam;

      sum = 0.0;
      for (i = 0; i < n; i++)
         sum += Math.exp ((xmax - x[i]) * lam);
      par[0] = 1.0 / lam;
      par[1] = xmax - Math.log (sum / n) / lam;

      return par;
   }

   /**
    * Creates a new instance of an Gumbel distribution with parameters
    * @f$\beta@f$ and @f$\delta@f$ estimated using the maximum
    * likelihood method based on the @f$n@f$ observations @f$x[i]@f$, @f$i
    * = 0, 1, …, n-1@f$, <em>assuming that @f$\beta> 0@f$</em>.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static GumbelDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new GumbelDist (parameters[0], parameters[1]);
   }

   /**
    * Similar to  #getInstanceFromMLE, but <em>for the case @f$\beta<
    * 0@f$</em>.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static GumbelDist getInstanceFromMLEmin (double[] x, int n) {
      double parameters[] = getMLEmin (x, n);
      return new GumbelDist (parameters[0], parameters[1]);
   }

   /**
    * Returns the mean, @f$E[X] = \delta+ \gamma\beta@f$, of the Gumbel
    * distribution with parameters @f$\beta@f$ and @f$\delta@f$, where
    * @f$\gamma= 0.5772156649015329@f$ is the Euler-Mascheroni constant.
    *  @return the mean of the Extreme Value distribution @f$E[X] =
    * \delta+ \gamma* \beta@f$
    */
   public static double getMean (double beta, double delta) {
     if (beta == 0.0)
         throw new IllegalArgumentException ("beta = 0");

      return delta + Num.EULER * beta;
   }

   /**
    * Returns the variance @f$\mbox{Var}[X] = \pi^2 \beta^2\!/6@f$ of
    * the Gumbel distribution with parameters @f$\beta@f$ and
    * @f$\delta@f$.
    *  @return the variance of the Gumbel distribution @f$\mbox{Var}[X] =
    * ()\pi\beta)^2/6@f$
    */
   public static double getVariance (double beta, double delta) {
     if (beta == 0.0)
         throw new IllegalArgumentException ("beta = 0");

      return Math.PI * Math.PI * beta * beta / 6.0;
   }

   /**
    * Returns the standard deviation of the Gumbel distribution with
    * parameters @f$\beta@f$ and @f$\delta@f$.
    *  @return the standard deviation of the Gumbel distribution
    */
   public static double getStandardDeviation (double beta, double delta) {
     if (beta == 0.0)
         throw new IllegalArgumentException ("beta = 0");

      return  Math.sqrt(getVariance (beta, delta));
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
    * Sets the parameters @f$\beta@f$ and @f$\delta@f$ of this object.
    */
   public void setParams (double beta, double delta) {
     if (beta == 0)
         throw new IllegalArgumentException ("beta = 0");
      this.delta  = delta;
      this.beta = beta;
   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$\beta@f$,
    * @f$\delta@f$].
    */
   public double[] getParams () {
      double[] retour = {beta, delta};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : beta = " + beta + ", delta = " + delta;
   }

}