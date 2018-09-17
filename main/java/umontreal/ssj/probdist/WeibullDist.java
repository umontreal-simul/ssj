/*
 * Class:        WeibullDist
 * Description:  Weibull distribution
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
 * This class extends the class  @ref ContinuousDistribution for the
 * *Weibull* distribution @cite tJOH95a&thinsp; (page 628) with shape
 * parameter @f$\alpha> 0@f$, location parameter @f$\delta@f$, and scale
 * parameter @f$\lambda> 0@f$. The density function is
 * @anchor REF_probdist_WeibullDist_eq_fweibull
 * @f[
 *   f(x) = \alpha\lambda^{\alpha}(x-\delta)^{\alpha-1} e^{-(\lambda(x-\delta))^{\alpha}} \qquad\mbox{for }x>\delta, \tag{fweibull}
 * @f]
 * the distribution function is
 * @anchor REF_probdist_WeibullDist_eq_Fweibull
 * @f[
 *   F(x) = 1 - e^{-(\lambda(x - \delta))^{\alpha}} \qquad\mbox{for }x>\delta, \tag{Fweibull}
 * @f]
 * and the inverse distribution function is
 * @f[
 *   F^{-1}(u) = (-\ln(1-u))^{1/\alpha}/\lambda+ \delta\qquad\mbox{for } 0 \le u < 1.
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class WeibullDist extends ContinuousDistribution {
   private double alpha;
   private double lambda;
   private double delta;

   private static class Function implements MathFunction {
      private int n;
      private double xi[];
      private double lnXi[];
      private double sumLnXi = 0.0;
      private final double LN_EPS = Num.LN_DBL_MIN - Num.LN2;

      public Function (double x[], int n)
      {
         this.n = n;
         this.xi = new double[n];
         this.lnXi = new double[n];

         for (int i = 0; i < n; i++)
         {
            this.xi[i] = x[i];
            if (x[i] > 0.0)
               this.lnXi[i] = Math.log (x[i]);
            else
               this.lnXi[i] = LN_EPS;
            sumLnXi += this.lnXi[i];
         }
      }

      public double evaluate (double x)
      {
         if (x <= 0.0) return 1.0e200;
         double sumXiLnXi = 0.0;
         double sumXi = 0.0;
         double xalpha;

         for (int i = 0; i < n; i++)
         {
            xalpha = Math.pow (this.xi[i], x);
            sumXiLnXi += xalpha * lnXi[i];
            sumXi += xalpha;
         }

         return (x * (n * sumXiLnXi - sumLnXi * sumXi) - n * sumXi);
      }
   }

   /**
    * Constructs a `WeibullDist` object with parameters @f$\alpha@f$ =
    * `alpha`, @f$\lambda@f$ = 1, and @f$\delta@f$ = 0.
    */
   public WeibullDist (double alpha) {
      setParams (alpha, 1.0, 0.0);
   }

   /**
    * Constructs a `WeibullDist` object with parameters @f$\alpha=@f$
    * `alpha`, @f$\lambda@f$ = `lambda`, and @f$\delta@f$ = `delta`.
    */
   public WeibullDist (double alpha, double lambda, double delta) {
      setParams (alpha, lambda, delta);
   }


   public double density (double x) {
      return density (alpha, lambda, delta, x);
   }

   public double cdf (double x) {
      return cdf (alpha, lambda, delta, x);
   }

   public double barF (double x) {
      return barF (alpha, lambda, delta, x);
   }

   public double inverseF (double u) {
      return inverseF (alpha, lambda, delta, u);
   }

   public double getMean() {
      return WeibullDist.getMean (alpha, lambda, delta);
   }

   public double getVariance() {
      return WeibullDist.getVariance (alpha, lambda, delta);
   }

   public double getStandardDeviation() {
      return WeibullDist.getStandardDeviation (alpha, lambda, delta);
   }

/**
 * Computes the density function.
 */
public static double density (double alpha, double lambda,
                                 double delta, double x) {
      if (alpha <= 0.0)
        throw new IllegalArgumentException ("alpha <= 0");
      if (lambda <= 0.0)
        throw new IllegalArgumentException ("lambda <= 0");
      if (x <= delta)
         return 0.0;
      double y = Math.log(lambda*(x - delta)) * alpha;
      if (y >= 7.0)
         return 0.0;
      y = Math.exp(y);

      return alpha * (y / (x - delta)) * Math.exp(-y);
   }

   /**
    * Same as `density (alpha, 1, 0, x)`.
    */
   public static double density (double alpha, double x) {
      return density (alpha, 1.0, 0.0, x);
   }

   /**
    * Computes the distribution function.
    */
   public static double cdf (double alpha, double lambda,
                             double delta, double x) {
      if (alpha <= 0.0)
        throw new IllegalArgumentException ("alpha <= 0");
      if (lambda <= 0.0)
        throw new IllegalArgumentException ("lambda <= 0");
      if (x <= delta)
         return 0.0;
      if ((lambda*(x - delta) >= XBIG) && (alpha >= 1.0))
         return 1.0;
      double y = Math.log(lambda*(x - delta)) * alpha;
      if (y >= 3.65)
         return 1.0;
      y = Math.exp(y);
      return -Math.expm1 (-y);   // in JDK-1.5
   }

   /**
    * Same as `cdf (alpha, 1, 0, x)`.
    */
   public static double cdf (double alpha, double x) {
      return cdf (alpha, 1.0, 0.0, x);
   }

   /**
    * Computes the complementary distribution function.
    */
   public static double barF (double alpha, double lambda,
                              double delta, double x) {
      if (alpha <= 0)
        throw new IllegalArgumentException ("alpha <= 0");
      if (lambda <= 0)
        throw new IllegalArgumentException ("lambda <= 0");
      if (x <= delta)
         return 1.0;
      if (alpha >= 1.0 && x >= Num.DBL_MAX_EXP*2)
         return 0.0;

      double temp = Math.log (lambda*(x - delta)) * alpha;
      if (temp >= Num.DBL_MAX_EXP * Num.LN2)
         return 0.0;
      temp = Math.exp(temp);
      return Math.exp (-temp);
   }

   /**
    * Same as `barF (alpha, 1, 0, x)`.
    */
   public static double barF (double alpha, double x) {
      return barF (alpha, 1.0, 0.0, x);
   }

   /**
    * Computes the inverse of the distribution function.
    */
   public static double inverseF (double alpha, double lambda,
                                  double delta, double u) {
        double t;
        if (alpha <= 0.0)
            throw new IllegalArgumentException ("alpha <= 0");
        if (lambda <= 0.0)
            throw new IllegalArgumentException ("lambda <= 0");

        if (u < 0.0 || u > 1.0)
            throw new IllegalArgumentException ("u not in [0, 1]");
        if (u <= 0.0)
           return 0.0;
        if (u >= 1.0)
           return Double.POSITIVE_INFINITY;

        t = -Math.log1p (-u);
        if (Math.log (t)/Math.log (10) >= alpha*Num.DBL_MAX_10_EXP)
           throw new ArithmeticException
              ("inverse function cannot be positive infinity");

        return Math.pow (t, 1.0/alpha)/lambda + delta;
   }

   /**
    * Same as `inverseF (alpha, 1, 0, x)`.
    */
   public static double inverseF (double alpha, double x) {
      return inverseF (alpha, 1.0, 0.0, x);
   }

   private static double[] getMaximumLikelihoodEstimate (double[] x, int n,
                                                         double delta) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      if (delta != 0.0)
         throw new IllegalArgumentException ("delta must be equal to 0");
// Verifier cette fonction si delta != 0.

      final double LN_EPS = Num.LN_DBL_MIN - Num.LN2;
      double sumLn = 0.0;
      double sumLn2 = 0.0;
      double lnxi;
      for (int i = 0; i < n; i++) {
         if (x[i] <= delta)
            lnxi = LN_EPS;
         else
            lnxi = Math.log (x[i]);
         sumLn += lnxi;
         sumLn2 += lnxi * lnxi;
      }

      double alpha0 = Math.sqrt ((double) n / ((6.0 / (Math.PI * Math.PI)) *
                  (sumLn2 - sumLn * sumLn / (double) n)));
      double a = alpha0 - 20.0;
      if (a <= delta)
         a = delta + 1.0e-5;

      double param[] = new double[3];
      param[2] = 0.0;
      Function f = new Function (x, n);
      param[0] = RootFinder.brentDekker (a, alpha0 + 20.0, f, 1e-5);

      double sumXalpha = 0.0;
      for (int i = 0; i < n; i++)
         sumXalpha += Math.pow (x[i], param[0]);
      param[1] = Math.pow ((double) n / sumXalpha, 1.0 / param[0]);

      return param;
   }

/**
 * Estimates the parameters @f$(\alpha, \lambda)@f$ of the Weibull
 * distribution, assuming that @f$\delta= 0@f$, using the maximum likelihood
 * method, from the @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$.
 * The estimates are returned in a two-element array, in regular order:
 * [@f$\alpha@f$, @f$\lambda@f$].  The maximum likelihood estimators are
 * the values @f$(\hat{\alpha}@f$, @f$\hat{\lambda})@f$ that satisfy the
 * equations
 * @f{align*}{
 *    \frac{\sum_{i=1}^n x_i^{\hat{\alpha}} \ln(x_i)}{\sum_{i=1}^n x_i^{\hat{\alpha}}} - \frac{1}{\hat{\alpha}} 
 *    & 
 *    = 
 *    \frac{\sum_{i=1}^n \ln(x_i)}{n}
 *    \\ 
 *   \hat{\lambda} 
 *    & 
 *    = 
 *    \left( \frac{n}{\sum_{i=1}^n x_i^{\hat{\alpha}}} \right)^{1/\hat{\alpha}}
 * @f}
 * See @cite sLAW00a&thinsp; (page 303).
 *  @param x            the list of observations to use to evaluate
 *                      parameters
 *  @param n            the number of observations to use to evaluate
 *                      parameters
 *  @return returns the parameter [@f$\hat{\alpha}@f$, @f$\hat{\lambda}@f$,
 * @f$\hat{\delta}@f$ = 0]
 */
public static double[] getMLE (double[] x, int n)
   {
      return getMaximumLikelihoodEstimate (x, n, 0.0);
   }

   /**
    * Creates a new instance of a Weibull distribution with parameters
    * @f$\alpha@f$, @f$\lambda@f$ and @f$\delta= 0@f$ estimated using
    * the maximum likelihood method based on the @f$n@f$ observations
    * @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static WeibullDist getInstanceFromMLE (double[] x, int n) {
      double param[] = getMLE (x, n);
      return new WeibullDist (param[0], param[1], param[2]);
   }

   /**
    * Computes and returns the mean
    *  @f$E[X] = \delta+ \Gamma(1 + 1/\alpha)/\lambda@f$
    *  of the Weibull distribution with parameters @f$\alpha@f$,
    * @f$\lambda@f$ and @f$\delta@f$.
    *  @return the mean of the Weibull distribution @f$E[X] = \delta+
    * \Gamma(1 + 1/\alpha) / \lambda@f$
    */
   public static double getMean (double alpha, double lambda, double delta) {
      if (alpha <= 0.0)
        throw new IllegalArgumentException ("alpha <= 0");
      if (lambda <= 0.0)
        throw new IllegalArgumentException ("lambda <= 0");

      return (delta + Math.exp (Num.lnGamma(1.0 + 1.0 / alpha)) / lambda);
   }

   /**
    * Computes and returns the variance
    *  @f$\mbox{Var}[X] = | \Gamma(2/\alpha+ 1) - \Gamma^2(1/\alpha+
    * 1) | /\lambda^2@f$
    *  of the Weibull distribution with parameters @f$\alpha@f$,
    * @f$\lambda@f$ and @f$\delta@f$.
    *  @return the variance of the Weibull distribution @f$\mbox{Var}[X] =
    * 1 / \lambda^2 | \Gamma(2/\alpha+ 1) - \Gamma^2(1/\alpha+ 1)
    * |@f$
    */
   public static double getVariance (double alpha, double lambda,
                                     double delta) {
      double gAlpha;

      if (alpha <= 0.0)
        throw new IllegalArgumentException ("alpha <= 0");
      if (lambda <= 0.0)
        throw new IllegalArgumentException ("lambda <= 0");

      gAlpha = Math.exp (Num.lnGamma (1.0 / alpha + 1.0));

      return (Math.abs (Math.exp (Num.lnGamma(2 / alpha + 1)) - gAlpha * gAlpha) / (lambda * lambda));
   }

   /**
    * Computes and returns the standard deviation of the Weibull
    * distribution with parameters @f$\alpha@f$, @f$\lambda@f$ and
    * @f$\delta@f$.
    *  @return the standard deviation of the Weibull distribution
    */
   public static double getStandardDeviation (double alpha, double lambda,
                                              double delta) {
      return Math.sqrt (WeibullDist.getVariance (alpha, lambda, delta));
   }

   /**
    * Returns the parameter @f$\alpha@f$.
    */
   public double getAlpha() {
      return alpha;
   }

   /**
    * Returns the parameter @f$\lambda@f$.
    */
   public double getLambda() {
      return lambda;
   }

   /**
    * Returns the parameter @f$\delta@f$.
    */
   public double getDelta() {
      return delta;
   }

   /**
    * Sets the parameters @f$\alpha@f$, @f$\lambda@f$ and @f$\delta@f$
    * for this object.
    */
   public void setParams (double alpha, double lambda, double delta) {
      if (alpha <= 0.0)
        throw new IllegalArgumentException ("alpha <= 0");
      if (lambda <= 0.0)
        throw new IllegalArgumentException ("lambda <= 0");

      this.alpha  = alpha;
      this.lambda = lambda;
      this.delta  = delta;
      supportA = delta;
   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$\alpha@f$,
    * @f$\lambda@f$, @f$\delta@f$].
    */
   public double[] getParams () {
      double[] retour = {alpha, lambda, delta};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : alpha = " + alpha + ", lambda = " + lambda + ", delta = " + delta;
   }

}