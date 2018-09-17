/*
 * Class:        LogisticDist
 * Description:  logistic distribution
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
 * Extends the class  @ref ContinuousDistribution for the *logistic*
 * distribution (e.g., @cite tJOH95b&thinsp; (page 115)). It has location
 * parameter @f$\alpha@f$ and scale parameter @f$\lambda> 0@f$. The density
 * is
 * @anchor REF_probdist_LogisticDist_eq_flogistic
 * @f[
 *   f (x) = \frac{\lambda e^{-\lambda(x - \alpha)}}{(1 + e^{-\lambda(x - \alpha)})^2} \qquad\qquad\mbox{for } -\infty< x < \infty, \tag{flogistic}
 * @f]
 * and the distribution function is
 * @anchor REF_probdist_LogisticDist_eq_Flogistic
 * @f[
 *   F(x) = \frac{1}{1 + e^{-\lambda(x - \alpha)}}  \qquad\qquad\mbox{for } -\infty< x < \infty. \tag{Flogistic}
 * @f]
 * For @f$\lambda=1@f$ and @f$\alpha=0@f$, one can write
 * @f[
 *   F(x) = \frac{1 + \tanh({x/2})}{2}.
 * @f]
 * The inverse distribution function is given by
 * @f[
 *   F^{-1}(u) = \ln(u/(1-u))/\lambda+ \alpha\qquad\mbox{for } 0 \le u < 1.
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class LogisticDist extends ContinuousDistribution {
   private double alpha;
   private double lambda;

   private static class Optim implements Lmder_fcn
   {
      protected double[] xi;
      protected int n;

      public Optim (double[] x, int n) {
         this.n = n;
         this.xi = new double[n];
         System.arraycopy (x, 0, this.xi, 0, n);
      }

      public void fcn (int m, int n, double[] x, double[] fvec, double[][] fjac, int iflag[])
      {
         if (x[2] <= 0.0) {
             final double BIG = 1.0e100;
             fvec[1] = BIG;
             fvec[2] = BIG;
             fjac[1][1] = BIG;
             fjac[1][2] = 0.0;
             fjac[2][1] = 0.0;
             fjac[2][2] = BIG;
             return;
         }

         double sum;
         double prod;

         if (iflag[1] == 1)
         {
            sum = 0.0;
            for (int i = 0; i < n; i++)
               sum += (1.0 / (1.0 + Math.exp (x[2] * (xi[i] - x[1]))));
            fvec[1] = sum - n / 2.0;

            sum = 0.0;
            for (int i = 0; i < n; i++)
            {
               prod = x[2] * (xi[i] - x[1]);
               sum -= prod * Math.tanh(prod/2.0);
            }
            fvec[2] = sum - n;
         }
         else if (iflag[1] == 2)
         {
            sum = 0.0;
            for (int i = 0; i < n; i++)
            {
               prod = Math.exp (x[2] * (xi[i] - x[1]));
               sum -= x[2] * prod / ((1 + prod) * (1 + prod));
            }
            fjac[1][1] = sum;

            sum = 0.0;
            for (int i = 0; i < n; i++)
            {
               prod = Math.exp (x[2] * (xi[i] - x[1]));
               sum -= (xi[i] - x[1])  * prod / ((1 + prod) * (1 + prod));
            }
            fjac[1][2] = sum;

            sum = 0.0;
            for (int i = 0; i < n; i++)
            {
               prod = Math.exp (x[2] * (xi[i] - x[1]));
               sum -= (x[2] * ((-1.0 + prod) * (1.0 + prod) - (2.0 * (x[2] * (xi[i] - x[1])) * prod))) / ((1.0 + prod) * (1.0 + prod));
            }
            fjac[2][1] = sum;

            sum = 0.0;
            for (int i = 0; i < n; i++)
            {
               prod = Math.exp (x[2] * (xi[i] - x[1]));
               sum -= ((x[1] - xi[1])  * ((-1.0 + prod) * (1.0 + prod) - (2.0 * (x[2] * (xi[i] - x[1])) * prod))) / ((1.0 + prod) * (1.0 + prod));
            }
            fjac[2][2] = sum;
         }
      }
   }

   /**
    * Constructs a `LogisticDist` object with default parameters
    * @f$\alpha= 0@f$ and @f$\lambda=1@f$.
    */
   public LogisticDist() {
      setParams (0.0, 1.0);
   }

   /**
    * Constructs a `LogisticDist` object with parameters @f$\alpha@f$ =
    * `alpha` and @f$\lambda@f$ = `lambda`.
    */
   public LogisticDist (double alpha, double lambda) {
      setParams (alpha, lambda);
   }


   public double density (double x) {
      return density (alpha, lambda, x);
   }

   public double cdf (double x) {
      return cdf (alpha, lambda, x);
   }

   public double barF (double x) {
      return barF (alpha, lambda, x);
   }

   public double inverseF (double u) {
      return inverseF (alpha, lambda, u);
   }

   public double getMean() {
      return LogisticDist.getMean (alpha, lambda);
   }

   public double getVariance() {
      return LogisticDist.getVariance (alpha, lambda);
   }

   public double getStandardDeviation() {
      return LogisticDist.getStandardDeviation (alpha, lambda);
   }

/**
 * Computes the density function @f$f(x)@f$.
 */
public static double density (double alpha, double lambda, double x) {
      if (lambda <= 0)
        throw new IllegalArgumentException ("lambda <= 0");
      double z = lambda * (x - alpha);
      if (z >= -100.0) {
         double v = Math.exp(-z);
         return lambda * v / ((1.0 + v)*(1.0 + v));
      }
      return lambda * Math.exp(z);
   }

   /**
    * Computes the distribution function @f$F(x)@f$.
    */
   public static double cdf (double alpha, double lambda, double x) {
      if (lambda <= 0)
        throw new IllegalArgumentException ("lambda <= 0");
      double z = lambda * (x - alpha);
      if (z >= -100.0)
         return 1.0 / (1.0 + Math.exp(-z));
      return Math.exp(z);
   }

   /**
    * Computes the complementary distribution function @f$1-F(x)@f$.
    */
   public static double barF (double alpha, double lambda, double x) {
      if (lambda <= 0)
        throw new IllegalArgumentException ("lambda <= 0");
      double z = lambda * (x - alpha);
      if (z <= 100.0)
         return 1.0 / (1.0 + Math.exp(z));
      return Math.exp(-z);
   }

   /**
    * Computes the inverse distribution function @f$F^{-1}(u)@f$.
    */
   public static double inverseF (double alpha, double lambda, double u) {
        if (lambda <= 0)
           throw new IllegalArgumentException ("lambda <= 0");
        if (u < 0.0 || u > 1.0)
           throw new IllegalArgumentException ("u not in [0, 1]");
        if (u >= 1.0)
            return Double.POSITIVE_INFINITY;
        if (u <= 0.0)
            return Double.NEGATIVE_INFINITY;

        return Math.log (u/(1.0 - u))/lambda + alpha;
   }

   /**
    * Estimates the parameters @f$(\alpha, \lambda)@f$ of the logistic
    * distribution using the maximum likelihood method, from the @f$n@f$
    * observations @f$x[i]@f$, @f$i = 0, 1,…, n-1@f$. The estimates are
    * returned in a two-element array, in regular order: [@f$\alpha@f$,
    * @f$\lambda@f$].  The maximum likelihood estimators are the values
    * @f$(\hat{\alpha}, \hat{\lambda})@f$ that satisfy the equations:
    * @f{align*}{
    *    \sum_{i=1}^n \frac{1}{1 + e^{\hat{\lambda} (x_i - \hat{\alpha})}} 
    *    & 
    *    = 
    *    \frac{n}{2}
    *    \\ 
    *    \sum_{i=1}^n \hat{\lambda} (x_i - \hat{\alpha}) \frac{1 - e^{\hat{\lambda} (x_i - \hat{\alpha})}}{1 + e^{\hat{\lambda} (x_i - \hat{\alpha})}} 
    *    & 
    *    = 
    *    n.
    * @f}
    * See @cite mEVA00a&thinsp; (page 128).
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameter [@f$\hat{\alpha}@f$,
    * @f$\hat{\lambda}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double sum = 0.0;
      for (int i = 0; i < n; i++)
         sum += x[i];

      double[] param = new double[3];
      param[1] = sum / (double) n;

      sum = 0.0;
      for (int i = 0; i < n; i++)
         sum += ((x[i] - param[1]) * (x[i] - param[1]));

      param[2] = Math.sqrt (Math.PI * Math.PI * n / (3.0 * sum));

      double[] fvec = new double [3];
      double[][] fjac = new double[3][3];
      int[] iflag = new int[2];
      int[] info = new int[2];
      int[] ipvt = new int[3];
      Optim system = new Optim (x, n);

      Minpack_f77.lmder1_f77 (system, 2, 2, param, fvec, fjac, 1e-5, info, ipvt);

      double parameters[] = new double[2];
      parameters[0] = param[1];
      parameters[1] = param[2];

      return parameters;
   }

   /**
    * Creates a new instance of a logistic distribution with parameters
    * @f$\alpha@f$ and @f$\lambda@f$ estimated using the maximum
    * likelihood method based on the @f$n@f$ observations @f$x[i]@f$, @f$i
    * = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static LogisticDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new LogisticDist (parameters[0], parameters[1]);
   }

   /**
    * Computes and returns the mean @f$E[X] = \alpha@f$ of the logistic
    * distribution with parameters @f$\alpha@f$ and @f$\lambda@f$.
    *  @return the mean of the logistic distribution @f$E[X] = \alpha@f$
    */
   public static double getMean (double alpha, double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return alpha;
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] = \pi^2
    * /(3\lambda^2)@f$ of the logistic distribution with parameters
    * @f$\alpha@f$ and @f$\lambda@f$.
    *  @return the variance of the logistic distribution @f$\mbox{Var}[X]
    * = 1 / 3 \pi^2 * (1 / \lambda^2)@f$
    */
   public static double getVariance (double alpha, double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return ((Math.PI * Math.PI / 3) * (1 / (lambda * lambda)));
   }

   /**
    * Computes and returns the standard deviation of the logistic
    * distribution with parameters @f$\alpha@f$ and @f$\lambda@f$.
    *  @return the standard deviation of the logistic distribution
    */
   public static double getStandardDeviation (double alpha, double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return (Math.sqrt(1.0 / 3.0) * Math.PI / lambda);
   }

   /**
    * Return the parameter @f$\alpha@f$ of this object.
    */
   public double getAlpha() {
      return alpha;
   }

   /**
    * Returns the parameter @f$\lambda@f$ of this object.
    */
   public double getLambda() {
      return lambda;
   }

   /**
    * Sets the parameters @f$\alpha@f$ and @f$\lambda@f$ of this object.
    */
   public void setParams (double alpha, double lambda) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      this.alpha  = alpha;
      this.lambda = lambda;
   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$\alpha@f$,
    * @f$\lambda@f$].
    */
   public double[] getParams () {
      double[] retour = {alpha, lambda};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : alpha = " + alpha + ", lambda = " + lambda;
   }

}