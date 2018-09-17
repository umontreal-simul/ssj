/*
 * Class:        ExtremeValueDist
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
 * <strong>This class has been replaced by  @ref GumbelDist </strong>.
 *
 * Extends the class  @ref ContinuousDistribution for the *extreme value* (or
 * <em>Gumbel</em>) distribution @cite tJOH95b&thinsp; (page 2), with
 * location parameter @f$\alpha@f$ and scale parameter @f$\lambda> 0@f$. It
 * has density
 * @anchor REF_probdist_ExtremeValueDist_eq_fextremevalue
 * @f[
 *   f (x) = \lambda e^{-\lambda(x-\alpha)} e^{-e^{-\lambda(x-\alpha)}}, \qquad\qquad\mbox{for } -\infty< x < \infty, \tag{fextremevalue}
 * @f]
 * distribution function
 * @anchor REF_probdist_ExtremeValueDist_eq_Fextreme
 * @f[
 *   F(x) = e^{-e^{-\lambda(x - \alpha)}} \qquad\qquad\mbox{for } -\infty< x < \infty, \tag{Fextreme}
 * @f]
 * and inverse distribution function
 * @f[
 *   F^{-1}(u) = -\ln(-\ln(u))/\lambda+ \alpha, \qquad\mbox{for } 0 \le u \le1.
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
@Deprecated
public class ExtremeValueDist extends ContinuousDistribution {
   private double alpha;
   private double lambda;

   private static class Function implements MathFunction {
      protected int n;
      protected double mean;
      protected double[] x;

      public Function (double[] x, int n, double mean) {
         this.n = n;
         this.mean = mean;
         this.x = new double[n];

         System.arraycopy(x, 0, this.x, 0, n);
      }

      public double evaluate (double lambda) {
         if (lambda <= 0.0) return 1.0e200;
         double exp = 0.0;
         double sumXiExp = 0.0;
         double sumExp = 0.0;

         for (int i = 0; i < n; i++)
         {
            exp = Math.exp (-x[i] * lambda);
            sumExp += exp;
            sumXiExp += x[i] * exp;
         }

         return ((mean - 1.0 / lambda) * sumExp - sumXiExp);
      }
   }

   /**
    * <strong>THIS CLASS HAS BEEN REPLACED BY  @ref GumbelDist </strong>.
    * Constructs a `ExtremeValueDist` object with parameters @f$\alpha@f$
    * = 0 and @f$\lambda@f$ = 1.
    */
   public ExtremeValueDist() {
      setParams (0.0, 1.0);
   }

   /**
    * <strong>THIS CLASS HAS BEEN REPLACED BY  @ref GumbelDist </strong>.
    * Constructs a `ExtremeValueDist` object with parameters @f$\alpha@f$
    * = `alpha` and @f$\lambda@f$ = `lambda`.
    */
   public ExtremeValueDist (double alpha, double lambda) {
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
      return ExtremeValueDist.getMean (alpha, lambda);
   }

   public double getVariance() {
      return ExtremeValueDist.getVariance (alpha, lambda);
   }

   public double getStandardDeviation() {
      return ExtremeValueDist.getStandardDeviation (alpha, lambda);
   }

/**
 * Computes the density function.
 */
public static double density (double alpha, double lambda, double x) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      final double z = lambda*(x - alpha);
      if (z <= -10.0)
         return 0.0;
      double t = Math.exp (-z);
      return lambda * t * Math.exp (-t);
   }

   /**
    * <strong>THIS CLASS HAS BEEN REPLACED BY  @ref GumbelDist </strong>.
    * Computes the distribution function.
    */
   public static double cdf (double alpha, double lambda, double x) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      final double z = lambda*(x - alpha);
      if (z <= -10.0)
         return 0.0;
      if (z >= XBIG)
         return 1.0;
      return Math.exp (-Math.exp (-z));
   }

   /**
    * Computes the complementary distribution function.
    */
   public static double barF (double alpha, double lambda, double x) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      final double z = lambda*(x - alpha);
      if (z <= -10.0)
         return 1.0;
      if (z >= XBIGM)
         return 0.0;
      return -Math.expm1 (-Math.exp (-z));
   }

   /**
    * Computes the inverse distribution function.
    */
   public static double inverseF (double alpha, double lambda, double u) {
       if (u < 0.0 || u > 1.0)
          throw new IllegalArgumentException ("u not in [0, 1]");
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
       if (u >= 1.0)
           return Double.POSITIVE_INFINITY;
       if (u <= 0.0)
           return Double.NEGATIVE_INFINITY;

       return -Math.log (-Math.log (u))/lambda+alpha;
   }

   /**
    * Estimates the parameters @f$(\alpha,\lambda)@f$ of the extreme
    * value distribution using the maximum likelihood method, from the
    * @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1,…, n-1@f$. The
    * estimates are returned in a two-element array, in regular order:
    * [@f$\alpha@f$, @f$\lambda@f$].  The maximum likelihood estimators
    * are the values @f$(\hat{\alpha}, \hat{\lambda})@f$ that satisfy
    * the equations:
    * @f{align*}{
    *    \hat{\lambda} 
    *    & 
    *    = 
    *    \bar{x}_n - \frac{\sum_{i=1}^n x_i  e^{- \hat{\lambda} x_i}}{\sum_{i=1}^n e^{-\hat{\lambda} x_i}}
    *    \\ 
    *    \hat{\alpha} 
    *    & 
    *    = 
    *    - \frac{1}{\hat{\lambda}} \ln\left( \frac{1}{n} \sum_{i=1}^n e^{-\hat{\lambda} x_i} \right),
    * @f}
    * where @f$\bar{x}_n@f$ is the average of @f$x[0],…,x[n-1]@f$
    * @cite tEVA00a&thinsp; (page 89).
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{\alpha}@f$,
    * @f$\hat{\lambda}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double parameters[] = new double[2];

      double sum = 0.0;
      for (int i = 0; i < n; i++)
         sum += x[i];
      double mean = sum / (double) n;

      sum = 0.0;
      for (int i = 0; i < n; i++)
         sum += (x[i] - mean) * (x[i] - mean);
      double variance = sum / ((double) n - 1.0);

      double lambda0 = Math.PI / Math.sqrt (6 * variance);

      Function f = new Function (x, n, mean);

      double a;
      if ((a = lambda0 - 10.0) < 0)
         a = 1e-15;
      parameters[1] = RootFinder.brentDekker (a, lambda0 + 10.0, f, 1e-7);

      double sumExp = 0.0;
      for (int i = 0; i < n; i++)
         sumExp += Math.exp (- x[i] * parameters[1]);
      parameters[0] = - Math.log (sumExp / (double) n) / parameters[1];

      return parameters;
   }

   /**
    * Same as  #getMLE.
    */
   @Deprecated
   public static double[] getMaximumLikelihoodEstimate (double[] x, int n) {
      return getMLE(x, n);
   }

   /**
    * Creates a new instance of an extreme value distribution with
    * parameters @f$\alpha@f$ and @f$\lambda@f$ estimated using the
    * maximum likelihood method based on the @f$n@f$ observations
    * @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static ExtremeValueDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new ExtremeValueDist (parameters[0], parameters[1]);
   }

   /**
    * Computes and returns the mean, @f$E[X] = \alpha+
    * \gamma/\lambda@f$, of the extreme value distribution with
    * parameters @f$\alpha@f$ and @f$\lambda@f$, where @f$\gamma=
    * 0.5772156649@f$ is the Euler-Mascheroni constant.
    *  @return the mean of the Extreme Value distribution @f$E[X] =
    * \alpha+ \gamma/ \lambda@f$
    */
   public static double getMean (double alpha, double lambda) {
     if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return (alpha + Num.EULER / lambda);
   }

   /**
    * Computes and returns the variance, @f$\mbox{Var}[X] =
    * \pi^2/(6\lambda^2)@f$, of the extreme value distribution with
    * parameters @f$\alpha@f$ and @f$\lambda@f$.
    *  @return the variance of the extreme value distribution
    * @f$\mbox{Var}[X] = 1/6 \pi^2 1/\lambda^2@f$
    */
   public static double getVariance (double alpha, double lambda) {
     if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return ((1.0 / 6.0 * Math.PI * Math.PI) * (1.0 / (lambda * lambda)));
   }

   /**
    * Computes and returns the standard deviation of the extreme value
    * distribution with parameters @f$\alpha@f$ and @f$\lambda@f$.
    *  @return the standard deviation of the extreme value distribution
    */
   public static double getStandardDeviation (double alpha, double lambda) {
     if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return (Math.sqrt(1.0 / 6.0) * Math.PI / lambda);
   }

   /**
    * Returns the parameter @f$\alpha@f$ of this object.
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