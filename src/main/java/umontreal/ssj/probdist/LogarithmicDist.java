/*
 * Class:        LogarithmicDist
 * Description:  logarithmic distribution
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
 * Extends the class  @ref DiscreteDistributionInt for the *logarithmic*
 * distribution. It has shape parameter @f$\theta@f$, where @f$0 <
 * \theta<1@f$. Its mass function is
 * @anchor REF_probdist_LogarithmicDist_eq_flogar
 * @f[
 *   p(x) = \frac{-\theta^x}{x\log(1- \theta)}  \qquad\mbox{for } x = 1,2,3,…\tag{flogar}
 * @f]
 * Its distribution function is
 * @f[
 *   F(x) = \frac{-1}{\log(1 - \theta)}\sum_{i=1}^x \frac{\theta^i}{i}, \qquad\mbox{ for } x = 1, 2, 3, …
 * @f]
 * and is 0 for @f$ x\le0@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_discrete
 */
public class LogarithmicDist extends DiscreteDistributionInt {

   private double theta;
   private double t;

   private static class Function implements MathFunction {
      protected double mean;

      public Function (double mean) {
         this.mean = mean;
      }

      public double evaluate (double x) {
         if (x <= 0.0 || x >= 1.0) return 1.0e200;
         return (x + mean * (1.0 - x) * Math.log1p (-x));
      }
   }

   /**
    * Constructs a logarithmic distribution with parameter @f$\theta= @f$
    * `theta`.
    */
   public LogarithmicDist (double theta) {
      setTheta (theta);
   }


   public double prob (int x) {
      if (x < 1)
         return 0;
      return t*Math.pow (theta, x)/x;
   }

   public double cdf (int x) {
      if (x < 1)
         return 0;
      double res = prob (1);
      double term = res;
      for (int i = 2; i <= x; i++) {
         term *= theta;
         res += term/i;
      }
      return res;
   }

   public double barF (int x) {
      if (x <= 1)
         return 1.0;
      double res = prob (x);
      double term = res;
      int i = x + 1;
      while (term > EPSILON) {
         term *= theta*(i-1)/i;
         res += term;
      }
      return res;
   }

   public int inverseFInt (double u) {
      return inverseF (theta, u);
   }

   public double getMean() {
      return LogarithmicDist.getMean (theta);
   }

   public double getVariance() {
      return LogarithmicDist.getVariance (theta);
   }

   public double getStandardDeviation() {
      return LogarithmicDist.getStandardDeviation (theta);
   }

/**
 * Computes the logarithmic probability @f$p(x)@f$ given in (
 * {@link REF_probdist_LogarithmicDist_eq_flogar flogar} ) .
 */
public static double prob (double theta, int x) {
      if (theta <= 0 || theta >= 1)
         throw new IllegalArgumentException ("theta not in range (0,1)");
      if (x < 1)
         return 0;
      return -1.0/Math.log1p(-theta) * Math.pow (theta, x)/x;
   }

   /**
    * Computes the distribution function @f$F(x)@f$.
    */
   public static double cdf (double theta, int x) {
      if (theta <= 0 || theta >= 1)
         throw new IllegalArgumentException ("theta not in range (0,1)");
      if (x < 1)
         return 0;
      double res = prob (theta, 1);
      double term = res;
      for (int i = 2; i <= x; i++) {
         term *= theta;
         res += term/i;
      }
      return res;
   }

   /**
    * Computes the complementary distribution function. *WARNING:* The
    * complementary distribution function is defined as @f$\bar{F}(x) =
    * P[X \ge x]@f$.
    */
   public static double barF (double theta, int x) {
      if (theta <= 0 || theta >= 1)
         throw new IllegalArgumentException ("theta not in range (0,1)");
      if (x <= 1)
         return 1.0;
      double res = prob (theta, x);
      double term = res;
      int i = x + 1;
      while (term > EPSILON) {
         term *= theta*(i-1)/i;
         res += term;
      }
      return res;
   }


   public static int inverseF (double theta, double u) {
      throw new UnsupportedOperationException();
   }

/**
 * Estimates the parameter @f$\theta@f$ of the logarithmic distribution
 * using the maximum likelihood method, from the @f$n@f$ observations
 * @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$. The estimate is returned in element 0
 * of the returned array.  The maximum likelihood estimator
 * @f$\hat{\theta}@f$ satisfies the equation (see @cite mEVA00a&thinsp;
 * (page 122))
 * @f{align*}{
 *    \bar{x}_n = \frac{-\hat{\theta}}{(1 - \hat{\theta}) \ln(1 - \hat{\theta})}
 * @f}
 * where @f$\bar{x}_n@f$ is the average of @f$x[0], …, x[n-1]@f$.
 *  @param x            the list of observations used to evaluate parameters
 *  @param n            the number of observations used to evaluate
 *                      parameters
 *  @return returns the parameter [@f$\hat{\theta}@f$]
 */
public static double[] getMLE (int[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double parameters[];
      parameters = new double[1];
      double sum = 0.0;
      for (int i = 0; i < n; i++) {
         sum += x[i];
      }

      double mean = (double) sum / (double) n;

      Function f = new Function (mean);
      parameters[0] = RootFinder.brentDekker (1e-15, 1.0-1e-15, f, 1e-7);

      return parameters;
   }

   /**
    * Creates a new instance of a logarithmic distribution with parameter
    * @f$\theta@f$ estimated using the maximum likelihood method based on
    * the @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static LogarithmicDist getInstanceFromMLE (int[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new LogarithmicDist (parameters[0]);
   }

   /**
    * Computes and returns the mean
    * @f[
    *   E[X] = \frac{-\theta}{(1 - \theta)\ln(1 - \theta)}
    * @f]
    * of the logarithmic distribution with parameter @f$\theta= @f$
    * `theta`.
    *  @return the mean of the logarithmic distribution @f$E[X]
    * = -\theta/ ((1 - \theta) ln(1 - \theta))@f$
    */
   public static double getMean (double theta) {
      if (theta <= 0.0 || theta >= 1.0)
         throw new IllegalArgumentException ("theta not in range (0,1)");

      return ((-1 / Math.log1p(-theta)) * (theta / (1 - theta)));
   }

   /**
    * Computes and returns the variance
    * @f[
    *   \mbox{Var}[X] = \frac{-\theta(\theta+ \ln(1 - \theta))}{[(1 - \theta) \ln(1 - \theta)]^2}
    * @f]
    * of the logarithmic distribution with parameter @f$\theta=@f$
    * `theta`.
    *  @return the variance of the logarithmic distribution
    * @f$\mbox{Var}[X] = -\theta(\theta+ ln(1 - \theta)) / ((1 -
    * \theta)^2 (ln(1 - \theta))^2)@f$
    */
   public static double getVariance (double theta) {
      if (theta <= 0.0 || theta >= 1.0)
         throw new IllegalArgumentException ("theta not in range (0,1)");

      double v = Math.log1p(-theta);
      return ((-theta * (theta + v)) / ((1 - theta) * (1 - theta) * v * v));
   }

   /**
    * Computes and returns the standard deviation of the logarithmic
    * distribution with parameter @f$\theta= @f$ `theta`.
    *  @return the standard deviation of the logarithmic distribution
    */
   public static double getStandardDeviation (double theta) {
      return Math.sqrt (LogarithmicDist.getVariance (theta));
   }

   /**
    * Returns the @f$\theta@f$ associated with this object.
    */
   public double getTheta() {
      return theta;
   }

   /**
    * Sets the @f$\theta@f$ associated with this object.
    */
   public void setTheta (double theta) {
      if (theta <= 0 || theta >= 1)
         throw new IllegalArgumentException ("theta not in range (0,1)");
      this.theta = theta;
      t = -1.0/Math.log1p (-theta);
      supportA = 1;
   }

   /**
    * Return a table containing the parameters of the current
    * distribution.
    */
   public double[] getParams () {
      double[] retour = {theta};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : theta = " + theta;
   }

}