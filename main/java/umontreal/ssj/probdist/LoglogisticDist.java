/*
 * Class:        LoglogisticDist
 * Description:  log-logistic distribution
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
import umontreal.ssj.util.Misc;
import optimization.*;

/**
 * Extends the class  @ref ContinuousDistribution for the *Log-Logistic*
 * distribution with shape parameter @f$\alpha> 0@f$ and scale parameter
 * @f$\beta> 0@f$. Its density is
 * @anchor REF_probdist_LoglogisticDist_eq_floglogistic
 * @f[
 *   f(x) = \frac{\alpha(x / \beta)^{\alpha- 1}}{\beta[1 + (x / \beta)^{\alpha}]^2} \qquad\qquad\mbox{for } x > 0 \tag{floglogistic}
 * @f]
 * and its distribution function is
 * @anchor REF_probdist_LoglogisticDist_eq_Floglogistic
 * @f[
 *   F(x) = \frac{1}{1 + (\frac{x}{\beta})^{-\alpha}} \qquad\qquad\mbox{for } x > 0. \tag{Floglogistic}
 * @f]
 * The complementary distribution is
 * @anchor REF_probdist_LoglogisticDist_eq_Fbarloglogistic
 * @f[
 *   \bar{F}(x) = \frac{1}{1 + (\frac{x}{\beta})^{\alpha}} \qquad\qquad\mbox{for } x > 0. \tag{Fbarloglogistic}
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class LoglogisticDist extends ContinuousDistribution {
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
         if ((p[1] <= 0.0) || (p[2] <= 0.0))
            return 1e200;

         double sum = 0.0;
         for (int i = 0; i < n; i++) {
	    double tmp = density (p[1], p[2], xi[i]);
            if (tmp > 0.0)
	        sum -= Math.log (tmp);
            else
	        sum += 709.0;    // log (Double.MIN_VALUE)
	 }
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
    * Constructs a log-logistic distribution with parameters @f$\alpha@f$
    * and @f$\beta@f$.
    */
   public LoglogisticDist (double alpha, double beta) {
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

   public double inverseF (double u) {
      return inverseF (alpha, beta, u);
   }

   public double getMean() {
      return getMean (alpha, beta);
   }

   public double getVariance() {
      return getVariance (alpha, beta);
   }

   public double getStandardDeviation() {
      return getStandardDeviation (alpha, beta);
   }

   /**
    * Computes the density function (
    * {@link REF_probdist_LoglogisticDist_eq_floglogistic
    * floglogistic} ) for a log-logisitic distribution with parameters
    * @f$\alpha@f$ and @f$\beta@f$.
    */
   public static double density (double alpha, double beta, double x) {
      double denominateur;

      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (x <= 0.0 || x >= Double.MAX_VALUE / 2.0)
         return 0.0;

      if (x <= beta) {
         double v = Math.pow (x / beta, alpha);
         denominateur = 1.0 + v;
         denominateur *= denominateur * beta;
         return alpha * v * beta / (x * denominateur);
      } else {
         double v = Math.pow (beta / x, alpha);
         denominateur = 1.0 + v;
         denominateur *= denominateur * beta;
         return alpha * v * beta / (x * denominateur);
      }
   }

   /**
    * Computes the distribution function (
    * {@link REF_probdist_LoglogisticDist_eq_Floglogistic
    * Floglogistic} ) of the log-logistic distribution with parameters
    * @f$\alpha@f$ and @f$\beta@f$.
    */
   public static double cdf (double alpha, double beta, double x) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (x <= 0.0)
         return 0.0;
      if (x >= Double.MAX_VALUE / 2.0)
         return 1.0;
      double z = x/beta;
      if (z >= 1.0)
         return 1.0 / (1.0 + Math.pow (1.0/z, alpha));
      double v = Math.pow (z, alpha);
      return v/(v + 1.0);
   }

   /**
    * Computes the complementary distribution function (
    * {@link REF_probdist_LoglogisticDist_eq_Fbarloglogistic
    * Fbarloglogistic} ) of the log-logistic distribution with parameters
    * @f$\alpha@f$ and @f$\beta@f$.
    */
   public static double barF (double alpha, double beta, double x) {
      double power;

      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (x <= 0.0)
         return 1.0;
      if (x >= Double.MAX_VALUE / 2.0)
         return 0.0;

      double z = x/beta;
      if (z <= 1.0)
         return 1.0 / (1.0 + Math.pow (z, alpha));
      double v = Math.pow (1.0/z, alpha);
      return v/(v + 1.0);
   }

   /**
    * Computes the inverse of the log-logistic distribution with
    * parameters @f$\alpha@f$ and @f$\beta@f$.
    */
   public static double inverseF (double alpha, double beta, double u) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u not in (0, 1]");
      if (u >= 1.0)
         return Double.POSITIVE_INFINITY;
      if (u <= 0.0)
         return 0.0;

      if (u <= 0.5)
         return (beta * Math.pow (u / (1.0 - u), 1.0 / alpha));
      else
         return (beta / Math.pow ((1.0 - u)/ u, 1.0 / alpha));
   }

   /**
    * Estimates the parameters @f$(\alpha,\beta)@f$ of the log-logistic
    * distribution using the maximum likelihood method, from the @f$n@f$
    * observations @f$x[i]@f$, @f$i = 0, 1,…, n-1@f$. The estimates are
    * returned in a two-element array, in regular order: [@f$\alpha@f$,
    * @f$\beta@f$].  The estimate of the parameters is given by
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

      param[2] = EmpiricalDist.getMedian (x, n);

      if (param[2] < 0) throw new IllegalArgumentException ("median < 0");
      if (param[2] <= 0) param[2] = 1.0;

      int m = Math.round ((float) n / 4.0f);
      double q1 = Misc.quickSelect (x, n, m);

      if (q1 < 0) throw new IllegalArgumentException ("x[i] < 0");
      if (q1 > 0)
          param[1] = Math.log (3) / (Math.log(param[2]) - Math.log(q1));
      else
          param[1] = 1.0;

      Uncmin_f77.optif0_f77 (2, param, system, xpls, fpls, gpls, itrcmd, a, udiag);

      for (int i = 0; i < 2; i++)
         parameters[i] = xpls[i+1];

      return parameters;
   }

   /**
    * Creates a new instance of a log-logistic distribution with
    * parameters @f$\alpha@f$ and @f$\beta@f$ estimated using the
    * maximum likelihood method based on the @f$n@f$ observations
    * @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static LoglogisticDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new LoglogisticDist (parameters[0], parameters[1]);
   }

   /**
    * Computes and returns the mean
    *  @f$E[X] = \beta\theta \mbox{cosec}(\theta), \mbox{ where }
    * \theta= \pi/\alpha,@f$
    *  of the log-logistic distribution with parameters @f$\alpha@f$ and
    * @f$\beta@f$.
    *  @return the mean of the log-logistic distribution @f$E[X] =
    * \beta\theta \mbox{cosec}(\theta), \mbox{ where } \theta= \pi/
    * \alpha@f$
    */
   public static double getMean (double alpha, double beta) {
      double theta;

      if (alpha <= 1.0)
         throw new IllegalArgumentException ("alpha <= 1");
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");

      theta = Math.PI / alpha;

      return (beta * theta / Math.sin (theta));
   }

   /**
    * Computes and returns the variance
    *  @f$\mbox{Var}[X] = \beta^2 \theta(2 \mbox{cosec}(2 \theta) -
    * \theta[\mbox{cosec}(\theta)]^2), \mbox{ where } \theta=
    * {\pi/\alpha},@f$
    *  of the log-logistic distribution with parameters @f$\alpha@f$ and
    * @f$\beta@f$.
    *  @return the variance of the log-logistic distribution
    * @f$\mbox{Var}[X] = \beta^2 \theta(2 \mbox{cosec}(2 \theta) -
    * \theta[\mbox{cosec}(\theta)]^2), \mbox{ where } \theta= \pi/
    * \alpha@f$
    */
   public static double getVariance (double alpha, double beta) {
      double theta;

      if (alpha <= 2.0)
         throw new IllegalArgumentException ("alpha <= 2");
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");

      theta = Math.PI / alpha;

      return (beta * beta * theta * ((2.0 / Math.sin (2.0 * theta)) -
                   (theta / (Math.sin (theta) * Math.sin (theta)))));
   }

   /**
    * Computes and returns the standard deviation of the log-logistic
    * distribution with parameters @f$\alpha@f$ and @f$\beta@f$.
    *  @return the standard deviation of the log-logistic distribution
    */
   public static double getStandardDeviation (double alpha, double beta) {
      return Math.sqrt (getVariance (alpha, beta));
   }

   /**
    * Return the parameter @f$\alpha@f$ of this object.
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
    * Sets the parameters @f$\alpha@f$ and @f$\beta@f$ of this object.
    */
   public void setParams (double alpha, double beta) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");

      this.alpha  = alpha;
      this.beta = beta;
      supportA = 0.0;
   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$\alpha@f$,
    * @f$\beta@f$].
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