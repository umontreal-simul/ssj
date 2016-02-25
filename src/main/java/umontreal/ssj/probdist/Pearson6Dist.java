/*
 * Class:        Pearson6Dist
 * Description:  Pearson type VI distribution
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
import optimization.*;

/**
 * Extends the class  @ref ContinuousDistribution for the *Pearson type VI*
 * distribution with shape parameters @f$\alpha_1 > 0@f$ and @f$\alpha_2 >
 * 0@f$, and scale parameter @f$\beta> 0@f$. The density function is given
 * by
 * @anchor REF_probdist_Pearson6Dist_eq_fpearson6
 * @f[
 *   f(x) =\left\{\begin{array}{ll}
 *    \displaystyle\frac{\left(x/{\beta}\right)^{\alpha_1 - 1}}{\beta\mathcal{B}(\alpha_1, \alpha_2)(1 + x/{\beta})^{\alpha_1 + \alpha_2}} 
 *    & 
 *    \quad\mbox{for } x > 0, 
 *    \\ 
 *    0 
 *    & 
 *    \quad\mbox{otherwise,} 
 *   \end{array} \right. \tag{fpearson6}
 * @f]
 * where @f$\mathcal{B}@f$ is the beta function. The distribution function is
 * given by
 * @anchor REF_probdist_Pearson6Dist_eq_Fpearson6
 * @f[
 *   F(x) = F_B\left(\frac{x}{x + \beta}\right) \qquad\mbox{for } x > 0, \tag{Fpearson6}
 * @f]
 * and @f$F(x) = 0@f$ otherwise, where @f$F_B(x)@f$ is the distribution
 * function of a beta distribution with shape parameters @f$\alpha_1@f$ and
 * @f$\alpha_2@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class Pearson6Dist extends ContinuousDistribution {
   protected double alpha1;
   protected double alpha2;
   protected double beta;
   protected double logBeta; // Ln (Beta (alpha1, alpha2))

   private static class Optim implements Uncmin_methods {
      private int n;
      private double[] x;

      public Optim (double[] x, int n) {
         this.n = n;
         this.x = new double[n];
         System.arraycopy (x, 0, this.x, 0, n);
      }

      public double f_to_minimize (double[] param) {

         if ((param[1] <= 0.0) || (param[2] <= 0.0) || (param[3] <= 0.0))
            return 1e200;

         double sumLogY = 0.0;
         double sumLog1_Y = 0.0;
         for (int i = 0; i < n; i++)
         {
            if (x[i] > 0.0)
               sumLogY += Math.log (x[i] / param[3]);
            else
               sumLogY -= 709.0;
            sumLog1_Y += Math.log1p (x[i] / param[3]);
         }

         return (n * (Math.log (param[3]) + Num.lnBeta (param[1], param[2])) -
         (param[1] - 1.0) * sumLogY + (param[1] + param[2]) * sumLog1_Y);
      }

      public void gradient (double[] x, double[] g)
      {
      }

      public void hessian (double[] x, double[][] h)
      {
      }
   }

   /**
    * Constructs a `Pearson6Dist` object with parameters @f$\alpha_1@f$ =
    * `alpha1`, @f$\alpha_2@f$ = `alpha2` and @f$\beta@f$ = `beta`.
    */
   public Pearson6Dist (double alpha1, double alpha2, double beta) {
      setParam (alpha1, alpha2, beta);
   }


   public double density (double x) {
      if (x <= 0.0)
         return 0.0;
      return Math.exp ((alpha1 - 1.0) * Math.log (x / beta) - (logBeta +
            (alpha1 + alpha2) * Math.log1p (x / beta))) / beta;
   }

   public double cdf (double x) {
      return cdf (alpha1, alpha2, beta, x);
   }

   public double barF (double x) {
      return barF (alpha1, alpha2, beta, x);
   }

   public double inverseF (double u) {
      return inverseF (alpha1, alpha2, beta, u);
   }

   public double getMean () {
      return getMean (alpha1, alpha2, beta);
   }

   public double getVariance () {
      return getVariance (alpha1, alpha2, beta);
   }

   public double getStandardDeviation () {
      return getStandardDeviation (alpha1, alpha2, beta);
   }

/**
 * Computes the density function of a Pearson VI distribution with shape
 * parameters @f$\alpha_1@f$ and @f$\alpha_2@f$, and scale parameter
 * @f$\beta@f$.
 */
public static double density (double alpha1, double alpha2,
                                 double beta, double x) {
      if (alpha1 <= 0.0)
         throw new IllegalArgumentException("alpha1 <= 0");
      if (alpha2 <= 0.0)
         throw new IllegalArgumentException("alpha2 <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException("beta <= 0");
      if (x <= 0.0)
         return 0.0;

      return Math.exp ((alpha1 - 1.0) * Math.log (x / beta) -
         (Num.lnBeta (alpha1, alpha2) + (alpha1 + alpha2) * Math.log1p (x / beta))) / beta;
   }

   /**
    * Computes the distribution function of a Pearson VI distribution with
    * shape parameters @f$\alpha_1@f$ and @f$\alpha_2@f$, and scale
    * parameter @f$\beta@f$.
    */
   public static double cdf (double alpha1, double alpha2,
                             double beta, double x) {
      if (alpha1 <= 0.0)
         throw new IllegalArgumentException("alpha1 <= 0");
      if (alpha2 <= 0.0)
         throw new IllegalArgumentException("alpha2 <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException("beta <= 0");
      if (x <= 0.0)
         return 0.0;

      return BetaDist.cdf (alpha1, alpha2, x / (x + beta));
   }

   /**
    * Computes the complementary distribution function of a Pearson VI
    * distribution with shape parameters @f$\alpha_1@f$ and
    * @f$\alpha_2@f$, and scale parameter @f$\beta@f$.
    */
   public static double barF (double alpha1, double alpha2,
                              double beta, double x) {
      if (alpha1 <= 0.0)
         throw new IllegalArgumentException("alpha1 <= 0");
      if (alpha2 <= 0.0)
         throw new IllegalArgumentException("alpha2 <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException("beta <= 0");
      if (x <= 0.0)
         return 1.0;

      return BetaDist.barF (alpha1, alpha2, x / (x + beta));
   }

   /**
    * Computes the inverse distribution function of a Pearson VI
    * distribution with shape parameters @f$\alpha_1@f$ and
    * @f$\alpha_2@f$, and scale parameter @f$\beta@f$.
    */
   public static double inverseF (double alpha1, double alpha2,
                                  double beta, double u) {
      if (alpha1 <= 0.0)
         throw new IllegalArgumentException("alpha1 <= 0");
      if (alpha2 <= 0.0)
         throw new IllegalArgumentException("alpha2 <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException("beta <= 0");

      double y = BetaDist.inverseF (alpha1, alpha2, u);

      return ((y * beta) / (1.0 - y));
   }

   /**
    * Estimates the parameters @f$(\alpha_1,\alpha_2,\beta)@f$ of the
    * Pearson VI distribution using the maximum likelihood method, from
    * the @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1,…, n-1@f$. The
    * estimates are returned in a three-element array, in regular order:
    * [@f$\alpha_1, \alpha_2@f$, @f$\beta@f$].  The estimate of the
    * parameters is given by maximizing numerically the log-likelihood
    * function, using the Uncmin package @cite iSCHa, @cite iVERa&thinsp;.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{\alpha_1},
    * \hat{\alpha_2}, \hat{\beta}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double[] parameters = new double[3];
      double[] xpls = new double[4];
      double[] param = new double[4];
      double[] fpls = new double[4];
      double[] gpls = new double[4];
      int[] itrcmd = new int[2];
      double[][] h = new double[4][4];
      double[] udiag = new double[4];

      Optim system = new Optim (x, n);

      double mean = 0.0;
      double mean2 = 0.0;
      double mean3 = 0.0;
      for (int i = 0; i < n; i++)
      {
         mean += x[i];
         mean2 += x[i] * x[i];
         mean3 += x[i] * x[i] * x[i];
      }
      mean /= (double) n;
      mean2 /= (double) n;
      mean3 /= (double) n;

      double r1 = mean2 / (mean * mean);
      double r2 = mean2 * mean / mean3;

      param[1] = - (2.0 * (-1.0 + r1 * r2)) / (-2.0 + r1 + r1 * r2);
      if(param[1] <= 0) param[1] = 1;
      param[2] = (- 3.0 - r2 + 4.0 * r1 * r2) / (- 1.0 - r2 + 2.0 * r1 * r2);
      if(param[2] <= 0) param[2] = 1;
      param[3] = (param[2] - 1.0) * mean / param[1];
      if(param[3] <= 0) param[3] = 1;

      Uncmin_f77.optif0_f77 (3, param, system, xpls, fpls, gpls, itrcmd, h, udiag);

      for (int i = 0; i < 3; i++)
         parameters[i] = xpls[i+1];

      return parameters;
   }

   /**
    * Creates a new instance of a Pearson VI distribution with parameters
    * @f$\alpha_1@f$, @f$\alpha_2@f$ and @f$\beta@f$, estimated using
    * the maximum likelihood method based on the @f$n@f$ observations
    * @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static Pearson6Dist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new Pearson6Dist (parameters[0], parameters[1], parameters[2]);
   }

   /**
    * Computes and returns the mean @f$E[X] = (\beta\alpha_1) /
    * (\alpha_2 - 1)@f$ of a Pearson VI distribution with shape
    * parameters @f$\alpha_1@f$ and @f$\alpha_2@f$, and scale parameter
    * @f$\beta@f$.
    */
   public static double getMean (double alpha1, double alpha2,
                                 double beta) {
      if (alpha1 <= 0.0)
         throw new IllegalArgumentException("alpha1 <= 0");
      if (alpha2 <= 1.0)
         throw new IllegalArgumentException("alpha2 <= 1");
      if (beta <= 0.0)
         throw new IllegalArgumentException("beta <= 0");

      return ((beta * alpha1) / (alpha2 - 1.0));
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] = [\beta^2
    * \alpha_1 (\alpha_1 + \alpha_2 - 1)] / [(\alpha_2 -
    * 1)^2(\alpha_2 - 2)]@f$ of a Pearson VI distribution with shape
    * parameters @f$\alpha_1@f$ and @f$\alpha_2@f$, and scale parameter
    * @f$\beta@f$.
    */
   public static double getVariance (double alpha1, double alpha2,
                                     double beta) {
      if (alpha1 <= 0.0)
         throw new IllegalArgumentException("alpha1 <= 0");
      if (alpha2 <= 0.0)
         throw new IllegalArgumentException("alpha2 <= 2");
      if (beta <= 0.0)
         throw new IllegalArgumentException("beta <= 0");

      return (((beta * beta) * alpha1 * (alpha1 + alpha2 - 1.0)) /
((alpha2 - 1.0) * (alpha2 - 1.0) * (alpha2 - 2.0)));
   }

   /**
    * Computes and returns the standard deviation of a Pearson VI
    * distribution with shape parameters @f$\alpha_1@f$ and
    * @f$\alpha_2@f$, and scale parameter @f$\beta@f$.
    */
   public static double getStandardDeviation (double alpha1, double alpha2,
                                              double beta) {
      return Math.sqrt (getVariance (alpha1, alpha2, beta));
   }

   /**
    * Returns the @f$\alpha_1@f$ parameter of this object.
    */
   public double getAlpha1() {
      return alpha1;
   }

   /**
    * Returns the @f$\alpha_2@f$ parameter of this object.
    */
   public double getAlpha2() {
      return alpha2;
   }

   /**
    * Returns the @f$\beta@f$ parameter of this object.
    */
   public double getBeta() {
      return beta;
   }

   /**
    * Sets the parameters @f$\alpha_1@f$, @f$\alpha_2@f$ and
    * @f$\beta@f$ of this object.
    */
   public void setParam (double alpha1, double alpha2, double beta) {
      if (alpha1 <= 0.0)
         throw new IllegalArgumentException("alpha1 <= 0");
      if (alpha2 <= 0.0)
         throw new IllegalArgumentException("alpha2 <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException("beta <= 0");
      supportA = 0.0;
      this.alpha1 = alpha1;
      this.alpha2 = alpha2;
      this.beta = beta;
      logBeta = Num.lnBeta (alpha1, alpha2);
   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$\alpha_1@f$,
    * @f$\alpha_2@f$, @f$\beta@f$].
    */
   public double[] getParams () {
      double[] retour = {alpha1, alpha2, beta};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : alpha1 = " + alpha1 + ", alpha2 = " + alpha2 + ", beta = " + beta;
   }

}