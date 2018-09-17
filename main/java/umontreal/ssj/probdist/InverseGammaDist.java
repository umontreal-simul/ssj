/*
 * Class:        InverseGammaDist
 * Description:  inverse gamma distribution
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

/**
 * Extends the class  @ref ContinuousDistribution for the *inverse gamma*
 * distribution with shape parameter @f$\alpha> 0@f$ and scale parameter
 * @f$\beta> 0@f$, also known as the *Pearson type V* distribution. 
 * The density function is given by
 * @anchor REF_probdist_InverseGammaDist_eq_dinvgam
 * @f[
 *   f(x) = \left\{\begin{array}{ll}
 *    \displaystyle\frac{\beta^{\alpha}e^{-\beta/ x}}{x^{\alpha+ 1} \Gamma(\alpha)} 
 *    & 
 *    \quad\mbox{for } x > 0 
 *    \\ 
 *    0 
 *    & 
 *    \quad\mbox{otherwise,} 
 *   \end{array} \right. \tag{dinvgam}
 * @f]
 * where @f$\Gamma@f$ is the gamma function. The distribution function is
 * given by
 * @anchor REF_probdist_InverseGammaDist_eq_Finvgam
 * @f[
 *   F(x) = 1 - F_G\left(\frac{1}{x}\right) \qquad\mbox{for } x > 0, \tag{Finvgam}
 * @f]
 * and @f$F(x) = 0@f$ otherwise, where @f$F_G(x)@f$ is the distribution
 * function of a gamma distribution with shape parameter @f$\alpha@f$ and
 * scale parameter @f$\beta@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class InverseGammaDist extends ContinuousDistribution {
   protected double alpha;
   protected double beta;
   protected double logam;   // Ln (Gamma(alpha))

   /**
    * Constructs an `InverseGammaDist` object with parameters
    * @f$\alpha@f$ = `alpha` and @f$\beta@f$ = `beta`.
    */
   public InverseGammaDist (double alpha, double beta) {
      setParam (alpha, beta);
   }


   public double density (double x) {
      if (x <= 0.0)
         return 0.0;
      return Math.exp (alpha * Math.log (beta/x) - (beta / x) - logam) / x;
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

   public double getMean () {
      return getMean (alpha, beta);
   }

   public double getVariance () {
      return getVariance (alpha, beta);
   }

   public double getStandardDeviation () {
      return getStandardDeviation (alpha, beta);
   }

/**
 * Computes the density function of the inverse gamma distribution with shape
 * parameter @f$\alpha@f$ and scale parameter @f$\beta@f$.
 */
public static double density (double alpha, double beta, double x) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException("beta <= 0");
      if (x <= 0.0)
         return 0.0;

      return Math.exp (alpha * Math.log (beta/x) - (beta / x) - Num.lnGamma (alpha)) / x;
   }

   /**
    * Computes the cumulative probability function of the inverse gamma
    * distribution with shape parameter @f$\alpha@f$ and scale parameter
    * @f$\beta@f$.
    */
   public static double cdf (double alpha, double beta, double x) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException("beta <= 0");
      if (x <= 0.0)
         return 0.0;

      return GammaDist.barF (alpha, beta, 15, 1.0 / x);
   }

   /**
    * Computes the complementary distribution function of the inverse
    * gamma distribution with shape parameter @f$\alpha@f$ and scale
    * parameter @f$\beta@f$.
    */
   public static double barF (double alpha, double beta, double x) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException("beta <= 0");
      if (x <= 0.0)
         return 1.0;

      return GammaDist.cdf (alpha, beta, 15, 1.0 / x);
   }

   /**
    * Computes the inverse distribution function of the inverse gamma
    * distribution with shape parameter @f$\alpha@f$ and scale parameter
    * @f$\beta@f$.
    */
   public static double inverseF (double alpha, double beta, double u) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException("beta <= 0");

      return 1.0 / GammaDist.inverseF (alpha, beta, 15, 1 - u);
   }

   /**
    * Estimates the parameters @f$(\alpha,\beta)@f$ of the inverse gamma
    * distribution using the maximum likelihood method, from the @f$n@f$
    * observations @f$x[i]@f$, @f$i = 0, 1,…, n-1@f$. The estimates are
    * returned in a two-element array, in regular order: [@f$\alpha@f$,
    * @f$\beta@f$].  The equations of the maximum likelihood are the same
    * as the equations of the gamma distribution, with the sample @f$y_i =
    * 1/x_i@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{\alpha}, \hat{\beta}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      double[] y = new double[n];

      for (int i = 0; i < n; i++) {
	      if(x[i] > 0)
	         y[i] = 1.0 / x[i];
	      else
	         y[i] = 1.0E100;
      }

      return GammaDist.getMLE (y, n);
   }

   /**
    * Creates a new instance of the inverse gamma distribution with
    * parameters @f$\alpha@f$ and @f$\beta@f$ estimated using the
    * maximum likelihood method based on the @f$n@f$ observations
    * @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static InverseGammaDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new InverseGammaDist (parameters[0], parameters[1]);
   }

   /**
    * Returns the mean @f$E[X] = \beta/ (\alpha- 1)@f$ of the inverse
    * gamma distribution with shape parameter @f$\alpha@f$ and scale
    * parameter @f$\beta@f$.
    */
   public static double getMean (double alpha, double beta) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException("beta <= 0");

      return (beta / (alpha - 1.0));
   }

   /**
    * Returns the variance @f$\mbox{Var}[X] = \beta^2 / ((\alpha-
    * 1)^2(\alpha- 2))@f$ of the inverse gamma distribution with shape
    * parameter @f$\alpha@f$ and scale parameter @f$\beta@f$.
    */
   public static double getVariance (double alpha, double beta) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException("beta <= 0");

      return ((beta * beta) / ((alpha - 1.0) * (alpha - 1.0) * (alpha - 2.0)));
   }

   /**
    * Returns the standard deviation of the inverse gamma distribution
    * with shape parameter @f$\alpha@f$ and scale parameter @f$\beta@f$.
    */
   public static double getStandardDeviation (double alpha, double beta) {
      return Math.sqrt (getVariance (alpha, beta));
   }

   /**
    * Returns the @f$\alpha@f$ parameter of this object.
    */
   public double getAlpha() {
      return alpha;
   }

   /**
    * Returns the @f$\beta@f$ parameter of this object.
    */
   public double getBeta() {
      return beta;
   }

   /**
    * Sets the parameters @f$\alpha@f$ and @f$\beta@f$ of this object.
    */
   public void setParam (double alpha, double beta) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException("beta <= 0");
      supportA = 0.0;
      this.alpha = alpha;
      this.beta = beta;
      logam = Num.lnGamma (alpha);
   }

   /**
    * Returns a table containing the parameters of the current
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