/*
 * Class:        ParetoDist
 * Description:  Pareto distribution
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
 * Extends the class  @ref ContinuousDistribution for a distribution from the
 * *Pareto* family, with shape parameter @f$\alpha> 0@f$ and location
 * parameter @f$\beta> 0@f$ @cite tJOH95a&thinsp; (page 574). The density
 * for this type of Pareto distribution is
 * @anchor REF_probdist_ParetoDist_eq_fpareto
 * @f[
 *   f(x) = \frac{\alpha\beta^{\alpha}}{x^{\alpha+1}} \qquad\mbox{for }x \ge\beta, \tag{fpareto}
 * @f]
 * and 0 otherwise. The distribution function is
 * @anchor REF_probdist_ParetoDist_eq_Fpareto
 * @f[
 *   F(x) = 1 - \left(\beta/x\right)^{\alpha}\qquad\mbox{for }x\ge\beta, \tag{Fpareto}
 * @f]
 * and the inverse distribution function is
 * @f[
 *   F^{-1}(u) = \beta(1 - u)^{-1/\alpha} \qquad\mbox{for } 0 \le u < 1.
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class ParetoDist extends ContinuousDistribution {
   private double alpha;
   private double beta;

   /**
    * Constructs a `ParetoDist` object with parameters @f$\alpha=@f$
    * `alpha` and @f$\beta= 1@f$.
    */
   public ParetoDist (double alpha) {
      setParams (alpha, 1.0);
   }

   /**
    * Constructs a `ParetoDist` object with parameters @f$\alpha=@f$
    * `alpha` and @f$\beta= @f$ `beta`.
    */
   public ParetoDist (double alpha, double beta) {
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
      return ParetoDist.getMean (alpha, beta);
   }

   public double getVariance() {
      return ParetoDist.getVariance (alpha, beta);
   }

   public double getStandardDeviation() {
      return ParetoDist.getStandardDeviation (alpha, beta);
   }

/**
 * Computes the density function.
 */
public static double density (double alpha, double beta, double x) {
      if (alpha <= 0.0)
        throw new IllegalArgumentException ("alpha <= 0");
      if (beta <= 0.0)
        throw new IllegalArgumentException ("beta <= 0");

      return x < beta ? 0 : alpha*Math.pow (beta/x, alpha)/x;
   }

   /**
    * Computes the distribution function.
    */
   public static double cdf (double alpha, double beta, double x) {
      if (alpha <= 0.0)
        throw new IllegalArgumentException ("alpha <= 0");
      if (beta <= 0.0)
        throw new IllegalArgumentException ("beta <= 0");
      if (x <= beta)
         return 0.0;
      return 1.0 - Math.pow (beta/x, alpha);
   }

   /**
    * Computes the complementary distribution function.
    */
   public static double barF (double alpha, double beta, double x) {
      if (alpha <= 0)
        throw new IllegalArgumentException ("c <= 0");
      if (beta <= 0.0)
        throw new IllegalArgumentException ("beta <= 0");
      if (x <= beta)
         return 1.0;
      return Math.pow (beta/x, alpha);
   }

   /**
    * Computes the inverse of the distribution function.
    */
   public static double inverseF (double alpha, double beta, double u) {
      if (alpha <= 0)
        throw new IllegalArgumentException ("c <= 0");
      if (beta <= 0.0)
        throw new IllegalArgumentException ("beta <= 0");

      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u not in [0,1]");

      if (u <= 0.0)
         return beta;

      double t;
      t = -Math.log1p (-u);
      if ((u >= 1.0) || t/Math.log(10) >= alpha * Num.DBL_MAX_10_EXP)
         return Double.POSITIVE_INFINITY;

      return beta / Math.pow (1 - u, 1.0/alpha);
   }

   /**
    * Estimates the parameters @f$(\alpha,\beta)@f$ of the Pareto
    * distribution using the maximum likelihood method, from the @f$n@f$
    * observations @f$x[i]@f$, @f$i = 0, 1,…, n-1@f$. The estimates are
    * returned in a two-element array, in regular order: [@f$\alpha@f$,
    * @f$\beta@f$].  The maximum likelihood estimators are the values
    * @f$(\hat{\alpha}, \hat{\beta})@f$ that satisfy the equations:
    * @f{align*}{
    *    \hat{\beta} 
    *    & 
    *    = 
    *    \min_i \{x_i\}
    *    \\ 
    *    \hat{\alpha} 
    *    & 
    *    = 
    *    \frac{n}{\sum_{i=1}^n \ln\left(\frac{x_i}{\hat{\beta}\Rule{0.0pt}{5.5pt}{0.0pt}}\right)}.
    * @f}
    * See @cite tJOH95a&thinsp; (page 581).
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{\alpha}@f$,
    * @f$\hat{\beta}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double [] parameters = new double[2];
      parameters[1] = Double.POSITIVE_INFINITY;
      for (int i = 0; i < n; i++) {
         if (x[i] < parameters[1])
            parameters[1] = x[i];
      }

      double sum = 0.0;
      for (int i = 0; i < n; i++) {
         if (x[i] > 0.0)
            sum += Math.log (x[i] / parameters[1]);
         else
            sum -= 709.0;
      }
      parameters[0] = n / sum;
      return parameters;
   }

   /**
    * Creates a new instance of a Pareto distribution with parameters
    * @f$\alpha@f$ and @f$\beta@f$ estimated using the maximum
    * likelihood method based on the @f$n@f$ observations @f$x[i]@f$, @f$i
    * = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static ParetoDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new ParetoDist (parameters[0], parameters[1]);
   }

   /**
    * Computes and returns the mean @f$E[X] = \alpha\beta/(\alpha-
    * 1)@f$ of the Pareto distribution with parameters @f$\alpha@f$ and
    * @f$\beta@f$.
    *  @return the mean of the Pareto distribution
    */
   public static double getMean (double alpha, double beta) {
      if (alpha <= 1.0)
         throw new IllegalArgumentException("alpha <= 1");
      if (beta <= 0.0)
        throw new IllegalArgumentException ("beta <= 0");

      return ((alpha * beta) / (alpha - 1.0));
   }

   /**
    * Computes and returns the variance
    *  @f$\mbox{Var}[X] = \frac{\alpha\beta^2}{(\alpha- 2)(\alpha-
    * 1)}@f$
    *  of the Pareto distribution with parameters @f$\alpha@f$ and
    * @f$\beta@f$.
    *  @return the variance of the Pareto distribution @f$\mbox{Var}[X] =
    * \alpha\beta^2 / [(\alpha- 2)(\alpha- 1)]@f$
    */
   public static double getVariance (double alpha, double beta) {
      if (alpha <= 2)
         throw new IllegalArgumentException("alpha <= 2");
      if (beta <= 0.0)
        throw new IllegalArgumentException ("beta <= 0");

      return ((alpha * beta * beta) / ((alpha - 2.0) * (alpha - 1.0) * (alpha - 1.0)));
   }

   /**
    * Computes and returns the standard deviation of the Pareto
    * distribution with parameters @f$\alpha@f$ and @f$\beta@f$.
    *  @return the standard deviation of the Pareto distribution
    */
   public static double getStandardDeviation (double alpha, double beta) {
      return Math.sqrt (ParetoDist.getVariance (alpha, beta));
   }

   /**
    * Returns the parameter @f$\alpha@f$.
    */
   public double getAlpha() {
      return alpha;
   }

   /**
    * Returns the parameter @f$\beta@f$.
    */
   public double getBeta() {
      return beta;
   }

   /**
    * Sets the parameter @f$\alpha@f$ and @f$\beta@f$ for this object.
    */
   public void setParams (double alpha, double beta) {
      if (alpha <= 0.0)
        throw new IllegalArgumentException ("alpha <= 0");
      if (beta <= 0.0)
        throw new IllegalArgumentException ("beta <= 0");

      this.alpha = alpha;
      this.beta = beta;
      supportA = beta;
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