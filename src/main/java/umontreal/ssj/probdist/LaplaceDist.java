/*
 * Class:        LaplaceDist
 * Description:  Laplace distribution
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
 * Extends the class  @ref ContinuousDistribution for the *Laplace*
 * distribution (see, e.g., @cite tJOH95b&thinsp; (page 165)). It has
 * location parameter @f$\mu@f$ and scale parameter @f$\beta> 0@f$. The
 * density function is given by
 * @anchor REF_probdist_LaplaceDist_eq_flaplace
 * @f[
 *   f(x) = \frac{e^{-|x - \mu|/\beta}}{2\beta} \qquad\mbox{ for }-\infty< x < \infty. \tag{flaplace}
 * @f]
 * The distribution function is
 * @f[
 *   F (x) = \left\{\begin{array}{ll}
 *    \frac{1}{2} e^{(x - \mu)/\beta} 
 *    & 
 *    \mbox{ if } x\le\mu, 
 *    \\ 
 *    1 - \frac{1}{2} e^{(\mu- x)/\beta} 
 *    & 
 *    \mbox{ otherwise, } 
 *   \end{array}\right.
 * @f]
 * and its inverse is
 * @f[
 *   F^{-1} (u) = \left\{\begin{array}{ll}
 *    \beta\log(2u) + \mu
 *    & 
 *    \mbox{ if } 0\le u\le\frac{1}{2}, 
 *    \\ 
 *    \mu- \beta\log(2(1-u)) 
 *    & 
 *    \mbox{ otherwise. } 
 *   \end{array}\right.
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class LaplaceDist extends ContinuousDistribution {
   private double mu;
   private double beta;

   /**
    * Constructs a `LaplaceDist` object with default parameters @f$\mu=
    * 0@f$ and @f$\beta= 1@f$.
    */
   public LaplaceDist() {
      mu = 0;
      beta = 1;
   }

   /**
    * Constructs a `LaplaceDist` object with parameters @f$\mu@f$ = `mu`
    * and @f$\beta@f$ = `beta`.
    */
   public LaplaceDist (double mu, double beta) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");

      this.mu = mu;
      this.beta = beta;
   }


   public double density (double x) {
      return density (mu, beta, x);
   }

   public double cdf (double x) {
      return cdf (mu, beta, x);
   }

   public double barF (double x) {
      return barF (mu, beta, x);
   }

   public double inverseF (double u) {
      return inverseF (mu, beta, u);
   }

   public double getMean() {
      return LaplaceDist.getMean (mu, beta);
   }

   public double getVariance() {
      return LaplaceDist.getVariance (mu, beta);
   }

   public double getStandardDeviation() {
      return LaplaceDist.getStandardDeviation (mu, beta);
   }

/**
 * Computes the Laplace density function.
 */
public static double density (double mu, double beta, double x) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      return Math.exp (-Math.abs (x - mu)/beta)/(2.0*beta);
   }

   /**
    * Computes the Laplace distribution function.
    */
   public static double cdf (double mu, double beta, double x) {
      if (x <= mu)
         return Math.exp ((x - mu)/beta)/2.0;
      else
         return 1.0 - Math.exp ((mu - x)/beta)/2.0;
   }

   /**
    * Computes the Laplace complementary distribution function.
    */
   public static double barF (double mu, double beta, double x) {
      if (x <= mu)
         return 1.0 - Math.exp ((x - mu)/beta)/2.0;
      else
         return Math.exp ((mu - x)/beta)/2.0;
   }

  //====================================================
  // code taken and adapted from unuran
  // file /distribution/c_laplca_gen.c
  //====================================================

/**
 * Computes the inverse Laplace distribution function.
 */
public static double inverseF (double mu, double beta, double u) {
     // transform to random variate
     if (u < 0.0 || u > 1.0)
        throw new IllegalArgumentException ("u should be in [0,1]");
     if (u <= 0.0)
        return Double.NEGATIVE_INFINITY;
     if (u >= 1.0)
        return Double.POSITIVE_INFINITY;

     double x = (u>0.5) ? -Math.log(2.-2*u) : Math.log(2*u);
     return mu + beta*x;
   }

   /**
    * Estimates the parameters @f$(\mu, \beta)@f$ of the Laplace
    * distribution using the maximum likelihood method, from the @f$n@f$
    * observations @f$x[i]@f$, @f$i = 0, 1,…, n-1@f$. The estimates are
    * returned in a two-element array, in regular order: [@f$\mu@f$,
    * @f$\beta@f$].  The maximum likelihood estimators are the values
    * @f$(\hat{\mu}, \hat{\beta})@f$ that satisfy the equations:
    * @f{align*}{
    *    \hat{\mu} 
    *    & 
    *    = 
    *    \mbox{the median of } \{x_1,…,x_n\}
    *    \\ 
    *   \hat{\beta} 
    *    & 
    *    = 
    *    \frac{1}{n} \sum_{i=1}^n |x_i - \hat{\mu}|.
    * @f}
    * See @cite tJOH95b&thinsp; (page 172).
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{\mu}@f$,
    * @f$\hat{\beta}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double parameters[];
      parameters = new double[2];

      parameters[0] = EmpiricalDist.getMedian (x, n);

      double sum = 0.0;
      for (int i = 0; i < n; i++)
         sum += Math.abs (x[i] - parameters[0]);
      parameters[1] = sum / (double) n;

      return parameters;
   }

   /**
    * Creates a new instance of a Laplace distribution with parameters
    * @f$\mu@f$ and @f$\beta@f$ estimated using the maximum likelihood
    * method based on the @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1, …,
    * n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static LaplaceDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new LaplaceDist (parameters[0], parameters[1]);
   }

   /**
    * Computes and returns the mean @f$E[X] = \mu@f$ of the Laplace
    * distribution with parameters @f$\mu@f$ and @f$\beta@f$.
    *  @return the mean of the Laplace distribution @f$E[X] = \mu@f$
    */
   public static double getMean (double mu, double beta) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");

      return mu;
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] = 2 \beta^2@f$
    * of the Laplace distribution with parameters @f$\mu@f$ and
    * @f$\beta@f$.
    *  @return the variance of the Laplace distribution @f$\mbox{Var}[X] =
    * 2 \beta^2@f$
    */
   public static double getVariance (double mu, double beta) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");

      return (2.0 * beta * beta);
   }

   /**
    * Computes and returns the standard deviation of the Laplace
    * distribution with parameters @f$\mu@f$ and @f$\beta@f$.
    *  @return the standard deviation of the Laplace distribution
    */
   public static double getStandardDeviation (double mu, double beta) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");

      return (Num.RAC2 * beta);
   }

   /**
    * Returns the parameter @f$\mu@f$.
    */
   public double getMu() {
      return mu;
   }

   /**
    * Returns the parameter @f$\beta@f$.
    */
   public double getBeta() {
      return beta;
   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$\mu@f$,
    * @f$\beta@f$].
    */
   public double[] getParams () {
      double[] retour = {mu, beta};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : mu = " + mu+ ", beta = " + beta;
   }

}