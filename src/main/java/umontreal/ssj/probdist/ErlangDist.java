/*
 * Class:        ErlangDist
 * Description:  Erlang distribution
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

/**
 * Extends the class  @ref GammaDist for the special case of the *Erlang*
 * distribution with shape parameter @f$k > 0@f$ and scale parameter
 * @f$\lambda> 0@f$. This distribution is a special case of the gamma
 * distribution for which the shape parameter @f$k=\alpha@f$ is an integer.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class ErlangDist extends GammaDist {

   /**
    * Constructs a `ErlangDist` object with parameters @f$k@f$ = `k` and
    * @f$\lambda=1@f$.
    */
   public ErlangDist (int k) {
      super (k);
   }

   /**
    * Constructs a `ErlangDist` object with parameters @f$k@f$ = `k` and
    * @f$\lambda@f$ = `lambda`.
    */
   public ErlangDist (int k, double lambda) {
      super (k, lambda);
   }

   /**
    * Computes the density function.
    */
   public static double density (int k, double lambda, double x) {
      return density ((double)k, lambda, x);
   }

   /**
    * Computes the distribution function using the gamma distribution
    * function.
    */
   public static double cdf (int k, double lambda, int d, double x) {
      return cdf ((double)k, d, lambda*x);
   }

   /**
    * Computes the complementary distribution function.
    */
   public static double barF (int k, double lambda, int d, double x) {
      return barF ((double)k, d, lambda*x);
   }

   /**
    * Returns the inverse distribution function.
    */
   public static double inverseF (int k, double lambda, int d, double u) {
      return inverseF ((double)k, lambda, d, u);
   }

   /**
    * Estimates the parameters @f$(k,\lambda)@f$ of the Erlang
    * distribution using the maximum likelihood method, from the @f$n@f$
    * observations @f$x[i]@f$, @f$i = 0, 1,…, n-1@f$. The estimates are
    * returned in a two-element array, in regular order: [@f$k@f$,
    * @f$\lambda@f$].  The equations of the maximum likelihood are the
    * same as for the gamma distribution. The @f$k@f$ parameter is the
    * rounded value of the @f$\alpha@f$ parameter of the gamma
    * distribution, and the @f$\lambda@f$ parameter is equal to the
    * @f$\beta@f$ parameter of the gamma distribution.
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{k}@f$,
    * @f$\hat{\lambda}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      double parameters[] = GammaDist.getMLE (x, n);
      parameters[0] = Math.round (parameters[0]);
      return parameters;
   }

   /**
    * Creates a new instance of an Erlang distribution with parameters
    * @f$k@f$ and @f$\lambda@f$ estimated using the maximum likelihood
    * method based on the @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1, …,
    * n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static ErlangDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new ErlangDist ((int) parameters[0], parameters[1]);
   }

   /**
    * Computes and returns the mean, @f$E[X] = k/\lambda@f$, of the
    * Erlang distribution with parameters @f$k@f$ and @f$\lambda@f$.
    *  @return the mean of the Erlang distribution @f$E[X] = k /
    * \lambda@f$
    */
   public static double getMean (int k, double lambda) {
      if (k <= 0)
         throw new IllegalArgumentException ("k <= 0");
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      return (k / lambda);
   }

   /**
    * Computes and returns the variance, @f$\mbox{Var}[X] =
    * k/\lambda^2@f$, of the Erlang distribution with parameters @f$k@f$
    * and @f$\lambda@f$.
    *  @return the variance of the Erlang distribution @f$\mbox{Var}[X] =
    * k / \lambda^2@f$
    */
   public static double getVariance (int k, double lambda) {
      if (k <= 0)
         throw new IllegalArgumentException ("k <= 0");
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return (k / (lambda * lambda));
   }

   /**
    * Computes and returns the standard deviation of the Erlang
    * distribution with parameters @f$k@f$ and @f$\lambda@f$.
    *  @return the standard deviation of the Erlang distribution
    */
   public static double getStandardDeviation (int k, double lambda) {
      if (k <= 0)
         throw new IllegalArgumentException ("k <= 0");
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");

      return (Math.sqrt (k) / lambda);
   }

   /**
    * Returns the parameter @f$k@f$ for this object.
    */
   public int getK() {
      return (int) getAlpha();
   }

   /**
    * Sets the parameters @f$k@f$ and @f$\lambda@f$ of the distribution
    * for this object. Non-static methods are computed with a rough target
    * of `d` decimal digits of precision.
    */
   public void setParams (int k, double lambda, int d) {
      super.setParams (k, lambda, d);
   }

   /**
    * Return a table containing parameters of the current distribution.
    * This table is put in regular order: [@f$k@f$, @f$\lambda@f$].
    */
   public double[] getParams () {
      return super.getParams();
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : k = " + (int)super.getAlpha() + ", lambda = " + super.getLambda();
   }

}