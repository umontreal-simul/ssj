/*
 * Class:        RayleighDist
 * Description:  Rayleigh distribution
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
 * This class extends the class  @ref ContinuousDistribution for the
 * *Rayleigh* distribution @cite tEVA00a&thinsp; with location parameter
 * @f$a@f$, and scale parameter @f$\beta> 0@f$. The density function is
 * @anchor REF_probdist_RayleighDist_eq_frayleigh
 * @f[
 *   f(x) = \frac{(x-a)}{\beta^2}  e^{-(x-a)^2/(2\beta^2)} \qquad\mbox{for } x \ge a, \tag{frayleigh}
 * @f]
 * and @f$f(x) = 0@f$ for @f$x < a@f$. The distribution function is
 * @anchor REF_probdist_RayleighDist_eq_Frayleigh
 * @f[
 *   F(x) = 1 - e^{-(x - a)^2/(2\beta^2)} \qquad\mbox{for } x \ge a, \tag{Frayleigh}
 * @f]
 * and the inverse distribution function is
 * @anchor REF_probdist_RayleighDist_eq_Invrayleigh
 * @f[
 *   F^{-1}(u) = x = a + \beta\sqrt{-2\ln(1-u)} \qquad\mbox{for } 0 \le u \le1. \tag{Invrayleigh}
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class RayleighDist extends ContinuousDistribution {
   private double a;
   private double beta;

   /**
    * Constructs a `RayleighDist` object with parameters @f$a = 0@f$ and
    * @f$\beta@f$ = `beta`.
    */
   public RayleighDist (double beta) {
      setParams (0.0, beta);
   }

   /**
    * Constructs a `RayleighDist` object with parameters @f$a =@f$ `a`,
    * and @f$\beta@f$ = `beta`.
    */
   public RayleighDist (double a, double beta) {
      setParams (a, beta);
   }


   public double density (double x) {
      return density (a, beta, x);
   }

   public double cdf (double x) {
      return cdf (a, beta, x);
   }

   public double barF (double x) {
      return barF (a, beta, x);
   }

   public double inverseF (double u) {
      return inverseF (a, beta, u);
   }

   public double getMean() {
      return RayleighDist.getMean (a, beta);
   }

   public double getVariance() {
      return RayleighDist.getVariance (beta);
   }

   public double getStandardDeviation() {
      return RayleighDist.getStandardDeviation (beta);
   }

/**
 * Computes the density function (
 * {@link REF_probdist_RayleighDist_eq_frayleigh frayleigh} ).
 *  @param a            the location parameter
 *  @param beta         the scale parameter
 *  @param x            the value at which the density is evaluated
 *  @return the density function
 */
public static double density (double a, double beta, double x) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (x <= a)
         return 0.0;
      final double Z = (x - a)/beta;
      return Z/beta * Math.exp(-Z*Z/2.0);
   }

   /**
    * Same as `density (0, beta, x)`.
    *  @param beta         the scale parameter
    *  @param x            the value at which the density is evaluated
    *  @return returns the density function
    */
   public static double density (double beta, double x) {
      return density (0.0, beta, x);
   }

   /**
    * Computes the distribution function (
    * {@link REF_probdist_RayleighDist_eq_Frayleigh
    * Frayleigh} ).
    *  @param a            the location parameter
    *  @param beta         the scale parameter
    *  @param x            the value at which the distribution is
    *                      evaluated
    *  @return returns the distribution function
    */
   public static double cdf (double a, double beta, double x) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (x <= a)
         return 0.0;
      final double Z = (x - a)/beta;
      if (Z >= 10.0)
         return 1.0;
      return -Math.expm1(-Z*Z/2.0);
   }

   /**
    * Same as `cdf (0, beta, x)`.
    *  @param beta         the scale parameter
    *  @param x            the value at which the distribution is
    *                      evaluated
    *  @return returns the distribution function
    */
   public static double cdf (double beta, double x) {
      return cdf (0.0, beta, x);
   }

   /**
    * Computes the complementary distribution function.
    *  @param a            the location parameter
    *  @param beta         the scale parameter
    *  @param x            the value at which the complementary
    *                      distribution is evaluated
    *  @return returns the complementary distribution function
    */
   public static double barF (double a, double beta, double x) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (x <= a)
         return 1.0;
      double z = (x - a)/beta;
      if (z >= 44.0)
         return 0.0;
      return Math.exp(-z*z/2.0);
   }

   /**
    * Same as `barF (0, beta, x)`.
    *  @param beta         the scale parameter
    *  @param x            the value at which the complementary
    *                      distribution is evaluated
    *  @return returns the complementary distribution function
    */
   public static double barF (double beta, double x) {
      return barF (0.0, beta, x);
   }

   /**
    * Computes the inverse of the distribution function (
    * {@link REF_probdist_RayleighDist_eq_Invrayleigh
    * Invrayleigh} ).
    *  @param a            the location parameter
    *  @param beta         the scale parameter
    *  @param u            the value at which the inverse distribution is
    *                      evaluated
    *  @return returns the inverse of the distribution function
    */
   public static double inverseF (double a, double beta, double u) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (u < 0.0 || u > 1.0)
          throw new IllegalArgumentException ("u not in [0, 1]");
      if (u <= 0.0)
         return a;
      if (u >= 1.0)
         return Double.POSITIVE_INFINITY;

      return a + beta * Math.sqrt(-2.0 * Math.log1p(-u));
   }

   /**
    * Same as `inverseF (0, beta, u)`.
    *  @param beta         the scale parameter
    *  @param u            the value at which the inverse distribution is
    *                      evaluated
    *  @return returns the inverse of the distribution function
    */
   public static double inverseF (double beta, double u) {
      return inverseF (0.0, beta, u);
   }

   /**
    * Estimates the parameter @f$\beta@f$ of the Rayleigh distribution
    * using the maximum likelihood method, assuming that @f$a@f$ is known,
    * from the @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$.
    * The estimate is returned in a one-element array:
    * [@f$\hat{\beta}@f$].  The maximum likelihood estimator is the value
    * @f$\hat{\beta}@f$ that satisfies the equation
    * @f[
    *   \hat{\beta} = \sqrt{\frac{1}{2n}\sum_{i=1}^n x_i^2}
    * @f]
    * @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @param a            the location parameter
    *  @return returns the parameter [@f$\hat{\beta}@f$]
    */
   public static double[] getMLE (double[] x, int n, double a) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double somme = 0;
      for (int i = 0 ; i < n ; ++i) somme += (x[i]-a)*(x[i]-a);

      double [] parametres = new double [1];
      parametres[0] = Math.sqrt(somme/(2.0*n));
      return parametres;
   }

   /**
    * Creates a new instance of a Rayleigh distribution with parameters
    * @f$a@f$ and @f$\hat{\beta}@f$. This last is estimated using the
    * maximum likelihood method based on the @f$n@f$ observations
    * @f$x[i]@f$, @f$i = 0, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @param a            the location parameter
    */
   public static RayleighDist getInstanceFromMLE (double[] x, int n,
                                                  double a) {
      double parameters[] = getMLE (x, n, a);
      return new RayleighDist (parameters[0], parameters[1]);
   }

   /**
    * Returns the mean @f$a + \beta\sqrt{\pi/2}@f$ of the Rayleigh
    * distribution with parameters @f$a@f$ and @f$\beta@f$.
    *  @param a            the location parameter
    *  @param beta         the scale parameter
    *  @return the mean of the Rayleigh distribution
    */
   public static double getMean (double a, double beta) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      return (a + beta * Math.sqrt(Math.PI/2.0));
   }

   /**
    * Returns the variance of the Rayleigh distribution with parameter
    * @f$\beta@f$.
    *  @param beta         the scale parameter
    *  @return the variance of the Rayleigh distribution
    */
   public static double getVariance (double beta) {
      if (beta == 0.0)
        throw new IllegalArgumentException ("beta = 0");
      return (2.0 - 0.5*Math.PI) * beta * beta;
   }

   /**
    * Returns the standard deviation @f$\beta\sqrt{2 - \pi/2}@f$ of the
    * Rayleigh distribution with parameter @f$\beta@f$.
    *  @param beta         the scale parameter
    *  @return the standard deviation of the Rayleigh distribution
    */
   public static double getStandardDeviation (double beta) {
      return Math.sqrt (RayleighDist.getVariance (beta));
   }

   /**
    * Returns the parameter @f$a@f$.
    *  @return the location parameter @f$a@f$
    */
   public double getA() {
      return a;
   }

   /**
    * Returns the parameter @f$\beta@f$.
    *  @return the scale parameter @f$beta@f$
    */
   public double getSigma() {
      return beta;
   }

   /**
    * Sets the parameters @f$a@f$ and @f$\beta@f$ for this object.
    *  @param a            the location parameter
    *  @param beta         the scale parameter
    */
   public void setParams (double a, double beta) {
      if (beta <= 0.0)
        throw new IllegalArgumentException ("beta <= 0");
      this.a  = a;
      this.beta  = beta;
      supportA = a;
   }

   /**
    * Return an array containing the parameters of the current
    * distribution in the order: [@f$a@f$, @f$\beta@f$].
    *  @return [@f$a@f$, @f$\beta@f$]
    */
   public double[] getParams () {
      double[] retour = {a, beta};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    *  @return a `String` containing information about the current
    * distribution
    */
   public String toString () {
      return getClass().getSimpleName() + " : a = " + a + ", beta = " + beta;
   }

}