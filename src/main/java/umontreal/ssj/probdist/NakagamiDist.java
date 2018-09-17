/*
 * Class:        NakagamiDist
 * Description:  Nakagami distribution
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
package  umontreal.ssj.probdist;
import umontreal.ssj.util.*;
import umontreal.ssj.functions.MathFunction;

/**
 * Extends the class  @ref ContinuousDistribution for the *Nakagami*
 * distribution with location parameter @f$a@f$, scale parameter @f$\lambda>
 * 0@f$ and shape parameter @f$c > 0@f$. The density is
 * @anchor REF_probdist_NakagamiDist_eq_fnakagami
 * @f[
 *   f(x) = \frac{2\lambda^c}{\Gamma(c)} (x-a)^{2c-1} e^{-{\lambda}(x-a)^2} \qquad\mbox{for } x > a,\tag{fnakagami}
 * @f]
 * @f[
 *   f(x) = 0 \qquad\mbox{ for } x \le a,
 * @f]
 * where @f$\Gamma@f$ is the gamma function.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class NakagamiDist extends ContinuousDistribution {
   protected double a;              // Location parameter
   protected double lambda;         // Scale parameter
   protected double c;              // Shape parameter
   private double factor;
   private double ratio;            // Gamma(c + 1/2)/Gamma(c)

   /**
    * Constructs a `NakagamiDist` object with parameters @f$a =@f$ `a`,
    * @f$\lambda=@f$ `lambda` and @f$c =@f$ `c`.
    */
   public NakagamiDist (double a, double lambda, double c) {
      setParams (a, lambda, c);
   }


   public double density (double x) {
      if (x <= a) return 0.0;
      return 2.0 * Math.exp( factor
                             + Math.log(x-a)*(2.0*c-1.0)
                             - lambda*(x-a)*(x-a) );
   }

   public double cdf (double x) {
      return cdf (a, lambda, c, x);
   }

   public double barF (double x) {
      return barF (a, lambda, c, x);
   }

   public double inverseF (double u) {
      return inverseF (a, lambda, c, u);
   }

   public double getMean() {
      return a + ratio/Math.sqrt(lambda);
   }

   public double getVariance() {
      return (c - ratio*ratio)/lambda;
   }

   public double getStandardDeviation() {
      return Math.sqrt(getVariance ());
   }

/**
 * Computes the density function of the *Nakagami* distribution.
 *  @param a            the location parameter
 *  @param lambda       the scale parameter
 *  @param c            the shape parameter
 *  @param x            the value at which the density is evaluated
 *  @return returns the density function
 */
public static double density (double a, double lambda, double c,
                                 double x) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (c <= 0.0)
         throw new IllegalArgumentException ("c <= 0");
      if (x <= a)
         return 0.0;

      return 2.0 * Math.exp( Math.log(lambda)*c - Num.lnGamma(c)
                             +  Math.log(x-a)*(2.0*c-1.0)
                             - lambda*(x-a)*(x-a) );
   }

   /**
    * Computes the distribution function.
    *  @param a            the location parameter
    *  @param lambda       the scale parameter
    *  @param c            the shape parameter
    *  @param x            the value at which the distribution is
    *                      evaluated
    *  @return returns the cdf function
    */
   public static double cdf (double a, double lambda, double c, double x) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (c <= 0.0)
         throw new IllegalArgumentException ("c <= 0");
      if (x <= a)
         return 0.0;

      return GammaDist.cdf(c, 12, lambda*(x-a)*(x-a));
   }

   /**
    * Computes the complementary distribution function.
    *  @param a            the location parameter
    *  @param lambda       the scale parameter
    *  @param c            the shape parameter
    *  @param x            the value at which the complementary
    *                      distribution is evaluated
    *  @return returns the complementary distribution function
    */
   public static double barF (double a, double lambda, double c, double x) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (c <= 0.0)
         throw new IllegalArgumentException ("c <= 0");
      if (x <= a)
         return 1.0;
      return GammaDist.barF(c, 12, lambda*(x-a)*(x-a));
   }

   /**
    * Computes the inverse of the distribution function.
    *  @param a            the location parameter
    *  @param lambda       the scale parameter
    *  @param c            the shape parameter
    *  @param u            the value at which the inverse distribution is
    *                      evaluated
    *  @return returns the inverse distribution function
    */
   public static double inverseF (double a, double lambda, double c,
                                  double u) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (c <= 0.0)
         throw new IllegalArgumentException ("c <= 0");
      if (u > 1.0 || u < 0.0)
         throw new IllegalArgumentException ("u not in [0,1]");
      if (u <= 0.0) return a;
      if (u >= 1.0)
         return Double.POSITIVE_INFINITY;
      double res = GammaDist.inverseF(c, 12, u);
      return a + Math.sqrt(res/lambda);
   }

   /**
    * Computes and returns the mean
    * @f[
    *   E[X] = a + \frac{1}{\sqrt{\lambda}}\; \frac{\Gamma(c+1/2)}{\Gamma(c)}.
    * @f]
    * @param a            the location parameter
    *  @param lambda       the scale parameter
    *  @param c            the shape parameter
    *  @return returns the mean
    */
   public static double getMean (double a, double lambda, double c) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (c <= 0.0)
         throw new IllegalArgumentException ("c <= 0");
      return a + Num.gammaRatioHalf(c) / Math.sqrt(lambda);
   }

   /**
    * Computes and returns the variance
    * @f[
    *   \mbox{Var}[X] = \frac{1}{\lambda} \left[c - \left(\frac{\Gamma(c+1/2)}{\Gamma(c)}\right)^2\right].
    * @f]
    * @param a            the location parameter
    *  @param lambda       the scale parameter
    *  @param c            the shape parameter
    *  @return returns the variance
    */
   public static double getVariance (double a, double lambda, double c) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (c <= 0.0)
         throw new IllegalArgumentException ("c <= 0");
      double rat = Num.gammaRatioHalf(c);
      return (c - rat*rat) / lambda;
   }

   /**
    * Computes the standard deviation of the Nakagami distribution with
    * parameters @f$a@f$, @f$\lambda@f$ and @f$c@f$.
    *  @param a            the location parameter
    *  @param lambda       the scale parameter
    *  @param c            the shape parameter
    *  @return returns the standard deviation
    */
   public static double getStandardDeviation (double a, double lambda,
                                              double c) {
      return Math.sqrt (NakagamiDist.getVariance (a, lambda, c));
   }

   /**
    * Returns the location parameter @f$a@f$ of this object.
    *  @return returns the location parameter
    */
   public double getA() {
      return a;
   }

   /**
    * Returns the scale parameter @f$\lambda@f$ of this object.
    *  @return returns the scale parameter
    */
   public double getLambda() {
      return lambda;
   }

   /**
    * Returns the shape parameter @f$c@f$ of this object.
    *  @return returns the shape parameter
    */
   public double getC() {
      return c;
   }

    /**
     * Sets the parameters @f$a@f$, @f$\lambda@f$ and @f$c@f$ of this
     * object.
     *  @param a            the location parameter
     *  @param lambda       the scale parameter
     *  @param c            the shape parameter
     */
    public void setParams (double a, double lambda, double c) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (c <= 0.0)
         throw new IllegalArgumentException ("c <= 0");
      this.a = a;
      this.lambda = lambda;
      this.c = c;
      factor = (Math.log(lambda)*c - Num.lnGamma(c));
      ratio = Num.gammaRatioHalf(c);
    }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$a@f$,
    * @f$\lambda@f$, @f$c@f$].
    *  @return returns the parameters [@f$a@f$, @f$\lambda@f$, @f$c@f$]
    */
   public double[] getParams () {
      double[] retour = {a, lambda, c};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    *  @return returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : a = " + a + ", lambda = " + lambda
                                  + ", c = " + c;
   }

}