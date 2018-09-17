/*
 * Class:        PowerDist
 * Description:  power distribution
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
 * Extends the class  @ref ContinuousDistribution for the *power*
 * distribution @cite tEVA00a&thinsp; (page 161) with shape parameter @f$c >
 * 0@f$, over the interval @f$[a,b]@f$, where @f$a < b@f$.
 * This distribution has density
 * @anchor REF_probdist_PowerDist_eq_fpower
 * @f[
 *   f(x) = \frac{c(x-a)^{c - 1}}{(b - a)^c}, \qquad\mbox{for } a \le x \le b, \tag{fpower}
 * @f]
 * and @f$f(x) = 0@f$ elsewhere. Its distribution function is
 * @anchor REF_probdist_PowerDist_eq_Fpower
 * @f[
 *   F(x) = \frac{(x - a)^c}{(b - a)^c}, \qquad\mbox{for } a \le x \le b, \tag{Fpower}
 * @f]
 * with @f$F(x) = 0@f$ for @f$x \le a@f$ and @f$F(x) = 1@f$ for @f$x
 * \ge b@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class PowerDist extends ContinuousDistribution {
   private double a;
   private double b;
   private double c;

   /**
    * Constructs a `PowerDist` object with parameters @f$a =@f$ `a`, @f$b
    * =@f$ `b` and @f$c =@f$ `c`.
    */
   public PowerDist (double a, double b, double c) {
      setParams (a, b, c);
   }

   /**
    * Constructs a `PowerDist` object with parameters @f$a = 0@f$, @f$b
    * =@f$ `b` and @f$c =@f$ `c`.
    */
   public PowerDist (double b, double c) {
      setParams (0.0, b, c);
   }

   /**
    * Constructs a `PowerDist` object with parameters @f$a = 0@f$, @f$b
    * =1@f$ and @f$c =@f$ `c`.
    */
   public PowerDist (double c) {
      setParams (0.0, 1.0, c);
   }


   public double density (double x) {
      return density (a, b, c, x);
   }

   public double cdf (double x) {
      return cdf (a, b, c, x);
   }

   public double barF (double x) {
      return barF (a, b, c, x);
   }

   public double inverseF (double u) {
      return inverseF (a, b, c, u);
   }

   public double getMean() {
      return PowerDist.getMean (a, b, c);
   }

   public double getVariance() {
      return PowerDist.getVariance (a, b, c);
   }

   public double getStandardDeviation() {
      return PowerDist.getStandardDeviation (a, b, c);
   }

/**
 * Computes the density function (
 * {@link REF_probdist_PowerDist_eq_fpower fpower} ).
 *  @param a            left limit of interval
 *  @param b            right limit of interval
 *  @param c            shape parameter
 *  @param x            the value at which the density is evaluated
 *  @return returns the density function
 */
public static double density (double a, double b, double c, double x) {
      if (c <= 0.0)
         throw new IllegalArgumentException ("c <= 0");
      if (x <= a)
         return 0.0;
      if (x >= b)
         return 0.0;
      double z = (x-a)/(b-a);
      return c*Math.pow(z, c-1.0) / (b-a);
   }

   /**
    * Computes the distribution function (
    * {@link REF_probdist_PowerDist_eq_Fpower Fpower} ).
    *  @param a            left limit of interval
    *  @param b            right limit of interval
    *  @param c            shape parameter
    *  @param x            the value at which the distribution is
    *                      evaluated
    *  @return returns the distribution function
    */
   public static double cdf (double a, double b, double c, double x) {
      if (c <= 0.0)
         throw new IllegalArgumentException ("c <= 0");
      if (x <= a)
         return 0.0;
      if (x >= b)
         return 1.0;
      return Math.pow((x-a)/(b-a), c);
   }

   /**
    * Computes the complementary distribution function.
    *  @param a            left limit of interval
    *  @param b            right limit of interval
    *  @param c            shape parameter
    *  @param x            the value at which the complementary
    *                      distribution is evaluated
    *  @return returns the complementary distribution function
    */
   public static double barF (double a, double b, double c, double x) {
      if (c <= 0.0)
         throw new IllegalArgumentException ("c <= 0");
      if (x <= a)
         return 1.0;
      if (x >= b)
         return 0.0;
      return 1.0 - Math.pow((x-a)/(b-a),c);
   }

   /**
    * Computes the inverse of the distribution function.
    *  @param a            left limit of interval
    *  @param b            right limit of interval
    *  @param c            shape parameter
    *  @param u            the value at which the inverse distribution is
    *                      evaluated
    *  @return returns the inverse of the distribution function
    */
   public static double inverseF (double a, double b, double c, double u) {
      if (c <= 0.0)
         throw new IllegalArgumentException ("c <= 0");
      if (u < 0.0 || u > 1.0)
          throw new IllegalArgumentException ("u not in [0, 1]");
      if (u == 0.0)
         return a;
      if (u == 1.0)
         return b;

      return a + (b-a) * Math.pow(u,1.0/c);
   }

   /**
    * Estimates the parameter @f$c@f$ of the power distribution from the
    * @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$, using the
    * maximum likelihood method and assuming that @f$a@f$ and @f$b@f$ are
    * known. The estimate is returned in a one-element array: [@f$c@f$].
    * The maximum likelihood estimator is the value @f$\hat{c}@f$ that
    * satisfies the equation
    * @f{align*}{
    *    \frac{1}{\hat{c}} = -\frac{1}{n} \sum_{i=1}^n \ln\left(\frac{x_i - a}{b - a} \right)
    * @f}
    * @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @param a            left limit of interval
    *  @param b            right limit of interval
    *  @return returns the shape parameter [@f$\hat{c}@f$]
    */
   public static double[] getMLE (double[] x, int n, double a, double b) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double d = b - a;
      double somme = 0;
      for (int i = 0 ; i < n ; ++i) somme += Math.log((x[i] - a)/d);

      double [] parametres = new double [1];
      parametres[0] = -1.0 / (somme/n);
      return parametres;
   }

   /**
    * Creates a new instance of a power distribution with parameters
    * @f$a@f$ and @f$b@f$, with @f$c@f$ estimated using the maximum
    * likelihood method based on the @f$n@f$ observations @f$x[i]@f$, @f$i
    * = 0, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @param a            left limit of interval
    *  @param b            right limit of interval
    */
   public static PowerDist getInstanceFromMLE (double[] x, int n,
                                               double a, double b) {
      double parameters[] = getMLE (x, n, a, b);
      return new PowerDist (a, b, parameters[0]);
   }

   /**
    * Returns the mean @f$a + (b-a)c/(c+1)@f$ of the power distribution
    * with parameters @f$a@f$, @f$b@f$ and @f$c@f$.
    *  @param a            left limit of interval
    *  @param b            right limit of interval
    *  @param c            shape parameter
    *  @return returns the mean
    */
   public static double getMean (double a, double b, double c) {
      return a + (b-a) * c / (c+1.0);
   }

   /**
    * Computes and returns the variance @f$(b-a)^2 c / [(c+1)^2(c+2)]@f$
    * of the power distribution with parameters @f$a@f$, @f$b@f$ and
    * @f$c@f$.
    *  @param a            left limit of interval
    *  @param b            right limit of interval
    *  @param c            shape parameter
    *  @return returns the variance
    */
   public static double getVariance (double a, double b, double c) {
      return (b-a)*(b-a)*c / ((c+1.0)*(c+1.0)*(c+2.0));
   }

   /**
    * Computes and returns the standard deviation of the power
    * distribution with parameters @f$a@f$, @f$b@f$ and @f$c@f$.
    *  @return the standard deviation of the power distribution
    */
   public static double getStandardDeviation (double a, double b, double c) {
      return Math.sqrt (PowerDist.getVariance (a, b, c));
   }

   /**
    * Returns the parameter @f$a@f$.
    *  @return the left limit of interval @f$a@f$
    */
   public double getA() {
      return a;
   }

   /**
    * Returns the parameter @f$b@f$.
    *  @return the right limit of interval @f$b@f$
    */
   public double getB() {
      return b;
   }

   /**
    * Returns the parameter @f$c@f$.
    *  @return the shape parameter @f$c@f$
    */
   public double getC() {
      return c;
   }

   /**
    * Sets the parameters @f$a@f$, @f$b@f$ and @f$c@f$ for this object.
    *  @param a            left limit of interval
    *  @param b            right limit of interval
    *  @param c            shape parameter
    */
   public void setParams (double a, double b, double c) {
      this.a  = a;
      this.b  = b;
      this.c  = c;
   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$a@f$, @f$b@f$,
    * @f$c@f$].
    *  @return [@f$a@f$, @f$b, @f$c]
    */
   public double[] getParams () {
      double[] retour = {a, b, c};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    *  @return a `String` containing information about the current
    * distribution
    */
   public String toString () {
      return getClass().getSimpleName() + " : a = " + a + " : b = " + b + " : c = " + c;
   }

}