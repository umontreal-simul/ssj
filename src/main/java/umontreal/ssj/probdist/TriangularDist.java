/*
 * Class:        TriangularDist
 * Description:  triangular distribution
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

import java.util.Arrays;

/**
 * Extends the class  @ref ContinuousDistribution for the *triangular*
 * distribution (see @cite tJOH95b&thinsp; (page 297) and
 * @cite sLAW00a&thinsp; (page 317)) with domain @f$[a,b]@f$ and *mode* (or
 * shape parameter) @f$m@f$, where @f$a\le m\le b@f$. The density function is
 * @anchor REF_probdist_TriangularDist_eq_ftrian
 * @f[
 *   f(x) = \left\{\begin{array}{ll}
 *    \frac{2(x-a)}{(b-a)(m-a)} 
 *    & 
 *    \mbox{ if } a\le x\le m, 
 *    \\ 
 *    \frac{2(b-x)}{(b-a)(b-m)} 
 *    & 
 *    \mbox{ if } m\le x\le b, 
 *    \\ 
 *    0 
 *    & 
 *    \mbox{ elsewhere, } 
 *   \end{array}\right. \tag{ftrian}
 * @f]
 * the distribution function is
 * @f[
 *   F (x) = \left\{\begin{array}{ll}
 *    0 
 *    & 
 *    \mbox{ for } x < a, 
 *    \\ 
 *    \frac{(x - a)^2}{(b - a)(m - a)} 
 *    & 
 *    \mbox{ if } a\le x\le m, 
 *    \\ 
 *    1 - \frac{(b - x)^2}{(b - a)(b - m)} 
 *    & 
 *    \mbox{ if } m\le x\le b, 
 *    \\ 
 *    1 
 *    & 
 *    \mbox{ for } x > b, 
 *   \end{array}\right.
 * @f]
 * and the inverse distribution function is given by
 * @f[
 *   F^{-1}(u) = \left\{\begin{array}{ll}
 *    a + \sqrt{(b - a)(m - a)u} 
 *    & 
 *    \mbox{ if } 0\le u\le\frac{m-a}{b-a}, 
 *    \\ 
 *    b - \sqrt{(b - a)(b - m)(1 - u)} 
 *    & 
 *    \mbox{ if } \frac{m-a}{b-a}\le u \le1. 
 *   \end{array}\right.
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class TriangularDist extends ContinuousDistribution {
   private double a;
   private double b;
   private double m;

   /**
    * Constructs a `TriangularDist` object with default parameters
    * @f$a=0@f$, @f$b=1@f$, and @f$m=0.5@f$.
    */
   public TriangularDist() {
      setParams (0.0, 1.0, 0.5);
   }

   /**
    * Constructs a `TriangularDist` object with parameters @f$a = 0@f$ ,
    * @f$b = 1@f$ and @f$m@f$ = `m`.
    */
   public TriangularDist (double m) {
      setParams (0.0, 1.0, m);
   }

   /**
    * Constructs a `TriangularDist` object with parameters @f$a@f$,
    * @f$b@f$ and @f$m@f$.
    */
   public TriangularDist (double a, double b, double m) {
      setParams (a, b, m);
   }


   public double density (double x) {
      return density (a, b, m, x);
   }
       
   public double cdf (double x) {
      return cdf (a, b, m, x);
   }

   public double barF (double x) {
      return barF (a, b, m, x);
   }
     
   public double inverseF (double u){
      return inverseF (a, b, m, u);
   }

   public double getMean() {
      return TriangularDist.getMean (a, b, m);
   }

   public double getVariance() {
      return TriangularDist.getVariance (a, b, m);
   }

   public double getStandardDeviation() {
      return TriangularDist.getStandardDeviation (a, b, m);
   }

/**
 * Computes the density function.
 */
public static double density (double a, double b, double m, double x) {
      if (m < a || m > b)
         throw new IllegalArgumentException ("m is not in [a,b]");
      if (x < a || x > b)
         return 0.0;
      else if (x <= m && m != a)
         return 2.0*(x - a)/((b - a)*(m - a));
      else
         return 2.0*(b - x)/((b - a)*(b - m));
   }

   /**
    * Computes the distribution function.
    */
   public static double cdf (double a, double b, double m, double x) {
      if (m < a || m > b)
         throw new IllegalArgumentException ("m is not in [a,b]");
      if (x <= a)
         return 0.0;
      else if (x <= m && m != a)
         return (x - a)*(x - a)/((b - a)*(m - a));
      else if (x < b)
         return 1.0 - (b - x)*(b - x)/((b - a)*(b - m));
      else
         return 1.0;
   }

   /**
    * Computes the complementary distribution function.
    */
   public static double barF (double a, double b, double m, double x) {
      if (m < a || m > b)
         throw new IllegalArgumentException ("m is not in [a,b]");
      if (x <= a)
         return 1.0;
      else if (x <= m && m != a)
         return 1.0 - (x - a)*(x - a)/((b - a)*(m - a));
      else if (x < b)
         return (b - x)*(b - x)/((b - a)*(b - m));
      else
         return 0.0;
   }

   /**
    * Computes the inverse distribution function.
    */
   public static double inverseF (double a, double b, double m, double u) {
      if (m < a || m > b)
         throw new IllegalArgumentException ("m is not in [a,b]");
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u is not in [0,1]");
      if (u <= 0.0)
         return a;
      if (u >= 1.0)
         return b;
       // the code is taken and adapted from unuran
       // file /distributions/c_triangular_gen.c
       double h = (m - a)/(b - a);
       return u <= h && m != a ? a + Math.sqrt ((b - a)*(m - a)*u) 
                : b - Math.sqrt ((b - a)*(b - m)*(1 - u));
   }

   /**
    * Estimates the parameter @f$m@f$ of the triangular distribution using
    * the maximum likelihood method, from the @f$n@f$ observations
    * @f$x[i]@f$, @f$i = 0, 1,…, n-1@f$. The estimated parameter is
    * returned in a one-element array: [@f$\hat{m}@f$]. See @cite tOLI72a,
    * @cite tHUA07a, @cite tKOT04a&thinsp;.
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @param a            lower limit of range
    *  @param b            upper limit of range
    *  @return returns the parameter [@f$m@f$]
    */
   public static double[] getMLE (double[] x, int n, double a, double b) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      double[] Y = new double[n];   // sorted x[i]
      System.arraycopy (x, 0, Y, 0, n);
      Arrays.sort (Y);

      int rmax = -1;
      double prodmax = -1.0e300;
      final double ba = b - a;
      double z;
      int i;
      for (int r = 0; r < n; r++) {
         z = (Y[r] - a) / ba;
         if ((z <= (double)r/n) || (z >= (double)(r + 1)/n))
            continue;    // MLE cannot be there
         double prod = 1.0;
         double d = Y[r] - a;
         for (i = 0; i < r; i++)
            prod *= (Y[i] - a)/d;
         
         d = b - Y[r];
         for (i = r+1; i < n; i++)
            prod *= (b - Y[i])/d;

         if (prod > prodmax) {
            prodmax = prod;
            rmax = r;
         }
      }
      
      if (rmax < 0)
         throw new UnsupportedOperationException (
            "   data cannot fit a triangular distribution");

      double[] param = new double[1];
      param[0] = Y[rmax];
      return param;
   }

   /**
    * Creates a new instance of a triangular distribution with parameters
    * `a` and `b`. @f$m@f$ is estimated using the maximum likelihood
    * method based on the @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1, …,
    * n-1@f$.
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @param a            lower limit of range
    *  @param b            upper limit of range
    */
   public static TriangularDist getInstanceFromMLE (double[] x, int n,
                                                    double a, double b) {
      double param[] = getMLE (x, n, a, b);
      return new TriangularDist (a, b, param[0]);
   }

   /**
    * Computes and returns the mean @f$E[X] = (a + b + m)/3@f$ of the
    * triangular distribution with parameters @f$a@f$, @f$b@f$, @f$m@f$.
    *  @return the mean of the triangular distribution
    */
   public static double getMean (double a, double b, double m) {
      if ((a == 0.0 && b == 1.0) && (m < 0 || m > 1))
         throw new IllegalArgumentException ("m is not in [0,1]");
      else if (m < a || m > b) 
         throw new IllegalArgumentException ("m is not in [a,b]");

      return ((a + b + m) / 3.0);
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] = (a^2 + b^2 +
    * m^2 - ab - am - bm)/18@f$ of the triangular distribution with
    * parameters @f$a@f$, @f$b@f$, @f$m@f$.
    *  @return the variance of the triangular distribution
    */
   public static double getVariance (double a, double b, double m) {
      if ((a == 0.0 && b == 1.0) && (m < 0 || m > 1))
         throw new IllegalArgumentException ("m is not in [0,1]");
      else if (m < a || m > b) 
         throw new IllegalArgumentException ("m is not in [a,b]");

      return ((a * a + b * b + m * m - a * b - a * m - b * m) / 18.0);
   }

   /**
    * Computes and returns the standard deviation of the triangular
    * distribution with parameters @f$a@f$, @f$b@f$, @f$m@f$.
    *  @return the standard deviation of the triangular distribution
    */
   public static double getStandardDeviation (double a, double b, double m) {
      return Math.sqrt (TriangularDist.getVariance (a, b, m));
   }

   /**
    * Returns the value of @f$a@f$ for this object.
    */
   public double getA() {
      return a;
   }

   /**
    * Returns the value of @f$b@f$ for this object.
    */
   public double getB() {
      return b;
   }

   /**
    * Returns the value of @f$m@f$ for this object.
    */
   public double getM() {
      return m;
   }

   /**
    * Sets the value of the parameters @f$a@f$, @f$b@f$ and @f$m@f$ for
    * this object.
    */
   public void setParams (double a, double b, double m) {
      if ((a == 0.0 && b == 1.0) && (m < 0 || m > 1))
         throw new IllegalArgumentException ("m is not in [0,1]");
      else if (a >= b)
         throw new IllegalArgumentException ("a >= b");
      else if (m < a || m > b) 
         throw new IllegalArgumentException ("m is not in [a,b]");
      this.a = a;
      this.b = b;
      this.m = m;
      supportA = a;
      supportB = b;
   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$a@f$, @f$b@f$,
    * @f$m@f$].
    */
   public double[] getParams () {
      double[] retour = {a, b, m};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : a = " + a + ", b = " + b + ", m = " + m;
   }

}