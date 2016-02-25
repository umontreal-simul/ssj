/*
 * Class:        StudentDistDist
 * Description:  Student-t distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        March 2009
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
import umontreal.ssj.functions.MathFunction;

/**
 * Extends the class  @ref ContinuousDistribution for the *Student*
 * @f$t@f$-distribution @cite tJOH95b&thinsp; (page 362) with @f$n@f$ degrees
 * of freedom, where @f$n@f$ is a positive integer. Its density is
 * @anchor REF_probdist_StudentDist_eq_fstudent
 * @f[
 *   f (x) = \frac{\Gamma\left((n + 1)/2 \right)}{\Gamma(n/2) \sqrt{\pi n}} \left(1 + \frac{x^2}{n}\right)^{-(n+1)/2} \qquad\mbox{for } -\infty< x < \infty, \tag{fstudent}
 * @f]
 * where @f$\Gamma(x)@f$ is the gamma function defined in (
 * {@link REF_probdist_GammaDist_eq_Gamma Gamma} ).
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class StudentDist extends ContinuousDistribution {
   protected int n;
   private double factor;
   private static final int NLIM1 = 100000;
/*
   private static double cdfPeizer (int n, double x) {
      // Peizer-Pratt normal approximation for the cdf (n, u)
      // @cite tPEI68a
      double v = Math.log1p(x*x/n) / (n - 5.0/6.0);
      double z = -(n - 2.0/3.0 + 0.1/n) * Math.sqrt(v);
      double u = NormalDist.cdf01 (z);
      if (x >= 0.0)
         return 1.0 - u;
      return u;
   }

   private static double invPeizer (int n, double u) {
      // Peizer-Pratt normal approximation for the inverseF (n, u)
      // @cite tPEI68a
      double z = NormalDist.inverseF01 (u);
      double q = z / (n - 2.0/3.0 + 0.1/n);
      double v = q*q*(n - 5.0/6.0);
      double t = Math.sqrt(n * Math.expm1(v));
      if (u >= 0.5)
         return t;
      else
         return -t;
   }
*/

   private static double cdfGaver (int n, double x) {
      // Gaver-Kafadar normal approximation for the cdf
      // @cite tGAV84a
      double v = Math.log1p(x * x / n) / (n - 1.5);
      double z = -(n - 1) * Math.sqrt(v);
      double u = NormalDist.cdf01 (z);
      if (x >= 0.0)
         return 1.0 - u;
      return u;
   }


   private static double invGaver (int n, double u) {
      // Gaver-Kafadar normal approximation for the inverse
      // @cite tGAV84a
      double z = NormalDist.inverseF01 (u);
      double q = z / (n - 1.0);
      double v = q * q * (n - 1.5);
      double t = Math.sqrt(n * Math.expm1(v));
      if (u >= 0.5)
         return t;
      else
         return -t;
   }


   private static class Function implements MathFunction {
      private int n;
      private double[] xi;

      public Function (double[] x, int n) {
         this.n = n;
         this.xi = new double[n];
         System.arraycopy (x, 0, this.xi, 0, n);
      }

      public double evaluate (double x) {
         if (x <= 0.0)
            return 1e200;
         double sum = 0.0;
         for (int i = 0; i < n; i++)
            sum += Math.log (density ((int) Math.round (x), xi[i]));
         return sum;
      }
   }

   /**
    * Constructs a `StudentDist` object with `n` degrees of freedom.
    */
   public StudentDist (int n) {
     setN (n);
   }


   public double density (double x) {
      return factor*Math.pow (1.0 / (1.0 + x*x/n), (n + 1)/2.0);
   }

   public double cdf (double x) {
      return cdf (n, x);
   }

   public double barF (double x) {
      return barF (n, x);
   }

   public double inverseF (double u) {
      return inverseF (n, u);
   }

   public double getMean() {
      return StudentDist.getMean (n);
   }

   public double getVariance() {
      return StudentDist.getVariance (n);
   }

   public double getStandardDeviation() {
      return StudentDist.getStandardDeviation (n);
   }

/**
 * Computes the density function (
 * {@link REF_probdist_StudentDist_eq_fstudent fstudent} ) of a
 * Student @f$t@f$-distribution with @f$n@f$ degrees of freedom.
 */
public static double density (int n, double x) {
      double factor = Num.gammaRatioHalf(n/2.0)/ Math.sqrt (n*Math.PI);
      return factor*Math.pow (1.0 / (1.0 + x*x/n), (n + 1)/2.0);
   }

   /**
    * Computes the Student @f$t@f$-distribution function @f$u=F(x)@f$ with
    * @f$n@f$ degrees of freedom. Gives 13 decimal digits of precision for
    * @f$n \le10^5@f$. For @f$n > 10^5@f$, gives at least 6 decimal
    * digits of precision everywhere, and at least 9 decimal digits of
    * precision for all @f$u > 10^{-15}@f$.
    */
   public static double cdf (int n, double x) {
      if (n < 1)
        throw new IllegalArgumentException ("n < 1");
      if (n == 1)
         return CauchyDist.cdf(0, 1, x);

      if (x > 1.0e10)
         return 1.0;
      if (n > NLIM1)
         return cdfGaver(n, x);

      double r = Math.abs(x);
      if (r < 1.0e20)
         r = Math.sqrt (n + x*x);

      double z;
      if (x >= 0.0)
         z = 0.5*(1.0 + x/r);
      else
         z = 0.5*n/(r*(r - x));

      if (n == 2)
         return z;
      return BetaSymmetricalDist.cdf (0.5*n, 15, z);
   }

   /**
    * Same as  {@link #cdf(int,double) cdf(n, x)}.
    */
   @Deprecated
   public static double cdf2 (int n, int d, double x) {
      if (d <= 0)
         throw new IllegalArgumentException ("student2:   d <= 0");
      return cdf (n, x);
   }

   /**
    * Computes the complementary distribution function @f$v =
    * \bar{F}(x)@f$ with @f$n@f$ degrees of freedom. Gives 13 decimal
    * digits of precision for @f$n \le10^5@f$. For @f$n > 10^5@f$, gives
    * at least 6 decimal digits of precision everywhere, and at least 9
    * decimal digits of precision for all @f$v > 10^{-15}@f$.
    */
   public static double barF (int n, double x) {
      if (n < 1)
        throw new IllegalArgumentException ("n < 1");
      if (n == 1)
         return CauchyDist.barF(0, 1, x);

      if (n == 2) {
         double z = Math.abs(x);
         if (z < 1.0e20)
            z = Math.sqrt(2.0 + x*x);
         if (x <= 0.) {
            if (x < -1.0e10)
               return 1.0;
            return 0.5* (1.0 - x / z);
         } else
            return 1.0 / (z * (z + x));
      }

      return cdf (n, -x);
   }

   /**
    * Returns the inverse @f$x = F^{-1}(u)@f$ of Student
    * @f$t@f$-distribution function with @f$n@f$ degrees of freedom. Gives
    * 13 decimal digits of precision for @f$n \le10^5@f$, and at least 9
    * decimal digits of precision for @f$n > 10^5@f$.
    */
   public static double inverseF (int n, double u) {
        if (n < 1)
            throw new IllegalArgumentException ("Student:   n < 1");
        if (u > 1.0 || u < 0.0)
            throw new IllegalArgumentException ("Student:   u not in [0, 1]");
        if (u <= 0.0)
           return Double.NEGATIVE_INFINITY;
        if (u >= 1.0)
           return Double.POSITIVE_INFINITY;

        if (1 == n)
           return CauchyDist.inverseF(0, 1, u);

        if (2 == n)
           return (2.0*u - 1.0) / Math.sqrt(2.0*u*(1.0 - u));

        if (n > NLIM1)
           return invGaver(n, u);
        double z = BetaSymmetricalDist.inverseF (0.5*n, u);
        return (z - 0.5) * Math.sqrt(n / (z*(1.0 - z)));
   }

   /**
    * Estimates the parameter @f$n@f$ of the Student @f$t@f$-distribution
    * using the maximum likelihood method, from the @f$m@f$ observations
    * @f$x[i]@f$, @f$i = 0, 1,…, m-1@f$. The estimate is returned in a
    * one-element array.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param m            the number of observations to use to evaluate
    *                      parameters
    *  @return returns the parameter [@f$\hat{n}@f$]
    */
   public static double[] getMLE (double[] x, int m) {
      double sum = 0.0;
      double[] parameters = new double[1];

      if (m <= 0)
         throw new IllegalArgumentException ("m <= 0");

      double var = 0.0;
      for (int i = 0; i < m; i++)
         var += x[i] * x[i];
      var /= (double) m;

      Function f = new Function (x, m);

      double n0 = Math.round ((2.0 * var) / (var - 1.0));
      double fn0 = f.evaluate (n0);
      double min = fn0;
      double fn1 = f.evaluate (n0 + 1.0);
      double fn_1 = f.evaluate (n0 - 1.0);

      parameters[0] = n0;

      if (fn_1 > fn0) {
         double n = n0 - 1.0;
         double y;
         while (((y = f.evaluate (n)) > min) && (n >= 1.0)) {
            min = y;
            parameters[0] = n;
            n -= 1.0;
         }

      } else if (fn1 > fn0) {
         double n = n0 + 1.0;
         double y;
         while ((y = f.evaluate (n)) > min) {
            min = y;
            parameters[0] = n;
            n += 1.0;
         }
      }
      return parameters;
   }

   /**
    * Creates a new instance of a Student @f$t@f$-distribution with
    * parameter @f$n@f$ estimated using the maximum likelihood method
    * based on the @f$m@f$ observations @f$x[i]@f$, @f$i = 0, 1, …,
    * m-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param m            the number of observations to use to evaluate
    *                      parameters
    */
   public static StudentDist getInstanceFromMLE (double[] x, int m) {
      double parameters[] = getMLE (x, m);
      return new StudentDist ((int) parameters[0]);
   }

   /**
    * Returns the mean @f$E[X] = 0@f$ of the Student @f$t@f$-distribution
    * with parameter @f$n@f$.
    *  @return the mean of the Student @f$t@f$-distribution @f$E[X] = 0@f$
    */
   public static double getMean (int n) {
     if (n < 2)
        throw new IllegalArgumentException ("n <= 1");
      return 0;
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] = n/(n - 2)@f$ of
    * the Student @f$t@f$-distribution with parameter @f$n@f$.
    *  @return the variance of the Student @f$t@f$-distribution
    * @f$\mbox{Var}[X] = n / (n - 2)@f$
    */
   public static double getVariance (int n) {
      if (n < 3)
         throw new IllegalArgumentException("n <= 2");
      return (n / (n - 2.0));
   }

   /**
    * Computes and returns the standard deviation of the Student
    * @f$t@f$-distribution with parameter @f$n@f$.
    *  @return the standard deviation of the Student @f$t@f$-distribution
    */
   public static double getStandardDeviation (int n) {
      return Math.sqrt (StudentDist.getVariance (n));
   }

   /**
    * Returns the parameter @f$n@f$ associated with this object.
    */
   public int getN() {
      return n;
   }

   /**
    * Sets the parameter @f$n@f$ associated with this object.
    */
   public void setN (int n) {
     if (n <= 0)
        throw new IllegalArgumentException ("n <= 0");
      this.n = n;
      factor = Num.gammaRatioHalf(n/2.0) / Math.sqrt (n*Math.PI);
   }

   /**
    * Return a table containing the parameter of the current distribution.
    */
   public double[] getParams () {
      double[] retour = {n};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : n = " + n;
   }

}