/*
 * Class:        StudentDistQuick
 * Description:  Student t-distribution
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
// import umontreal.ssj.functions.MathFunction;

/**
 * Extends the class  @ref StudentDist for the *Student*
 * @f$t@f$-distribution. Uses methods that are faster but less precise than
 * @ref StudentDist.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class StudentDistQuick extends StudentDist {
    private static final int STUDENT_N1 = 20;
    private static final double STUDENT_X1 = 8.01;
    private static final int STUDENT_KMAX = 200;
    private static final double STUDENT_EPS = 0.5E-16;

   /**
    * Constructs a `StudentDistQuick` object with `n` degrees of freedom.
    */
   public StudentDistQuick (int n) {
      super (n);
   }

  /*  These public methods are necessary so that the methods cdf,
   *  barF and inverseF used are those of the present
   *  class and not those of the mother class.
   */

   public double cdf (double x) {
      return cdf (n, x);
   }

   public double barF (double x) {
      return barF (n, x);
   }

   public double inverseF (double u) {
      return inverseF (n, u);
   }

/**
 * Returns the approximation of @cite tKEN80a&thinsp; (page 96) of the
 * Student @f$t@f$-distribution function with @f$n@f$ degrees of freedom. Is
 * very poor in the tails but good in the central part of the range.
 */
public static double cdf (int n, double x) {
      if (n <= 2)
         return StudentDist.cdf(n, x);
      if (x == Double.NEGATIVE_INFINITY)
           return 0.0;
      if (x == Double.POSITIVE_INFINITY)
           return 1.0;
      double b, y, z, z2, prec, v;
      int k;

      // first case: small n and small x
      if (n <= STUDENT_N1 && x <= STUDENT_X1) {
         b = 1.0 + x*x/n;
         y = x/Math.sqrt ((double)n);
         z = 1.0;
         for (k = n - 2; k >= 2; k -= 2)
            z = 1.0 + z*(k - 1)/(k*b);
         if (n % 2 == 0) {
               v = (1.0 + z*y/Math.sqrt (b))/2.0;
         } else {
            if (y > -1.0)
               v = (0.5 + (Math.atan (y) + z*y/b)/Math.PI);
            else
               v = (Math.atan(-1.0/y) + z * y / b) / Math.PI;
         }
         if (v > 1.0e-18)
            return v;
         else
            return 0.0;
      }

      // second case: large n and small x
      else if (x < STUDENT_X1) {
         double a = n - 0.5;
         b = 48.0*a*a;
         z2 = a*Math.log1p (x*x/n);
         z = Math.sqrt (z2);
         y = (((((64.0*z2 + 788.0)*z2 + 9801.0)*z2 + 89775.0)*z2 +
               543375.0)*z2 + 1788885.0)*z/(210.0*b*b*b);
         y -= (((4.0*z2 + 33.0)*z2 + 240.0)*z2 +  855.0)*z/(10.0*b*b);
         y += z + (z2 + 3.0)*z/b;
         if (x >= 0.0)
            return NormalDist.barF01 (-y);
         else
            return NormalDist.barF01 (y);
      }

      // third case: large x
      else {
         // Compute the Student probability density
         b = 1.0 + x*x/n;
         y = Num.gammaRatioHalf (n/2.0);
         y *= 1.0/(Math.sqrt (Math.PI*n)*Math.pow (b, (n + 1)/2.0));

         y *= 2.0*Math.sqrt (n*b);
         z = y/n;
         k = 2;
         z2 = prec = 10.0;
         while (k < STUDENT_KMAX && prec > STUDENT_EPS) {
            y *= (k - 1)/(k*b);
            z += y/(n + k);
            prec = Math.abs (z - z2);
            z2 = z;
            k += 2;
         }
         if (k >= STUDENT_KMAX)
           System.err.println ("student: k >= STUDENT_KMAX");
         if (x >= 0.0)
            return 1.0 - z/2.0;
         else
            return z/2.0;
      }
   }

   /**
    * Computes the complementary distribution function @f$\bar{F}(x)@f$.
    */
   public static double barF (int n, double x) {
        if (n <= 2)
           return StudentDist.barF(n, x);
      return cdf (n, -x);
   }

   /**
    * Returns an approximation of @f$F^{-1}(u)@f$, where @f$F@f$ is the
    * Student @f$t@f$-distribution function with @f$n@f$ degrees of
    * freedom. Gives at least 5 decimal digits of precision when @f$n
    * \ge3@f$ (see @cite tHIL70a&thinsp;). Uses exact formulae for
    * @f$n=1@f$ and @f$n=2@f$.
    */
   public static double inverseF (int n, double u) {
        if (n <= 2)
           return StudentDist.inverseF(n, u);
      final double PI = Math.PI;
      double a, b, c, d, e, p, t, x, y;

      e = (double) n;
      if (u > 0.5)
         p = 2.0 * (1.0 - u);
      else
         p = 2.0 * u;

      a = 1. / (e - 0.5);
      b = 48. / (a * a);
      c = ((20700. / b * a - 98.) * a - 16.) * a + 96.36;
      d = e * Math.sqrt (a * PI / 2.) * ((94.5 / (b + c) - 3.) / b + 1.);
      y = Math.pow ((d * p), (2.0 / e));
      if (y > (a + 0.05)) {
         if (p == 1.0)
            x = 0.0;
         else
            x = NormalDist.inverseF01 (p * 0.5);
         y = x * x;
         if (n < 5)
            c = c + 0.3 * (e - 4.5) * (x + 0.6);

         c = (((0.05 * d * x - 5.) * x - 7.) * x - 2.) * x + b + c;
         y = (((((0.4 * y + 6.3) * y + 36.) * y + 94.5) /
               c - y - 3.) / b + 1.) * x;
         y = a * (y * y);
         y = Math.expm1 (y);

      } else {
         y = ((1. / (((e + 6.) / (e * y) - 0.089 * d - 0.822) *
                     (e + 2.) * 3.) + 0.5 / (e + 4.)) * y - 1.) *
               (e + 1.) / (e + 2.) + 1. / y;
      }

      t = Math.sqrt (e * y);
      if (u < 0.5)
         return -t;
      else
         return t;
    }

}