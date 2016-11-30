/*
 * Class:        RootFinder
 * Description:  Provides methods to solve non-linear equations.
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
package umontreal.ssj.util;
   import umontreal.ssj.functions.MathFunction;

/**
 * This class provides methods to solve non-linear equations.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class RootFinder {
   private static final double MINVAL = 5.0e-308;
   private RootFinder() {}

   /**
    * Computes a root @f$x@f$ of the function in `f` using the
    * Brent-Dekker method. The interval @f$[a, b]@f$ must contain the root
    * @f$x@f$. The calculations are done with an approximate relative
    * precision `tol`. Returns @f$x@f$ such that @f$f(x) = 0@f$.
    *  @param a            left endpoint of initial interval
    *  @param b            right endpoint of initial interval
    *  @param f            the function which is evaluated
    *  @param tol          accuracy goal
    *  @return the root @f$x@f$
    */
   public static double brentDekker (double a, double b,
                                     MathFunction f, double tol) {
      final double EPS = 0.5E-15;
      final int MAXITER = 120;    // Maximum number of iterations
      double c, d, e;
      double fa, fb, fc;
      final boolean DEBUG = false;

      // Special case I = [b, a]
      if (b < a) {
         double ctemp = a;
         a = b;
         b = ctemp;
      }

      // Initialization
      fa = f.evaluate (a);
      if (Math.abs (fa) <= MINVAL)
         return a;
      
      fb = f.evaluate (b);
      if (Math.abs (fb) <= MINVAL)
         return b;
      
      c = a;
      fc = fa;
      d = e = b - a;
      tol += EPS + Num.DBL_EPSILON; // in case tol is too small

      if (Math.abs (fc) < Math.abs (fb)) {
         a = b;
         b = c;
         c = a;
         fa = fb;
         fb = fc;
         fc = fa;
      }

      int i;
      for (i = 0; i < MAXITER; i++) {
         double s, p, q, r;
         double tol1 = tol + 4.0 * Num.DBL_EPSILON * Math.abs (b);
         double xm = 0.5 * (c - b);
         if (DEBUG) {
            double err = Math.abs(fa - fb);
            System.out.printf("[a, b] = [%g, %g]   fa = %g,   fb = %g   |fa - fb| = %.2g%n",
                    a, b, fa, fb, err);
         }

         if (Math.abs (fb) <= MINVAL) {
            return b;
         }
         if (Math.abs (xm) <= tol1) {
            if (Math.abs (b) > MINVAL)
               return b;
            else
               return 0;
         }

         if ((Math.abs (e) >= tol1) && (Math.abs (fa) > Math.abs (fb))) {
            if (a != c) {
               // Inverse quadratic interpolation
               q = fa / fc;
               r = fb / fc;
               s = fb / fa;
               p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
               q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            } else {
               // Linear interpolation
               s = fb / fa;
               p = 2.0 * xm * s;
               q = 1.0 - s;
            }

            // Adjust signs
            if (p > 0.0)
               q = -q;
            p = Math.abs (p);

            // Is interpolation acceptable ?
            if (((2.0 * p) >= (3.0 * xm * q - Math.abs (tol1 * q)))
                  || (p >= Math.abs (0.5 * e * q))) {
               d = xm;
               e = d;
            } else {
               e = d;
               d = p / q;
            }
         } else {
            // Bisection necessary
            d = xm;
            e = d;
         }

         a = b;
         fa = fb;
         if (Math.abs (d) > tol1)
            b += d;
         else if (xm < 0.0)
            b -= tol1;
         else
            b += tol1;
         fb = f.evaluate (b);
         if (fb * (Math.signum (fc)) > 0.0) {
            c = a;
            fc = fa;
            d = e = b - a;
         } else {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
         }
      }

      if (i >= MAXITER)
         System.err.println(" WARNING:  root finding does not converge");
      return b;
   }

   /**
    * Computes a root @f$x@f$ of the function in `f` using the *bisection*
    * method. The interval @f$[a, b]@f$ must contain the root @f$x@f$. The
    * calculations are done with an approximate relative precision `tol`.
    * Returns @f$x@f$ such that @f$f(x) = 0@f$.
    *  @param a            left endpoint of initial interval
    *  @param b            right endpoint of initial interval
    *  @param f            the function which is evaluated
    *  @param tol          accuracy goal
    *  @return the root @f$x@f$
    */
   public static double bisection (double a, double b,
                                   MathFunction f, double tol) {
      // Case I = [b, a]
      if (b < a) {
         double ctemp = a;
         a = b;
         b = ctemp;
      }
      double xa = a;
      double xb = b;
      double yb = f.evaluate (b);
      // do preliminary checks on the bounds
      if (Math.abs (yb) <= MINVAL)
         return b;
      double ya = f.evaluate (a);
      if (Math.abs (ya) <= MINVAL)
         return a;
      
      double x = 0, y = 0;
      final int MAXITER = 1200;   // Maximum number of iterations
      final boolean DEBUG = false;
      tol += Num.DBL_EPSILON; // in case tol is too small

      if (DEBUG)
         System.out.println
         ("\niter              xa                   xb              f(x)");

      boolean fini = false;
      int i = 0;
      while (!fini) {
         x = (xa + xb) / 2.0;
         y = f.evaluate (x);
         if ((Math.abs (y) <= MINVAL) ||
             (Math.abs (xb - xa) <= tol * Math.abs (x)) ||
             (Math.abs (xb - xa) <= MINVAL)) {
            if (Math.abs(x) > MINVAL)
               return x;
            else
               return 0;
         }
         if (y * ya < 0.0)
            xb = x;
         else
            xa = x;
         ++i;
         if (DEBUG)
            System.out.printf("%3d    %18.12g     %18.12g    %14.4g%n",
                              i, xa, xb, y);
         if (i > MAXITER) {
            System.out.println ("***** bisection:  SEARCH DOES NOT CONVERGE");
            fini = true;
         }
      }
      return x;
   }

}