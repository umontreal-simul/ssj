/*
 * Class:        NI1
 * Description:
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Université de Montréal
 * Organization: DIRO, Université de Montréal
 * @author       Nabil Channouf
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
package umontreal.ssj.probdistmulti.norta;
import umontreal.ssj.util.*;
import umontreal.ssj.probdist.*;

/**
 * Extends the class  @ref NortaInitDisc and implements the algorithm NI1. It
 * uses an algorithm based on Brent method for root-finding, which combines
 * root-bracketing, bisection and inverse quadratic interpolation. It calls
 * the method  #integ to compute the function @f$g_r@f$ given in (
 * {@link REF_probdistmulti_norta_NortaInitDisc_gr_M gr_M} ).
 * The search should be done in the interval @f$[-1,0]@f$ if
 * @f$r_X\in[-1,0]@f$, or @f$[0,1]@f$ if @f$r_X\in[0,1]@f$. At each
 * iteration, the algorithm halves the interval length and uses an accuracy
 * @f$\epsilon@f$ to find the root @f$\rho_Z@f$ of equation (
 * {@link REF_probdistmulti_norta_overview_fr fr} ).
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class NI1 extends NortaInitDisc
{
   private double tolerance; /* Desired accuracy for the root-finder
   				 algorithm (epsilon in paragraph
   				 "Method NI1" of section 3 in the paper).*/

   /**
    * Constructor with the target rank correlation `rX`, the two
    * discrete marginals `dist1` and `dist2`, the parameter for
    * truncation `tr` (see the constructor of class  @ref NortaInitDisc
    * ) and the specific parameter @f$\epsilon=@f$ `tolerance` defined
    * above for the algorithm NI1.
    */
   public NI1 (double rX,
               DiscreteDistributionInt dist1,
               DiscreteDistributionInt dist2,
               double tr,
               double tolerance) 
   {
      super(rX, dist1, dist2, tr);
      this.tolerance = tolerance;
      computeParams();
   }

   /**
    * Computes and returns the correlation @f$\rho_Z@f$ using the
    * algorithm NI1.
    */
   public double computeCorr()
   {
      final double ITMAX = 100;   // Maximum number of iterations
      final double EPS = 1.0e-15; // Machine accuracy
      double b; /* Latest iterate and closest approximation of the root
      		     (the returned solution at the end). */
      double a; // Previous iterate (previous value of b).
      double c; /* Previous or older iterate so that f(b) and f(c)
      	             have opposite signs (may coincide with a). */
      double fa, fb, fc; /* Evaluations of the function f at points
      			      a, b and c. */
      double tolerance1; /* Final tolerance wich involves the machine
      			      accuracy and the initial desired tolerance
      			      (epsilon). */
      /** a, b, c, fa, fb, fc, and tolerance correspond to a, b, c, f(a)
          f(b), f(c) and epsilon, respectively, in the paper
          (paragraph "Method NI1" of section 3). */
      double x1, x2; // Left and right endpoints of initial search interval
      double pp, q, rrr, s; /* Parameters to compute inverse quadratic
      				 interpolation. */
      double xm; // Bisection point.
      double min1, min2; /* Criterias to check whether inverse quadratic
      			      interpolation can be performed or not. */
      double e = 0.0, d = 0.0; /* Parameters to specify bounding interval
      				    and whether bisection or inverse
      				    quadratic interpolation was performed
      				    at one iteration before the last one. */
      // Precompute constants.
      double cc = rX * sd1 * sd2;
      double ccc = cc + mu1 * mu2;

      if (rX == 0)
         return 0.0;
      if (rX > 0) {  // Orient the search and initialize a, b, c
         x1 = 0.0;
         x2 = 1.0;
         a = x1;
         b = x2;
         c = x2;
         fa = - cc;
         fb = integ (b) - ccc;
      } else {
         x1 = -1.0;
         x2 = 0.0;
         a = x1;
         b = x2;
         c = x2;
         fa = integ (a) - ccc;
         fb = - cc;
      }
      fc = fb;

      for (int i = 1; i <= ITMAX; i++) {   // Begin the search
         if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            // Rename a, b, c and adjust bounding interval d
            c = a;
            fc = fa;
            e = d = b - a;
         }
         if (Math.abs (fc) < Math.abs (fb)) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
         }

         // Convergence check
         tolerance1 = 2.0 * EPS * Math.abs (b) + 0.5 * tolerance;
         xm = 0.5 * (c - b);
         if (Math.abs (xm) <= tolerance1 || fb == 0.0)
            return b;

         if (Math.abs (e) >= tolerance1 && Math.abs (fa) > Math.abs (fb)) {
            s = fb / fa;   // Attempt inverse quadratic interpolation
            if (a == c) {
               pp = 2.0 * xm * s;
               q = 1.0 - s;
            } else {
               q = fa / fc;
               rrr = fb / fc;
               pp = s * (2.0 * xm * q * (q - rrr) - (b - a) * (rrr - 1.0));
               q = (q - 1.0) * (rrr - 1.0) * (s - 1.0);
            }
            if (pp > 0.0)
               q = -q;            // Check whether in bounds
            pp = Math.abs (pp);
            min1 = 3.0 * xm * q - Math.abs (tolerance1 * q);
            min2 = Math.abs (e * q);
            if (2.0 * pp < (min1 < min2 ? min1 : min2)) {
               e = d;             // Accept interpolation
               d = pp / q;
            } else {        // Interpolation failed, use bisection
               d = xm;
               e = d;
            }
         } else {     // Bounds decreasing too slowly, use bisection
            d = xm;
            e = d;
         }
         a = b;                   // a becomes the best trial
         fa = fb;
         if (Math.abs (d) > tolerance1)
            b += d;              // Evaluate the new trial root
         else {
            if (xm > 0)
               b += Math.abs (tolerance1);
            else
               b += -Math.abs (tolerance1);
         }
         fb = integ (b) - ccc;
      }

      return b;
   }
   public String toString()
   {
      String desc = super.toString();
      desc += "tolerance : " + tolerance + "\n";
      return desc;
   }

}
