/*
 * Class:        NI3
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
import umontreal.ssj.probdist.*;

/**
 * Extends the class  @ref NortaInitDisc and implements the algorithm NI3. It
 * uses the function @f$g_r@f$ and its derivative @f$g’_r@f$, so it calls the
 * methods  #integ and  #deriv, given in (
 * {@link REF_probdistmulti_norta_NortaInitDisc_gr_M gr_M} )
 * and ( {@link REF_probdistmulti_norta_NortaInitDisc_grp_M
 * grp_M} ) and uses an adapted version of the Newton-Raphson algorithm. At
 * any iteration, if the solution falls outside the search interval, the
 * algorithm uses bisection and halves the interval length to guarantee
 * convergence. The initial solution is taken as @f$\rho_0=2 \sin(\pi r_X
 * /6)@f$, and then at each iteration @f$k@f$, @f$f_r(\rho_k)@f$ and
 * @f$f_r^’(\rho_k)@f$ are calculated and a solution is computed by the
 * recurrence formula:
 * @f[
 *   \rho_{k+1}=\rho_k-\frac{f_r(\rho_k)}{f_r^’(\rho_k)}.
 * @f]
 * The algorithm stops at iteration @f$k@f$ if
 * @f$|\rho_{k-1}-\rho_k|\leq\epsilon@f$. The function @f$f_r@f$ is the
 * one given in ( {@link REF_probdistmulti_norta_overview_fr
 * fr} ).
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class NI3 extends NortaInitDisc
{
   private double tolerance; /* Desired accuracy for the root-finder
   			        algorithm (epsilon in paragraph
   			        "Method NI3" of section 3 in paper).*/

   /**
    * Constructor with the target rank correlation `rX`, the two discrete
    * marginals `dist1` and `dist2`, the parameter for the truncation `tr`
    * (see the constructor of class  @ref NortaInitDisc ), and the
    * specific parameter @f$\epsilon=@f$ `tolerance` for the algorithm
    * NI3, as defined above.
    */
   public NI3 (double rX,
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
    * Computes and returns the correlation @f$\rho_Z@f$ using algorithm
    * NI3.
    */
   public double computeCorr()
   {
      final double ITMAX = 100; // Maximum number of iterations.
      double xl, xh;  /* Left and right endpoints of the root bracket at
			 all iterations. */
      double b = 0.0; // The returned solution.
      double f, df;   // Function and its derivative evaluations.
      double dx;      // Correction term.
      /** f, df, dx correspond to f(rho_k), f'(rho_k) and f(rho_k)/f'(rho_k),
          respectively, in the paper (paragraph "Method NI3"of section 3). */
      double dxold; // The root correction at one iteration before the last one.
      double temp;// The root at one iteration before the last one.
      double ccc = rX * sd1 * sd2 + mu1 * mu2; // Precompute constant.

      if (rX == 0.0)
         return 0.0;
      if (rX > 0.0) {              // Orient the search
         xl = 0.0;
         xh = 1.0;
      } else {
         xl = -1.0;
         xh = 0.0;
      }

      b = 2 * Math.sin (Math.PI * rX / 6); // Initial guess
      dxold = xh - xl;
      dx = dxold;
      f = integ (b) - ccc;
      df = deriv (b);

      for (int i = 1; i <= ITMAX; i++) { // Begin the search
         if ((((b - xh) * df - f) * ((b - xl) * df - f) > 0.0)
               || (Math.abs (2.0 * f) > Math.abs (dxold * df))) {
            // Do bisection if solution is out of range
            // or not decreasing fast enough
            dxold = dx;
            dx = 0.5 * (xh - xl);
            b = xl + dx;
            if (xl == b)      // Change in root is negligible.
               return b;      // Accept this root
         } else {
            dxold = dx;
            dx = f / df;
            temp = b;
            b -= dx;
            if (temp == b)
               return b;
         }
         if (Math.abs (dx) < tolerance)
            return b;          // Convergence check
         f = integ (b) - ccc;
         df = deriv (b);
         if (f < 0.0)
            xl = b;            // Maintain the brackets on the root
         else
            xh = b;
      }
      return b;
   }
   public String toString()
   {
    // To display the inputs.
      String desc = super.toString();
      desc += "tolerance : " + tolerance + "\n";
      return desc;
   }

}
