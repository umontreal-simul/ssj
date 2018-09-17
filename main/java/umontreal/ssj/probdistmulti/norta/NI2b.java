/*
 * Class:        NI2b
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
 * Extends the class  @ref NortaInitDisc and implements the algorithm NI2b.
 * It is a variant of NI2a. It uses the derivative, so it calls the method
 * #deriv to compute the function @f$g’_r@f$ given in (
 * {@link REF_probdistmulti_norta_NortaInitDisc_grp_M grp_M} )
 * and uses numerical integration with Simpson’s rules as well. But the
 * integration grid is either in the interval @f$[0,1-\delta]@f$ or
 * @f$[-1+\delta, 0]@f$, depending on the sign of @f$r_X@f$. Here the number
 * of subintervals of integration is fixed to @f$m@f$ and the algorithm stops
 * at iteration @f$k@f$ if the root is in subinterval
 * @f$[\rho_{k-1},\rho_k]@f$, and a quadratic interpolation is used to
 * compute the solution. For this, the method  #interpol of class
 * @ref umontreal.ssj.util.Misc (from package  util of SSJ
 * @cite iLEC04j&thinsp;) is used.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class NI2b extends NortaInitDisc
{
   private int m; /* Number of subintervals for the integration = max.
   		      number of iterations (also named m in the paper,
   		      paragraph "Method NI2" of section 3).*/
   private double delta; /* Small positive parameter to make sure that
                                rho_m is not too close to 1 or -1;
   			     (also named delta in the paper, paragraph
   			     "Method NI2" of section 3)*/

   /**
    * Constructor with the target rank correlation `rX`, the two discrete
    * marginals `dist1` and `dist2`, the parameter for the truncation `tr`
    * (see the constructor of class  @ref NortaInitDisc ), and the
    * specific parameters `m` and @f$\delta=@f$ `delta` for the algorithm
    * NI2b, as described above.
    */
   public NI2b(double rX,              
               DiscreteDistributionInt dist1, 
               DiscreteDistributionInt dist2,
               double tr, 
               int m, 
               double delta)
   {
      super(rX, dist1, dist2, tr);
      this.m = m;
      this.delta = delta;
      computeParams();
   }

   /**
    * Computes and returns the correlation @f$\rho_Z@f$ using the
    * algorithm NI2b.
    */
   public double computeCorr()
   {
      // x, y and c oefficients used for quadratic interpolation
      double[] x = new double[3];
      double[] y = new double[3];
      double[] c = new double[3];
      double xtemp = 0.0, temp = 0.0; /* Values of rho and the recursive
      					   quantity I_k, given in paragraph
      					   "Method NI2" of section 3 in
      					   the paper,2 iterations before the
      					   last one. */
      double xold, iold; /* Values of rho and I_k at one iteration before
      			      the last one.*/
      double xnew, inew; // rho and I_k at the last iteration.
      double dold, dmid, dnew; /* Values of the derivative function g'
      				   at points xold, xmid (xold+h) and xnew.
      				   They correspond to g'(rho_0+2kh-2h),
      				   g'(rho_0+2kh-h) and g'(rho_0+2kh) in the
      				   formula of I_k in the paper
      				   (paragraph "Method NI2" of section 3). */
      double h = (1 - delta) / (2 * m); /* Step size for the
      					 integration-grid spacing (2*h).
      					 It corresponds to h given in
      					 the third paragraph of
      					 section 4 in the paper */
      double b = 0.0; // The returned solution.
      double rho1 = 0.0; // The initial guess.
      double intg1 = mu1 * mu2; // Computes g_r(rho1).
      double gr = rX * sd1 * sd2 + mu1 * mu2; /* Target value; integ(rho)
      						   = gr is equivalent to
      						   rho = the solution. */
      // Pre-compute constants
      double hd3 = h / 3;
      double h2 = 2 * h;

      if (intg1 == gr)
         return rho1;

      if (0 < rX && rX < 1) { // Do search between 0 and rho_m=1-delta
         xold = rho1;
         dold = deriv(xold);
         iold = intg1;

         for (int i = 1; i <= m ;i++) { // Begin the search
            dmid = deriv(xold + h);
            xnew = xold + h2;
            dnew = deriv(xnew);
            inew = iold + hd3 * (dold + 4 * dmid + dnew);
            if (inew >= gr) { // The root is in current bracketing interval
               // Compute the parameters of quadratic interpolation
               x[0] = xtemp;
               x[1] = xold;
               x[2] = xnew;
               y[0] = temp;
               y[1] = iold;
               y[2] = inew;
               Misc.interpol(2, x, y, c);
               b = (c[2] * (xtemp + xold) - c[1] + Math.sqrt((c[1]
                    - c[2] * (xtemp + xold)) * (c[1] - c[2]
                    * (xtemp + xold)) - 4 * c[2] * (c[0] - c[1]
                    * xtemp + c[2] * xtemp * xold - gr))) / (2 * c[2]);
               return b;
            }
            xtemp = xold;
            temp = iold;
            xold = xnew;
            dold = dnew ;
            iold = inew;
         }
         // Integration up to 1-delta did not bracket root
         // return 1-delta/2 ( = midpoint of current bracketing interval)
         b = 1 - delta / 2;
      }

      if ( -1 < rX && rX < 0) { // Do search between rho_m=-1+delta and 0
         xold = rho1;
         dold = deriv(xold);
         iold = intg1;

         for (int i = 1; i <= m ;i++) { // Begin the search
            dmid = deriv(xold - h);
            xnew = xold - h2;
            dnew = deriv(xnew);
            inew = iold - hd3 * (dold + 4 * dmid + dnew);
            if (inew <= gr) { // The root is in current bracketing interval
               // Compute the parameters of quadratic interpolation
               x[0] = xnew;
               x[1] = xold;
               x[2] = xtemp;
               y[0] = inew;
               y[1] = iold;
               y[2] = temp;
               Misc.interpol(2, x, y, c);
               b = (c[2] * (xnew + xold) - c[1] + Math.sqrt((c[1]
                    - c[2] * (xnew + xold)) * (c[1] - c[2]
                    * (xnew + xold)) - 4 * c[2] * (c[0] - c[1]
                    * xnew + c[2] * xnew * xold - gr))) / (2 * c[2]);
               return b;
            }
            xtemp = xold;
            temp = iold;
            xold = xnew;
            dold = dnew ;
            iold = inew;
         }
         // Integration up to -1+delta did not bracket root
         b = -1 + delta / 2; // return 1-delta/2
         // ( = midpoint of current bracketing interval )
      }
      return b;
   }
   public String toString()
   {
      String desc = super.toString();
      desc += "m :  " + m + "\n";
      desc += "delta : " + delta + "\n";
      return desc;
   }

}
