/*
 * Class:        PoissonTIACGen
 * Description:  random variate generators having the Poisson distribution 
                 using the tabulated inversion combined with the acceptance
                 complement method
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
package umontreal.ssj.randvar;
import umontreal.ssj.probdist.*;
import umontreal.ssj.rng.*;

/**
 * This class implements random variate generators having the *Poisson*
 * distribution (see  @ref PoissonGen ). Uses the tabulated inversion
 * combined with the acceptance complement (*TIAC*) method of
 * @cite rAHR82b&thinsp;. The implementation is adapted from UNURAN
 * @cite iLEY02a&thinsp;.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_discrete
 */
public class PoissonTIACGen extends PoissonGen {
 
   private double[] pp    = new double[36];
   private int[]    llref = {0};
   private static double[] staticPP    = new double[36];
   private static int[]    staticllref = {0};
   // Used by TIAC, avoid creating a table upon each call.

   /**
    * Creates a Poisson random variate generator with parameter
    * @f$\lambda= @f$ `lambda`, using stream `s`.
    */
   public PoissonTIACGen (RandomStream s, double lambda) {
      super (s, null);
      init (lambda);
   }

   /**
    * Creates a new random variate generator using the Poisson
    * distribution `dist` and stream `s`.
    */
   public PoissonTIACGen (RandomStream s, PoissonDist dist) {
      super (s, dist);
      init (dist.getLambda ());
   }
public int nextInt() {
      return tiac (stream, lambda, pp, llref);
   }

/**
 * A static method for generating a random variate from a *Poisson*
 * distribution with parameter @f$\lambda@f$ = `lambda`.
 */
public static int nextInt (RandomStream s, double lambda) {
      return tiac (s, lambda, staticPP, staticllref);
   }
 

   private void init (double lam) {
      setParams (lam);
      for (int i = 0; i < pp.length; i++)
         pp[i] = 0.0;
   }


/* **************************************************************************
 * GENERATION METHOD : Tabulated Inversion combined with                     *
 *                     Acceptance Complement                                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * METHOD :   -  samples a random number from the Poisson distribution with  *
 *               parameter lambda > 0.                                       *
 *               Tabulated Inversion for  lambda < 10                        *
 *               Acceptance Complement for lambda >= 10.                     *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * CODE REFERENCE : The code is adapted from UNURAN                          *
 * UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * REFERENCE: - J.H. Ahrens, U. Dieter (1982): Computer generation of        * 
 *              Poisson deviates from modified normal distributions,         *
 *              ACM Trans. Math. Software 8, 163-179.                        *
 *                                                                           *
 * Implemented by R. Kremer, August 1990                                     *
 * Revised by E. Stadlober, April 1992                                       *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************
 *  UNURAN :  
 ****************************************************************************/
   private static int tiac (RandomStream s, double lambda, 
                            double[] pp, int[] llref) {
      final double  a0 = -0.5000000002;
      final double  a1 =  0.3333333343;
      final double  a2 = -0.2499998565;
      final double  a3 =  0.1999997049;
      final double  a4 = -0.1666848753;
      final double  a5 =  0.1428833286;
      final double  a6 = -0.1241963125;
      final double  a7 =  0.1101687109;
      final double  a8 = -0.1142650302;
      final double  a9 =  0.1055093006;
      // factorial for 0 <= k <= 9
      final int fac[] = {1,1,2,6,24,120,720,5040,40320,362880};

      double u;
      int i;
      if (lambda < 10) {
         int m = lambda > 1 ? (int)lambda : 1;
         int ll = llref[0];
         double p0, q, p;
         p0 = q = p = Math.exp (-lambda);
         int k = 0;       
         while (true) {
            u = s.nextDouble();     // Step u. Uniform sample 
            k = 0;
            if (u <= p0) 
               return k;

            // Step T. Table comparison 
            if (ll != 0) {               
               i = (u > 0.458) ? Math.min(ll,m) : 1;
               for (k = i; k <=ll; k++)
               if (u <= pp[k]) return k;
               if (ll == 35) continue;
            }

            // Step C. Creation of new prob. 
            for (k = ll +1; k <= 35; k++) {
                p *= lambda / (double)k;
                q += p;
                pp[k] = q;
                if (u <= q) {
                   ll = k;
                   llref[0] = ll;
                   return k;
                }
            }
            ll = 35;
            llref[0] = ll;
         }
      }
      else { // lambda > 10, we use  acceptance complement
         double sl = Math.sqrt (lambda);
         double d = 6*lambda*lambda;
         int l = (int)(lambda - 1.1484);

         // Step P. Preparations for steps Q and H
         double omega = 0.3989423 / sl;
         double b1 = 0.416666666667e-1/lambda;
         double b2 = 0.3*b1*b1;
         double c3 = 0.1428571*b1*b2;
         double c2 = b2 - 15.0*c3;
         double c1 = b1 - 6.0*b2 + 45.0*c3;
         double c0 = 1.0 - b1 + 3.0*b2 - 15.0*c3;
         double c = 0.1069/lambda;

         double t,g,lambda_k;
         double gx,gy,px,py,x,xx,delta,v;
         int sign;

         double e;
         int k;

         // Step N. Normal sample 
         // We don't use NormalGen because this could create bad side effects.
         // With NormalGen, the behavior of this method would depend
         // on the selected static method of NormalGen.
         t = NormalDist.inverseF (0, 1, s.nextDouble());
         g = lambda + sl*t;

         if (g >= 0) {
            k = (int)g;
            // Step I. Immediate acceptance
            if (k >= l) 
               return k;
            // Step S. Squeeze acceptance
            u = s.nextDouble();
            lambda_k = lambda - k;
            if (d*u >= lambda_k*lambda_k*lambda_k)
               return k;

            // FUNCTION F 
            if (k < 10) {
               px = -lambda;
               py = Math.exp (k*Math.log (lambda))/fac[k];
            }
            else {  // k >= 10
               delta = 0.83333333333e-1/(double)k;
               delta = delta - 4.8*delta*delta*delta*(1.-1./(3.5*k*k));
               v = (lambda_k)/(double)k;
               if (Math.abs (v) > 0.25)
                  px = k*Math.log (1. + v) - lambda_k - delta;
               else {
                  px = k*v*v;
                  px *= ((((((((a9*v+a8)*v+a7)*v+a6)*v+a5)*v+
                     a4)*v+a3)*v+a2)*v+a1)*v+a0;
                  px -= delta;
               }
               py = 0.3989422804 / Math.sqrt((double)k);
            }
            x = (0.5 - lambda_k)/sl;
            xx = x*x;
            gx = -0.5*xx;
            gy = omega*(((c3*xx + c2)*xx + c1)*xx + c0);
            // end FUNCTION F

            // Step Q. Quotient acceptance 
            if (gy*(1.0 - u) <= py*Math.exp (px - gx))
               return k;
         }

         // Step E. Double exponential sample
         while (true) {
            do {
               e = - Math.log (s.nextDouble());
               u = s.nextDouble();
               u = u + u - 1;
               sign = u < 0 ? -1 : 1;
               t = 1.8 + e*sign;
            } while (t <= -0.6744);
            k = (int)(lambda + sl*t);
            lambda_k = lambda - k;

            // FUNCTION F
            if (k < 10) {
               px = -lambda;
               py = Math.exp(k*Math.log (lambda))/fac[k];
            }
            else { // k >= 10
               delta = 0.83333333333e-1/(double)k;
               delta = delta - 4.8*delta*delta*delta*(1.0-1.0/(3.5*k*k));
               v = lambda_k/(double)k;
               if (Math.abs (v) > 0.25)
                  px = k*Math.log (1. + v) - lambda_k - delta;
               else {
                  px = k*v*v;
                  px *= ((((((((a9*v+a8)*v+a7)*v+a6)*v+a5)*v+
                           a4)*v+a3)*v+a2)*v+a1)*v+a0;
                  px -= delta;
               }
               py = 0.3989422804/Math.sqrt((double)k);
            }
            x = (0.5 - lambda_k)/sl;
            xx = x*x;
            gx = -0.5*xx;
            gy = omega*(((c3*xx + c2)*xx + c1)*xx + c0);
            // end FUNCTION F

            // Step H. Hat acceptance
            if (c*sign*u <= py*Math.exp (px + e) - gy*Math.exp (gx + e)) 
               return k;
         }
      }
   }

   static {
      for (int i = 0; i < staticPP.length; i++)
         staticPP[i] = 0.0;
   }
}