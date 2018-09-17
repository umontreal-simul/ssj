/*
 * Class:        InverseDistFromDensity
 * Description:  computing the inverse of an arbitrary continuous distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        June 2009
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
   import umontreal.ssj.functions.MathFunction;
import umontreal.ssj.util.Misc;
import umontreal.ssj.functions.MathFunctionUtil;

/**
 * Implements a method for computing the inverse of an *arbitrary continuous*
 * distribution function when only the probability density is known
 * @cite rDER09a&thinsp;. The cumulative probabilities (cdf) are pre-computed
 * by numerical quadrature of the density using Gauss-Lobatto integration
 * over suitably small intervals to satisfy the required precision, and these
 * values are kept in tables. Then the algorithm uses polynomial
 * interpolation over the tabulated values to get the inverse cdf. The user
 * can select the desired precision and the degree of the interpolating
 * polynomials.
 *
 * The algorithm may fail for some distributions for which the density
 * becomes infinite at a point (for ex. the Gamma and the Beta distributions
 * with @f$\alpha< 1@f$) if one chooses too high a precision (a too small
 * `eps`, for ex. @f$\epsilon\sim10^{-15}@f$). However, it should work also
 * for continuous densities with finite discontinuities.
 *
 * While the setup time for this class is relatively slow, the numerical
 * inversion is extremely fast and practically independent of the required
 * precision and of the specific distribution. For comparisons between the
 * times of standard inversion and inversion from this class as well as
 * comparisons between setup times, see the introduction in class
 * @ref umontreal.ssj.randvar.InverseFromDensityGen from package `randvar`.
 *
 * Thus if only a few inverses are needed, then using this class is not
 * efficient because of the slow set-up. But if one wants to call `inverseF`
 * thousands of times or more, then using this class will be very efficient.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_general
 */
public class InverseDistFromDensity extends ContinuousDistribution {
   private final boolean DEBUG = false;
   protected static final double HALF_PI = Math.PI/2.0;
   private final double epsc = 1.0e-5; // for small values of lc in tails
   private final double HUGE = Double.MAX_VALUE / 2.0;
   private double epsu0;      // initial u-resolution
   private double xc;         // mode, mean, or median of distribution
                              // xc is an x where the density is high
   MathFunction m_dens;       // probability density
   private String name;       // Name of class

   private final int K0 = 128;  // initial size of tables A, F, ...
   private int Kmax;          // final size of tables A, F, ... is (Kmax + 1)
   private double[] A;        // x-values
   private double[] F;        // corresponding u-values (the CDF)
   private double[][] X;      // interpolation x-values in [A[k], A[k+1]]
   private double[][] U;      // interpolation u-values in [A[k], A[k+1]]
   private double[][] C;      // interpolation coefficients in [A[k], A[k+1]]
   private int order;  // order of interpolation polynomial in each [A[k], A[k+1]]
   private int[] Index;       // for indexed search in F[k]
   private int Imax;          // final size of Index is (Imax + 1)

   private double bleft;      // computational left limit of density
   private double bright;     // computational right limit of density
   private double bl;         // left border of the computational domain
   private double br;         // right border of the computational domain
   private double llc;        // absolute local concavity in left tail
   private double rlc;        // absolute local concavity in right tail
/*
   For NormalDist(0,1), local concavity c = 1/z^2
   For CauchyDist(0,1), local concavity c = -1/2 + 1/(2z^2)
   For GammaDist(a,lam), local concavity c = (a-1)/(1 - a +lam*x)^2
*/
   private double lc1;        // left, c1 = f(p)/f'(p) or f(p) / (c*f'(p))
   private double lc2;        // left, c2 = |f'(p)|(1 + c) / (f(p)^2)
   private double lc3;        // left, c3 = c/(1 + c)
   private double rc1;        // right, c1 = f(p) / f'(p) or f(p) / (c*f'(p))
   private double rc2;        // right, c2 = |f'(p)|(1 + c) / (f(p)^2)
   private double rc3;        // right, c3 = c/(1 + c)
   private double epstail;    // = 0.05*epsu*I0;
   private boolean lcutF = false;    // cut-off flag, left tail
   private boolean rcutF = false;   // cut-off flag, right tail


   protected void printArray (double[] U) {
      System.out.print("      Tableau = (");
      for (int j = 0; j < U.length; j++)
         System.out.printf("  %f", U[j]);
      System.out.println("  )");
   }


   private class MaDensite implements MathFunction {
      private ContinuousDistribution cdist;

      public MaDensite(ContinuousDistribution dist) {
         cdist = dist;
         supportA = cdist.getXinf();
         supportB = cdist.getXsup();
     }

      public double evaluate (double x) {
         return cdist.density (x);
      }
   }


   private void init (double xc, double epsu, int n) {
      double[] zs = new double[n + 1];
      double[] ys = new double[n + 1];   // ksi[]
      double[] xs = new double[n + 1];
      double[] vs = new double[n + 1];
      double[] us = new double[n + 1];
      double[] cs = new double[n + 1];
      epsu = 0.9*epsu;
      findSupport(xc);

      double I0 = MathFunctionUtil.gaussLobatto (m_dens, bleft, bright, 1.0e-6);
      if (I0 > 1.05 || I0 < 0.95)
         throw new IllegalStateException("  NOT a probability density");
      epstail = 0.05*epsu*I0;
      epstail = Math.min(epstail, 1.e-10);
      epstail = Math.max(epstail, 1.e-15);
      double tol = epstail;
      findCutoff (bleft, epstail, false);    // left tail
      findCutoff (bright, epstail, true);    // right tail
      if (DEBUG)
         System.out.println("lcutF = " + lcutF + "\nrcutF = " + rcutF + "\n") ;

      reserve(0, n);
      A[0] = bl;
      if (lcutF)
         F[0] = epstail;
      else
         F[0] = 0;
      ys[0] = 0;
      final double HMIN = 1.0e-12;     // smallest integration step h
      double h = (br - bl) / K0;
      int j;
      int k = 0;
      calcChebyZ(zs, n);
      double eps = 0;
      if (DEBUG)
         System.out.println(
         "  k                 a_k                      F_k                h");

      while (A[k] < br) {
         while (h >= HMIN) {
            calcChebyX(zs, xs, n, h);
            calcU(m_dens, A[k], xs, us, n, tol);
            Misc.interpol(n, us, xs, cs);
            NTest(us, vs, n);
            // Evaluate Newton interpolating polynomial at vs[j].
            for (j = 1; j <= n; j++)
               ys[j] = Misc.evalPoly(n, us, cs, vs[j]);
            // NEval (cs, us, vs, ys, n);
            try {
               eps = calcEps(m_dens, A[k], ys, vs, n, tol);
            } catch (IllegalArgumentException e) {
               h = 0.5 * h;
               continue;
            }
            if (eps <= epsu)
               break;
            else
               h = 0.8 * h;
         }
         if (k + 1 >= A.length)
            reserve(k, n);
         copy (k, cs, us, xs, n);
         if (DEBUG)
            System.out.printf(
            " %d       %16.12f       %20.16g      %g%n", k, A[k], F[k], h);
         if (F[k] > 1.01)
            throw new IllegalStateException("Unable to compute CDF");
         k++;

         if (eps < epsu / 3.0)
            h = 1.3 * h;
         if (h < HMIN)
            h = HMIN;
         if (A[k] > br) {
            A[k] = br;
            F[k] = 1;
         }
      }

      if (DEBUG) {
         System.out.printf(
            " %d       %16.12f       %20.16g      %g%n",
            k, A[k], F[k], A[k] - A[k - 1]);
         System.out.println("\nFin du tableau");
      }
      Kmax = k;
      while (k > 0 && F[k] >= 1.) {
         F[k] = 1.;
         k--;
      }
      reserve(-Kmax, n);
      createIndex (Kmax);
   }

   /**
    * Given a continuous distribution `dist` with a well-defined density
    * method, this class will compute tables for the numerical inverse of
    * the distribution. The user may wish to set the left and the right
    * boundaries between which the density is non-zero by calling methods
    * umontreal.ssj.probdist.ContinuousDistribution.setXinf and
    * umontreal.ssj.probdist.ContinuousDistribution.setXsup of `dist`, for
    * better efficiency. Argument `xc` can be the mean, the mode or any
    * other @f$x@f$ for which the density is relatively large. The
    * @f$u@f$-resolution `eps` is the required absolute error in the cdf,
    * and `order` is the degree of the Newton interpolating polynomial
    * over each interval. An `order` of 3 or 5, and an `eps` of
    * @f$10^{-6}@f$ to @f$10^{-12}@f$ are usually good choices.
    * Restrictions: @f$3 \le\mathtt{order} \le12@f$.
    */
   public InverseDistFromDensity (ContinuousDistribution dist, double xc,
                                  double eps, int order) {
      setParams (dist, null, xc, eps, order);
      init (xc, eps, order);
   }

   /**
    * Given a continuous probability density `dens`, this class will
    * compute tables for the numerical inverse of the distribution. The
    * left and the right boundaries of the density are `xleft` and
    * `xright` (the density is 0 outside the interval <tt>[xleft,
    * xright]</tt>). See the description of the other constructor.
    */
   public InverseDistFromDensity (MathFunction dens, double xc, double eps,
                                  int order, double xleft, double xright) {
      supportA = xleft;
      supportB = xright;
      setParams (null, dens, xc, eps, order);
      init (xc, eps, order);
   }

   /**
    * Computes the probability density at @f$x@f$.
    */
   public double density (double x) {
      return m_dens.evaluate (x);
   }

   /**
    * Computes the distribution function at @f$x@f$.
    */
   public double cdf (double x) {
      throw new UnsupportedOperationException("cdf not implemented");
   }

   /**
    * Computes the inverse distribution function at @f$u@f$.
    */
   public double inverseF (double u) {
      if (u < 0.0 || u > 1.0)
          throw new IllegalArgumentException ("u not in [0,1]");
      if (u >= 1.0)
          return supportB;
      if (u <= 0.0)
          return supportA;
      if ((u < epstail) && lcutF)
         return uinvLeftTail (u);
      if ((u > 1.0 - epstail) && rcutF)
         return uinvRightTail (u);

      int k = searchIndex(u);
      double x = A[k] + Misc.evalPoly(order, U[k], C[k], u - F[k]);
      if (x <= supportA)
         return supportA;
      if (x >= supportB)
         return supportB;
      return x;
   }

   /**
    * Returns the `xc` given in the constructor.
    */
   public double getXc() {
      return xc;
   }

   /**
    * Returns the @f$u@f$-resolution `eps` associated with this object.
    */
   public double getEpsilon() {
      return epsu0;
   }

   /**
    * Returns the order associated with this object.
    */
   public int getOrder() {
      return order;
   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is returned as: [<tt>xc</tt>, `eps`,
    * <tt>order</tt>].
    */
   public double[] getParams() {
      double[] retour = {xc, epsu0, order};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString() {
      return name;
   }
   private void createIndex (int Kmax) {
      // create table for indexed search
      Imax = 2*Kmax;
      Index = new int [Imax + 1];
      Index[0] = 0;
      Index[Imax] = Kmax - 1;
      double u;
      int k = 1;
      for (int i = 1; i < Imax; i++) {
         u = (double) i / Imax;
         while (u >= F[k])
            k++;
         Index[i] = k-1;
      }
   }


   private int searchIndex (double u) {
      // search index of interval for interpolation in [F[k], F[k+1]]
      int i = (int) (Imax*u);
      int k = Index[i];
      while (u >= F[k]  && k < Kmax)
         k++;
      if (k <= 0)
         return 0;
      return k-1;
   }


   private void copy (int k, double[] cs, double[] us, double[] xs, int n) {
      // Copy parameters in interval [A[k], A[k+1]]
      for (int j = 0; j <= n; j++) {
         X[k][j] = xs[j];
         C[k][j] = cs[j];
         U[k][j] = us[j];
      }
      A[k+1] = A[k] + xs[n];
      F[k+1] = F[k] + us[n];
   }


   private double calcEps (MathFunction dens, double a, double[] Y,
                           double[] V, int n, double tol) {
      // Test if precision at test points Y is good enough
      // a is beginning of interval
      // Y are test points
      // V are values of CDF to compare with
      // n is order of interpolation
      // returns max eps
      // throw exception if Y[j] < Y[j-1] in gaussLobatto
      double eps = 0;
      double dif;
      double u = 0;
      for (int j = 1; j <= n; j++) {
         u += MathFunctionUtil.gaussLobatto (dens, a + Y[j-1], a + Y[j], tol);
         dif = Math.abs(u - V[j]);
         if (dif > eps)
            eps = dif;
      }
      return eps;
   }


   private void NEval (double[] C, double[] U, double[] T, double[] Y, int n) {
      // Evaluate Newton interpolating polynomial at T[j].
      // U are interpolation points
      // C are interpolation coefficients
      // Returns results in Y[j]
      int j;
      boolean fail = false;
      Y[0] = 0;
      for (j = 1; j <= n; j++) {
         Y[j] = Misc.evalPoly(n, U, C, T[j]);
         if (Y[j] < Y[j-1])
            fail = true;
      }

      if (fail) {
//       System.out.println("NEval");
         for (j = 1; j <= n; j++)
            Y[j] = Misc.evalPoly(1, U, C, T[j]);
      }
   }


   private void calcU (MathFunction dens, double a, double[] X, double[] U,
                       int n, double tol) {
      // compute CDF over n sub-intervals in [A[k], A[k+1]]
      // a is beginning of interval
      // X are x-values
      // U are values of CDF
      // precision is tol

      U[0] = 0;
      for (int j = 1; j <= n; j++)
         U[j] = U[j-1] +
            MathFunctionUtil.gaussLobatto (dens, a + X[j-1], a + X[j], tol);
   }


   private void reserve (int m, int n) {
      // Reserve memory for object
      A = reserve (A, m);
      F = reserve (F, m);
      C = reserve (C, m, n);
      U = reserve (U, m, n);
      X = reserve (X, m, n);
   }


   private double[] reserve (double[] T, int m) {
      if (m == 0) {
         // first call, just reserve memory.
         T = new double[K0 + 1];

      } else if (m < 0) {
         // Computation of table is complete. Table capacity is larger than
         // size: Resize table to exact size (-m + 1) and keep old values.
         m = -m;
         double[] tem = new double[m + 1];
         for (int i = 0; i <= m; i++)
            tem[i] = T[i];
         T = tem;

      } else {
         // Array is too short: reserve more memory and keep old values
         double[] tem = new double[2*m + 1];
         for (int i = 0; i <= m; i++)
            tem[i] = T[i];
         T = tem;
      }
      return T;
   }


   private double[][] reserve (double[][] T, int m, int n) {
      if (m == 0) {
         // first call, just reserve memory.
         T = new double[K0 + 1][n+1];

      } else if (m < 0) {
         // Computation of table is complete. Table capacity is larger than
         // size: Resize table to exact size (-m + 1) and keep old values.
         m = -m;
         double[][] tem = new double[m + 1][n+1];
         int j;
         for (int i = 0; i <= m; i++) {
            for (j = 0; j <= n; j++)
               tem[i][j] = T[i][j];
         }
         T = tem;

         } else {
         // Array is too short: reserve more memory and keep old values
         double[][] tem = new double[2*m + 1][n+1];
         int j;
         for (int i = 0; i <= m; i++) {
            for (j = 0; j <= n; j++)
               tem[i][j] = T[i][j];
         }
         T = tem;
      }
      return T;
   }


   private void NTest (double[] U, double[] T, int n) {
      // Routine 3 NTest in cite{rDER09a}
      // Compute test points T, given U
      int i, j, k;
      double s, sq, tem;
      T[0] = 0;
      for (k = 1; k <= n; k++) {
         T[k] = (U[k-1] + U[k]) / 2.;
         for (j = 0; j < 2 ; j++) {
            s = 0;
            sq = 0;
            for (i = 0; i <= n; i++) {
               tem = T[k] - U[i];
               if (tem == 0.)
                  break;
               tem = 1.0/tem;
               s += tem;
               sq += tem*tem;
            }
            if (sq != 0.)
               T[k] += s/sq;
         }
      }
   }


   private void calcChebyZ (double[] Z, int n) {
      // Eq. (3) in cite{rDER09a}. z_j = sin(j*phi)*sin((j+1)*phi)/cos(phi)
      // Compute normalized Chebyshev points in [0, 1]
      double phi = HALF_PI/(n+1);
      double c = Math.cos(phi);
      double y;
      double temp = 0;
      for (int j = 0; j < n; j++) {
         y = temp;
         temp = Math.sin((j+1)*phi);
         y *= temp;
         Z[j] = y/c;
      }
      Z[n] = 1;
   }


   private void calcChebyX (double[] Z, double[] X, int n, double h) {
      // Compute Chebyshev points in [0, h]
      for (int j = 1; j < n; j++)
         X[j] = h*Z[j];
      X[0] = 0;
      X[n] = h;
   }


   private double binSearch (double xa, double xb, double eps,
                                             boolean right) {
      /*
       * Binary search:
       *    find x such that   fa*epslow < f < fb*eps in the left tail
       *    find x such that   fa*eps > f > fb*epslow in the right tail
       *    where fa = density(xa), fb = density(xb),  f = density(x).
       * We find an x such that density(x) is a little smaller than eps
       */
      final double epslow = 0.1 * eps;
      double x = 0, y = 0;
      boolean fini = false;

      if (right) {    // right tail
         while (!fini) {
            x = 0.5 * (xa + xb);
            if ((xb - xa) < eps*Math.abs(x) || (xb - xa) < eps) {
               fini = true;
               if (x > supportB)
                  x = supportB;
            }
            y = m_dens.evaluate(x);
            if (y < epslow) {
               xb = x;
            } else if (y > eps) {
               xa = x;
            } else
               fini = true;
         }

      } else {   // left tail
         while (!fini) {
            x = 0.5 * (xa + xb);
            if ((xb - xa) < eps*Math.abs(x) || (xb - xa) < eps) {
               fini = true;
               if (x < supportA)
                  x = supportA;
            }
            y = m_dens.evaluate(x);
            if (y < epslow) {
               xa = x;
            } else if (y > eps) {
               xb = x;
            } else
               fini = true;
         }
      }
      if (DEBUG)
         System.out.printf(
         "binSearch   x =  %g    f =  %g     r =  %g%n", x, y, y/eps);

      return x;
   }


   private void findSupport (double xc) {
      /*
       * Find interval where density is non-negligible (above some epsilon):
       * find points bleft < xc < bright such that
       *      density(bleft) ~ density(bright) ~ 10^(-13)*density(xc)
       */
      boolean flagL = false;
      boolean flagR = false;
      final double DELTA = 1.0e-100;
      final double DELTAR = 1.0e-14;
      double x, y;
      double bl = supportA;
      double br = supportB;

      if (bl > Double.NEGATIVE_INFINITY) {
         // Density is 0 for finite x < bl
         y = m_dens.evaluate(bl);
         x = bl;
         if (y >= HUGE || y <= 0.0) {
            // density is infinite or 0 at bl; choose bl --> bl(1 + epsilon)
            x = bl + DELTAR * Math.abs(bl);
            if (x == 0)
               // bl is 0 --> choose bl = DELTA
               x = DELTA;
            y = m_dens.evaluate(x);
         }

         if (y >= HUGE)
            throw new UnsupportedOperationException
            ("Infinite density at left boundary");

         if (y >= 1.0e-50) {
            // f(bl) is large enough; we have found bl
            flagL = true;
            bl = x;
         }
      }

      if (br < Double.POSITIVE_INFINITY) {
         // Density is 0 for finite x > br
         y = m_dens.evaluate(br);
         x = br;
         if (y >= HUGE || y <= 0.0) {
            // density is infinite or 0 at br; choose br --> br(1 - epsilon)
            x = br - DELTAR * Math.abs(br);
            if (x == 0)
               // br is 0 --> choose br = -DELTA
               x = -DELTA;
            y = m_dens.evaluate(x);
         }

         if (y >= HUGE)
            throw new UnsupportedOperationException
            ("Infinite density at right boundary");

         if (y >= 1.0e-50) {
            // f(br) is large enough; we have found br
            flagR = true;
            br = x;
         }
      }

      bleft = bl;
      bright = br;
      if (flagL && flagR)
         return;

      // We have not found bl or br
      double h;
      y = m_dens.evaluate(xc);
      double epsy = 1.0e-13*y;
      double xa, xb;


      if (!flagR) {
         // Find br: start at xc; increase x until density is very small
         h = 1;
         xa = xc;
         xb = xc + h;
         while (m_dens.evaluate(xb) >= epsy) {
            xa = xb;
            h *= 2.0;
            xb += h;
         }
         // Now we have density(xa) > epsy > density(xb)

        if (xb > supportB) {
            // density = 0 outside [supportA, supportB]
            xb = supportB;
         }
         x = binSearch (xa, xb, epsy, true);
         bright = x;   // Have found br
    }

      if (!flagL) {
         h = 1;
         xb = xc;
         xa = xc - h;
         while (m_dens.evaluate(xa) >= epsy) {
            xb = xa;
            h *= 2.0;
            xa -= h;
         }
         // Now we have density(xa) < epsy < density(xb)

         if (xa < supportA) {
            // density = 0 outside [supportA, supportB]
            xa = supportA;
         }
         x = binSearch (xa, xb, epsy, false);
         bleft = x;   // Have found bl
     }
   }


   protected void setParams (ContinuousDistribution dist, MathFunction dens,
              double xc, double eps, int order) {
      // Sets the parameter of this object
      if (eps < 1.0e-15)
         throw new IllegalArgumentException ("eps < 10^{-15}");
      if (eps > 1.0e-3)
         throw new IllegalArgumentException ("eps > 10^{-3}");
      if (order < 3)
         throw new IllegalArgumentException ("order < 3");
      if (order > 12)
         throw new IllegalArgumentException ("order > 12");
      epsu0 = eps;
      this.xc = xc;
      this.order = order;

      StringBuffer sb = new StringBuffer ("InverseDistFromDensity: ");
      if (dist == null) {
         m_dens = dens;
      } else {
         m_dens = new MaDensite(dist);
         sb.append (dist.toString());
      }
      name = sb.toString();
   }


   private double uinvLeftTail (double u) {
      // Returns x = inverseF(u) in left tail

      double x = 0;
      if (llc <= epsc)
            x = bl + lc1 * Math.log (u*lc2);
      else
            x = bl + lc1 * (Math.pow (u*lc2, lc3) - 1.);
      if (x <= supportA)
         return supportA;
      return x;
   }


   private double uinvRightTail (double u) {
      // Returns x = inverseF(u) in right tail

      double x = 0;
      double v = 1. - u;
      if (rlc <= epsc)
            x = br + rc1 * Math.log (v*rc2);
      else
            x = br + rc1 * (Math.pow (v*rc2, rc3) - 1.);
      if (x >= supportB)
         return supportB;
      return x;
   }


   private void findCutoff (double x0, double eps, boolean right) {
      /*
       * Find cut-off points for the computational domain.
       * Find cut-off x in the tails such that cdf(x) = eps in the left
       *    tail, and eps = 1 - cdf(x) in the right tail.
       * Uses successive approximations starting at x = x0.
       * If right is true, case of the right tail; otherwise the left tail.
       * The program uses T_c-concavity of densities as described in
       * Leydold et al.
       */
      final double epsx = 1.0e-3;
      final double range = bright - bleft;
      double del;

      if (right) {
          del = m_dens.evaluate(bright) - m_dens.evaluate(bright - epsx);
          if ((supportB < Double.POSITIVE_INFINITY) &&
                (supportB - bright <= epsx ||
                 supportB - bright <= Math.abs(supportB)*epsx)) {
            // If density is non-negligible right up to domain limit supportB,
            // then cutoff is bright. There is no right tail. We want cutoff
            // at bright in case density(supportB) = infinite.
            if (del < 0) {
               // density decreases toward supportB;
               br = supportB;
            } else {
               // density increases toward supportB; may be infinite
               br = bright;
            }
            rcutF = false;
            return;
         } else {
            rcutF = true;   // There is a right tail
         }

      } else {
         del = m_dens.evaluate(bleft + epsx) - m_dens.evaluate(bleft);
         if ((supportA > Double.NEGATIVE_INFINITY) &&
             (bleft - supportA <= epsx ||
              bleft - supportA <= Math.abs(supportA)*epsx)) {
            // If density is non-negligible right down to domain limit supportA,
            // then cutoff is bleft. There is no left tail. We want cutoff
            // at bleft in case density(supportA) = infinite.
            if (del > 0) {
               // density decreases toward supportA
               bl = supportA;
            } else {
               // density increases toward supportA; may be infinite
               bl = bleft;
            }
            lcutF = false;
            return;
         } else {
            lcutF = true;   // There is a left tail
         }
      }

      double c = 0;
      double h = 1.0/64.0;      // step to compute derivative
      h = Math.max (h, (bright - bleft) / (1024));
      double x = x0, xnew;
      double y = 0, yl = 0, yr = 0, yprime = 0;
      double tem = 0;
      int iter = 0;
      final int ITERMAX = 30;
      boolean fini = false;

      while (!fini && iter < ITERMAX) {
         iter++;
         boolean ended = false;
         int it = 0;

         while (!ended && it < 10) {
            it++;
            if (x + h > supportB)
               h = supportB - x;
            if (x - h < supportA)
               h = x - supportA;
            yr = m_dens.evaluate(x + h);
            y = m_dens.evaluate(x);
            yl = m_dens.evaluate(x - h);
            if (!(yl == 0 || yr == 0 || y == 0))
               ended = true;
            else
               h /= 2;
         }

         c = yr / (yr - y) + yl / (yl - y) - 1.;  // the local concavity lc
         yprime = (yr - yl) / (2. * h);          // first derivative
         tem = Math.abs (y * y / ((c + 1.) * yprime));  // tail area of CDF

         if (Double.isNaN (tem))
            break;
         if (Math.abs (tem / eps - 1.) < 1.e-4)   // accuracy is good?
            break;
         if (Math.abs(c) <= epsc) {
            tem = eps * Math.abs(yprime) / (y * y);   // formula (10)
            if (tem <= 0)
               break;
            xnew = x + y / yprime * Math.log(tem);
         } else {
            tem = (1. + c) * eps * Math.abs(yprime) / (y * y); // formula(10)
            if (tem < 0)
               break;
            xnew = x + y / (c*yprime) * (Math.pow(tem, c / (1. + c)) - 1.);
         }

         if (DEBUG)
            System.out.printf(
            "Cutoff   x =  %g    y =  %g     c =  %g%n", xnew, y, c);

         if ((Math.abs(xnew - x) <= Math.abs(x)*epsx) ||
             (Math.abs(xnew - x) <= epsx))
            fini = true;     // found cut-off x
         else
            x = xnew;

         // Given good x, precompute some parameters in formula (10)
         if (right) {
            rlc = Math.abs(c);
            br = x;
            rc3 = c / (1 + c);
            rc2 = tem / eps;
            rc1 = y / yprime;
            if (Math.abs(c) > epsc)
               rc1 /= c;

         } else {
            llc = Math.abs(c);
            bl = x;
            lc3 = c / (1 + c);
            lc2 = tem / eps;
            lc1 = y / yprime;
            if (Math.abs(c) > epsc)
               lc1 /= c;
         }

       if (Math.abs(xnew - x) >= range)
           fini = true;
      }

      if (right) {
         if ((rc1 == 0 && rc2 == 0 && rc3 == 0)) {
            br = bright;
            rcutF = false;
         }
      } else {
          if ((lc1 == 0 && lc2 == 0 && lc3 == 0)) {
            bl = bleft;
            lcutF = false;
         }
      }
   }

}