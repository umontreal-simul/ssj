/*
 * Class:        KolmogorovSmirnovDist
 * Description:  Kolmogorov-Smirnov distribution
 * Environment:  Java
 * Software:     SSJ
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
import umontreal.ssj.util.*;
import umontreal.ssj.functions.MathFunction;

/**
 * Extends the class  @ref ContinuousDistribution for the
 * *Kolmogorov-Smirnov* distribution with parameter @f$n@f$
 * @cite tDUR73a&thinsp;. Given an empirical distribution @f$F_n@f$ with
 * @f$n@f$ independent observations and a continuous distribution @f$F(x)@f$,
 * the two-sided *Kolmogorov–Smirnov* statistic is defined as
 * @f[
 *   D_n = \sup_{-\infty\le x \le\infty} |F_n(x) - F(x)|  =  \max\{D_n^+, D_n^-\},
 * @f]
 * where @f$D_n^+@f$ and @f$D_n^-@f$ are the *Kolmogorov–Smirnov*+ and
 * *Kolmogorov–Smirnov*@f$-@f$ statistics as defined in equations
 * {@link REF_probdist_KolmogorovSmirnovPlusDist_eq_DNp DNp}
 * and  {@link REF_probdist_KolmogorovSmirnovPlusDist_eq_DNm
 * DNm} of this guide. This class implements a high
 * precision version of the *Kolmogorov–Smirnov* distribution @f$P[D_n
 * \le x]@f$; it is a Java translation of the @f$C@f$ program written in
 * @cite tMAR03a&thinsp;. According to its authors, it should give 13 decimal
 * digits of precision. It is extremely slow for large values of @f$n@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_edf
 */
public class KolmogorovSmirnovDist extends ContinuousDistribution {
   protected int n;
   protected static final int NEXACT = 500;

   // For the Durbin matrix algorithm
   private static final double NORM = 1.0e140;
   private static final double INORM = 1.0e-140;
   private static final int LOGNORM = 140;



   //========================================================================

   private static class Function implements MathFunction {
      protected int n;
      protected double u;

      public Function (int n, double u) {
         this.n = n;
         this.u = u;
      }

      public double evaluate (double x) {
         return u - cdf(n,x);
      }
   }

   /**
    * Constructs a *Kolmogorov–Smirnov* distribution with parameter
    * @f$n@f$. Restriction: @f$n \ge1@f$.
    */
   public KolmogorovSmirnovDist (int n) {
      setN (n);
   }


   public double density (double x) {
      return density (n, x);
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

   private static double dclem (int n, double x, double EPS) {
      return (cdf(n, x + EPS) - cdf(n, x - EPS)) / (2.0 * EPS);
   }

   protected static double densConnue (int n, double x) {
      if ((x >= 1.0) || (x <= 0.5 / n))
         return 0.0;
      if (n == 1)
         return 2.0;

      if (x <= 1.0 / n) {
         double w;
         double t = 2.0 * x * n - 1.0;
         if (n <= NEXACT) {
            w = 2.0 * n * n * Num.factoPow (n);
            w *= Math.pow (t, (double) (n - 1));
            return w;
         }
         w = Num.lnFactorial  (n) + (n-1) * Math.log (t/n);
         return 2*n*Math.exp (w);
      }

      if (x >= 1.0 - 1.0 / n)
         return 2.0 * n * Math.pow (1.0 - x, (double) (n - 1));

      return -1.0;
   }

/**
 * Computes the density for the *Kolmogorov–Smirnov* distribution with
 * parameter @f$n@f$.
 */
public static double density (int n, double x) {
      double Res = densConnue(n,x);
      if (Res != -1.0)
         return Res;

      double EPS = 1.0 / 200.0;
      final double D1 = dclem(n, x, EPS);
      final double D2 = dclem(n, x, 2.0 * EPS);
      Res = D1 + (D1 - D2) / 3.0;
      if (Res <= 0.0)
         return 0.0;
      return Res;
   }


   /*=========================================================================

   The following implements the Durbin matrix algorithm and was programmed
   in C by G. Marsaglia, Wai Wan Tsang and Jingbo Wong in C.

   I have translated their program in Java. Only small modifications have
   been made in their program; the most important is to prevent the return
   of NAN or infinite values in some regions. (Richard Simard)

   =========================================================================*/

   /*
    The C program to compute Kolmogorov's distribution

                K(n,d) = Prob(D_n < d),         where

         D_n = max(x_1-0/n,x_2-1/n...,x_n-(n-1)/n,1/n-x_1,2/n-x_2,...,n/n-x_n)

       with  x_1<x_2,...<x_n  a purported set of n independent uniform [0,1)
       random variables sorted into increasing order.
       See G. Marsaglia, Wai Wan Tsang and Jingbo Wong,
          J.Stat.Software, 8, 18, pp 1--4, (2003).
   */


   private static void mMultiply (double[] A, double[] B, double[] C, int m){
      int i, j, k;
      double s;
      for (i = 0; i < m; i++)
         for (j = 0; j < m; j++) {
            s = 0.0;
            for (k = 0; k < m; k++)
               s += A[i * m + k] * B[k * m + j];
            C[i * m + j] = s;
         }
   }


   private static void renormalize (double[] V, int m, int[] p)
   {
      int i;
      for (i = 0; i < m * m; i++)
         V[i] *= INORM;
      p[0] += LOGNORM;
   }


   private static void mPower (double[] A, int eA, double[] V, int[] eV,
                                              int m, int n)
   {
      int i;
      if (n == 1) {
         for (i = 0; i < m * m; i++)
            V[i] = A[i];
         eV[0] = eA;
         return;
      }
      mPower (A, eA, V, eV, m, n / 2);

      double[] B = new double[m * m];
      int[] pB = new int[1];

      mMultiply (V, V, B, m);
      pB[0] = 2 * (eV[0]);
      if (B[(m / 2) * m + (m / 2)] > NORM)
         renormalize (B, m, pB);

      if (n % 2 == 0) {
         for (i = 0; i < m * m; i++)
            V[i] = B[i];
         eV[0] = pB[0];
      } else {
         mMultiply (A, B, V, m);
         eV[0] = eA + pB[0];
      }

      if (V[(m / 2) * m + (m / 2)] > NORM)
         renormalize (V, m, eV);
   }


   protected static double DurbinMatrix (int n, double d)
   {
      int k, m, i, j, g, eH;
      double h, s;
      double[] H;
      double[] Q;
      int[] pQ;

      //Omit next two lines if you require >7 digit accuracy in the right tail
      if (false) {
         s = d * d * n;
         if (s > 7.24 || (s > 3.76 && n > 99))
            return 1 - 2 * Math.exp (-(2.000071 + .331 / Math.sqrt (n) +
                    1.409 / n) * s);
      }
      k = (int) (n * d) + 1;
      m = 2 * k - 1;
      h = k - n * d;
      H = new double[m * m];
      Q = new double[m * m];
      pQ = new int[1];

      for (i = 0; i < m; i++)
         for (j = 0; j < m; j++)
            if (i - j + 1 < 0)
               H[i * m + j] = 0;
            else
               H[i * m + j] = 1;
      for (i = 0; i < m; i++) {
         H[i * m] -= Math.pow (h, (double)(i + 1));
         H[(m - 1) * m + i] -= Math.pow (h, (double)(m - i));
      }
      H[(m - 1) * m] += (2 * h - 1 > 0 ? Math.pow (2 * h - 1, (double) m) : 0);
      for (i = 0; i < m; i++)
         for (j = 0; j < m; j++)
            if (i - j + 1 > 0)
               for (g = 1; g <= i - j + 1; g++)
                  H[i * m + j] /= g;
      eH = 0;
      mPower (H, eH, Q, pQ, m, n);
      s = Q[(k - 1) * m + k - 1];

      for (i = 1; i <= n; i++) {
         s = s * (double) i / n;
         if (s < INORM) {
            s *= NORM;
            pQ[0] -= LOGNORM;
         }
      }
      s *= Math.pow (10., (double) pQ[0]);
      return s;
   }


   protected static double cdfConnu (int n, double x) {
      // For nx^2 > 18, barF(n, x) is smaller than 5e-16
      if ((n * x * x >= 18.0) || (x >= 1.0))
         return 1.0;

      if (x <= 0.5 / n)
         return 0.0;

      if (n == 1)
         return 2.0 * x - 1.0;

      if (x <= 1.0 / n) {
         double w;
         double t = 2.0 * x * n - 1.0;
         if (n <= NEXACT) {
            w = Num.factoPow (n);
            return w * Math.pow (t, (double) n);
         }
         w = Num.lnFactorial(n) + n * Math.log (t/n);
         return Math.exp (w);
      }

      if (x >= 1.0 - 1.0 / n) {
         return 1.0 - 2.0 * Math.pow (1.0 - x, (double) n);
      }

      return -1.0;
   }

/**
 * Computes the *Kolmogorov–Smirnov* distribution function @f$F(x)@f$ with
 * parameter @f$n@f$ using Durbin’s matrix formula @cite tDUR73a&thinsp;. It
 * is a translation of the @f$C@f$ program in @cite tMAR03a&thinsp;;
 * according to its authors, it returns 13 decimal digits of precision. It is
 * extremely slow for large @f$n@f$.
 */
public static double cdf (int n, double x) {
      double Res = cdfConnu(n,x);
      if (Res != -1.0)
         return Res;

      return DurbinMatrix (n, x);
   }


   protected static double barFConnu (int n, double x) {
      final double w = n * x * x;

      if ((w >= 370.0) || (x >= 1.0))
         return 0.0;
      if ((w <= 0.0274) || (x <= 0.5 / n))
         return 1.0;
      if (n == 1)
         return 2.0 - 2.0 * x;

      if (x <= 1.0 / n) {
         double v;
         final double t = 2.0 * x*n - 1.0;
         if (n <= NEXACT) {
            v = Num.factoPow (n);
            return 1.0 - v * Math.pow (t, (double) n);
         }
         v = Num.lnFactorial(n) + n * Math.log (t/n);
         return 1.0 - Math.exp (v);
      }

      if (x >= 1.0 - 1.0 / n) {
         return 2.0 * Math.pow (1.0 - x, (double) n);
      }

      return -1.0;
   }

/**
 * Computes the complementary distribution function @f$\bar{F}(x)@f$ with
 * parameter @f$n@f$. Simply returns `1 - cdf(n,x)`. It is not precise in the
 * upper tail.
 */
public static double barF (int n, double x) {
      double h = barFConnu(n, x);
      if (h >= 0.0)
         return h;

      h = 1.0 - cdf(n, x);
      if (h >= 0.0)
         return h;
      return 0.0;
   }


   protected static double inverseConnue (int n, double u) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u must be in [0,1]");
      if (u == 1.0)
         return 1.0;

      if (u == 0.0)
         return 0.5/n;

      if (n == 1)
         return (u + 1.0) / 2.0;

      final double NLNN = n*Math.log (n);
      final double LNU = Math.log(u) - Num.lnFactorial (n);
      if (LNU <= -NLNN){
         double t = 1.0/n*(LNU);
         return 0.5 * (Math.exp(t) + 1.0/n);
      }

      if (u >= 1.0 - 2.0 / Math.exp (NLNN))
         return 1.0 - Math.pow((1.0-u)/2.0, 1.0/n);

      return -1.0;
   }

/**
 * Computes the inverse @f$x = F^{-1}(u)@f$ of the *Kolmogorov–Smirnov*
 * distribution @f$F(x)@f$ with parameter @f$n@f$.
 */
public static double inverseF (int n, double u) {
      double Res = inverseConnue(n,u);
      if (Res != -1.0)
         return Res;
      Function f = new Function (n,u);
      return RootFinder.brentDekker (0.5/n, 1.0, f, 1e-10);
   }

   /**
    * Returns the parameter @f$n@f$ of this object.
    */
   public int getN() {
      return n;
   }

   /**
    * Sets the parameter @f$n@f$ of this object.
    */
   public void setN (int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      this.n = n;
      supportA = 0.5 / n;
      supportB = 1.0;
   }

   /**
    * Returns an array containing the parameter @f$n@f$ of this object.
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
