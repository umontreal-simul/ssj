/*
 * Class:        AndersonDarlingDistQuick
 * Description:  Anderson-Darling distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
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
 * Extends the class  @ref AndersonDarlingDist for the *Anderson–Darling*
 * distribution (see @cite tAND52a, @cite tLEW61a, @cite tSTE86b&thinsp;).
 * This class implements a version faster and more precise in the tails than
 * class  @ref AndersonDarlingDist.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_edf
 */
public class AndersonDarlingDistQuick extends AndersonDarlingDist {

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

   //-------------------------------------------------------------------------

   private static double lower_Marsaglia (int n, double x) {
      // queue inférieure de l'algorithme de Marsaglia pour n = inf

      double p0 = (2.00012 + (.247105 - (.0649821 - (.0347962 -
                  (.011672 - .00168691 * x) * x) * x) * x) * x);
      p0 *= Math.exp (-1.2337141 / x) / Math.sqrt (x);
      return p0 >= 0 ? p0 : 0;
   }


   private static double diffcdf (int n, double x, double EPS) {
      return (cdf(n, x + EPS) - cdf(n, x - EPS)) / (2.0 * EPS);
   }

   /**
    * Constructs an *Anderson–Darling* distribution for a sample of size
    * @f$n@f$.
    */
   public AndersonDarlingDistQuick (int n) {
      super (n);
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

/**
 * Computes the density of the *Anderson–Darling* distribution with parameter
 * @f$n@f$.
 */
public static double density (int n, double x) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      if (n == 1)
         return density_N_1(x);

      if (x >= XBIG || x <= 0.0)
         return 0.0;
      final double EPS = 1.0 / 64.0;
      final double D1 = diffcdf(n, x, EPS);
      final double D2 = diffcdf(n, x, 2.0 * EPS);
      double res = D1 + (D1 - D2) / 3.0;
      return res >= 0. ? res : 0.;
   }

   // Tables for the approximation of the Anderson-Darling distribution
   private static double[] F2AD = new double[103];
   private static double[] CoAD = new double[103];

   static {
      F2AD[0] = 0.0;
      F2AD[1] = 1.7315E-10;
      F2AD[2] = 2.80781E-5;
      F2AD[3] = 1.40856E-3;
      F2AD[4] = 9.58772E-3;
      F2AD[5] = 2.960552E-2;
      F2AD[6] = 6.185146E-2;
      F2AD[7] = 1.0357152E-1;
      F2AD[8] = 1.5127241E-1;
      F2AD[9] = 2.0190317E-1;
      F2AD[10] = 2.5318023E-1;
      F2AD[11] = 3.0354278E-1;
      F2AD[12] = 3.5200015E-1;
      F2AD[13] = 3.9797537E-1;
      F2AD[14] = 4.4117692E-1;
      F2AD[15] = 4.8150305E-1;
      F2AD[16] = 5.1897375E-1;
      F2AD[17] = 5.5368396E-1;
      F2AD[18] = 5.8577199E-1;
      F2AD[19] = 6.1539864E-1;
      F2AD[20] = 6.4273362E-1;
      F2AD[21] = 6.6794694E-1;
      F2AD[22] = 6.9120359E-1;
      F2AD[23] = 7.126605E-1;
      F2AD[24] = 7.3246483E-1;
      F2AD[25] = 7.507533E-1;
      F2AD[26] = 7.6765207E-1;
      F2AD[27] = 7.8327703E-1;
      F2AD[28] = 7.9773426E-1;
      F2AD[29] = 8.1112067E-1;
      F2AD[30] = 8.2352466E-1;
      F2AD[31] = 8.3502676E-1;
      F2AD[32] = 8.4570037E-1;
      F2AD[33] = 8.5561231E-1;
      F2AD[34] = 8.6482346E-1;
      F2AD[35] = 8.7338931E-1;
      F2AD[36] = 8.8136046E-1;
      F2AD[37] = 8.8878306E-1;
      F2AD[38] = 8.9569925E-1;
      F2AD[39] = 9.0214757E-1;
      F2AD[40] = 9.081653E-1;
      F2AD[41] = 9.1378043E-1;
      F2AD[42] = 9.1902284E-1;
      F2AD[43] = 9.2392345E-1;
      F2AD[44] = 9.2850516E-1;
      F2AD[45] = 9.3279084E-1;
      F2AD[46] = 9.3680149E-1;
      F2AD[47] = 9.4055647E-1;
      F2AD[48] = 9.440736E-1;
      F2AD[49] = 9.4736933E-1;
      F2AD[50] = 9.5045883E-1;
      F2AD[51] = 9.5335611E-1;
      F2AD[52] = 9.5607414E-1;
      F2AD[53] = 9.586249E-1;
      F2AD[54] = 9.6101951E-1;
      F2AD[55] = 9.6326825E-1;
      F2AD[56] = 9.6538067E-1;
      F2AD[57] = 9.6736563E-1;
      F2AD[58] = 9.6923135E-1;
      F2AD[59] = 9.7098548E-1;
      F2AD[60] = 9.7263514E-1;
      F2AD[61] = 9.7418694E-1;
      F2AD[62] = 9.7564704E-1;
      F2AD[63] = 9.7702119E-1;
      F2AD[64] = 9.7831473E-1;
      F2AD[65] = 9.7953267E-1;
      F2AD[66] = 9.8067966E-1;
      F2AD[67] = 9.8176005E-1;
      F2AD[68] = 9.827779E-1;
      F2AD[69] = 9.8373702E-1;
      F2AD[70] = 9.8464096E-1;
      F2AD[71] = 9.8549304E-1;
      F2AD[72] = 9.8629637E-1;
      F2AD[73] = 9.8705386E-1;
      F2AD[74] = 9.8776824E-1;
      F2AD[75] = 9.8844206E-1;
      F2AD[76] = 9.8907773E-1;
      F2AD[77] = 9.8967747E-1;
      F2AD[78] = 9.9024341E-1;
      F2AD[79] = 9.9077752E-1;
      F2AD[80] = 9.9128164E-1;
      F2AD[81] = 9.9175753E-1;
      F2AD[82] = 9.9220682E-1;
      F2AD[83] = 9.9263105E-1;
      F2AD[84] = 9.9303165E-1;
      F2AD[85] = 9.9340998E-1;
      F2AD[86] = 9.9376733E-1;
      F2AD[87] = 9.9410488E-1;
      F2AD[88] = 9.9442377E-1;
      F2AD[89] = 9.9472506E-1;
      F2AD[90] = 9.9500974E-1;
      F2AD[91] = 9.9527876E-1;
      F2AD[92] = 9.95533E-1;
      F2AD[93] = 9.9577329E-1;
      F2AD[94] = 9.9600042E-1;
      F2AD[95] = 9.9621513E-1;
      F2AD[96] = 9.964181E-1;
      F2AD[97] = 0.99661;
      F2AD[98] = 9.9679145E-1;
      F2AD[99] = 9.9696303E-1;
      F2AD[100] = 9.9712528E-1;
      F2AD[101] = 9.9727872E-1;
      F2AD[102] = 9.9742384E-1;

      CoAD[0] = 0.0;
      CoAD[1] = 0.0;
      CoAD[2] = 0.0;
      CoAD[3] = 0.0;
      CoAD[4] = 0.0;
      CoAD[5] = -1.87E-3;
      CoAD[6] = 0.00898;
      CoAD[7] = 0.0209;
      CoAD[8] = 0.03087;
      CoAD[9] = 0.0377;
      CoAD[10] = 0.0414;
      CoAD[11] = 0.04386;
      CoAD[12] = 0.043;
      CoAD[13] = 0.0419;
      CoAD[14] = 0.0403;
      CoAD[15] = 0.038;
      CoAD[16] = 3.54804E-2;
      CoAD[17] = 0.032;
      CoAD[18] = 0.0293;
      CoAD[19] = 2.61949E-2;
      CoAD[20] = 0.0228;
      CoAD[21] = 0.0192;
      CoAD[22] = 1.59865E-2;
      CoAD[23] = 0.0129;
      CoAD[24] = 0.0107;
      CoAD[25] = 8.2464E-3;
      CoAD[26] = 0.00611;
      CoAD[27] = 0.00363;
      CoAD[28] = 1.32272E-3;
      CoAD[29] = -5.87E-4;
      CoAD[30] = -2.75E-3;
      CoAD[31] = -3.95248E-3;
      CoAD[32] = -5.34E-3;
      CoAD[33] = -6.892E-3;
      CoAD[34] = -8.10208E-3;
      CoAD[35] = -8.93E-3;
      CoAD[36] = -9.552E-3;
      CoAD[37] = -1.04605E-2;
      CoAD[38] = -0.0112;
      CoAD[39] = -1.175E-2;
      CoAD[40] = -1.20216E-2;
      CoAD[41] = -0.0124;
      CoAD[42] = -1.253E-2;
      CoAD[43] = -1.27076E-2;
      CoAD[44] = -0.0129;
      CoAD[45] = -1.267E-2;
      CoAD[46] = -1.22015E-2;
      CoAD[47] = -0.0122;
      CoAD[48] = -1.186E-2;
      CoAD[49] = -1.17218E-2;
      CoAD[50] = -0.0114;
      CoAD[51] = -1.113E-2;
      CoAD[52] = -1.08459E-2;
      CoAD[53] = -0.0104;
      CoAD[54] = -9.93E-3;
      CoAD[55] = -9.5252E-3;
      CoAD[56] = -9.24E-3;
      CoAD[57] = -9.16E-3;
      CoAD[58] = -8.8004E-3;
      CoAD[59] = -8.63E-3;
      CoAD[60] = -8.336E-3;
      CoAD[61] = -8.10512E-3;
      CoAD[62] = -7.94E-3;
      CoAD[63] = -7.71E-3;
      CoAD[64] = -7.55064E-3;
      CoAD[65] = -7.25E-3;
      CoAD[66] = -7.11E-3;
      CoAD[67] = -6.834E-3;
      CoAD[68] = -0.0065;
      CoAD[69] = -6.28E-3;
      CoAD[70] = -6.11008E-3;
      CoAD[71] = -5.86E-3;
      CoAD[72] = -5.673E-3;
      CoAD[73] = -5.35008E-3;
      CoAD[74] = -5.11E-3;
      CoAD[75] = -4.786E-3;
      CoAD[76] = -4.59144E-3;
      CoAD[77] = -4.38E-3;
      CoAD[78] = -4.15E-3;
      CoAD[79] = -4.07696E-3;
      CoAD[80] = -3.93E-3;
      CoAD[81] = -3.83E-3;
      CoAD[82] = -3.74656E-3;
      CoAD[83] = -3.49E-3;
      CoAD[84] = -3.33E-3;
      CoAD[85] = -3.20064E-3;
      CoAD[86] = -3.09E-3;
      CoAD[87] = -2.93E-3;
      CoAD[88] = -2.78136E-3;
      CoAD[89] = -2.72E-3;
      CoAD[90] = -2.66E-3;
      CoAD[91] = -2.56208E-3;
      CoAD[92] = -2.43E-3;
      CoAD[93] = -2.28E-3;
      CoAD[94] = -2.13536E-3;
      CoAD[95] = -2.083E-3;
      CoAD[96] = -1.94E-3;
      CoAD[97] = -1.82E-3;
      CoAD[98] = -1.77E-3;
      CoAD[99] = -1.72E-3;
      CoAD[100] = -1.71104E-3;
      CoAD[101] = -1.741E-3;
      CoAD[102] = -0.0016;
   }

/**
 * Computes the *Anderson–Darling* distribution function @f$F_n(x)@f$ at
 * @f$x@f$ for sample size @f$n@f$. For @f$0.2 < x < 5@f$, the asymptotic
 * distribution @f$F_{\infty}(x) = \lim_{n\to\infty} F_n(x)@f$ was first
 * computed by numerical integration; then a linear correction @f$O(1/n)@f$
 * obtained by simulation was added. For @f$5 < x@f$, the Grace-Wood
 * empirical approximation @cite tGRA12a&thinsp; is used. For @f$x < 0.2@f$,
 * the Marsaglias’ approximation @cite tMAR04a&thinsp; for @f$n=\infty@f$ is
 * used.
 *
 * For @f$n>6@f$, the method gives at least 3 decimal digits of precision
 * except for small @f$x@f$; for @f$n \le6@f$, it gives at least 2 decimal
 * digits of precision except for small @f$x@f$. For @f$n=1@f$, the exact
 * formula @f$F_1(x) = \sqrt{1 - 4e^{-x-1}}@f$, for @f$x\ge\ln(4) - 1@f$,
 * is used.
 */
public static double cdf (int n, double x) {
      if (n <= 0)
         throw new IllegalArgumentException ("   n <= 0");
      if (1 == n)
         return cdf_N_1 (x);
      if (x <= 0.0)
         return 0.0;
      if (x >= XBIG)
         return 1.0;
      if (x <= 0.2)
         return lower_Marsaglia (n, x);
      return 1.0 - barF (n, x);
   }

   /**
    * Computes the complementary distribution function @f$\bar{F}_n(x)@f$
    * with parameter @f$n@f$.
    */
   public static double barF (int n, double x) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      if (n == 1)
         return barF_N_1 (x);

      if (x <= 0.0)
         return 1.0;
      if (x >= XBIGM)
         return 0.0;

      double q;
      if (x > 5.) {
         // Grace-Wood approx. in upper tail
         double nd = n;
         q = (0.23945*Math.pow(nd, -0.9379) - 0.1201*Math.pow(nd, -0.96) -
             1.0002816)*x - 1.437*Math.pow(nd, -0.9379) +
             1.441*Math.pow(nd, -0.96) - 0.0633101;
         return Math.pow(x, -0.48897) * Math.exp(q);
      }

      if (x <= 0.2)
         return 1.0 - cdf (n, x);

      final double H = 0.05;  // the step of the interpolation table
      final int i = (int) (1 +  x / H);
      double res;
      q = x/H - i;

      // Newton backwards quadratic interpolation
      res = (F2AD[i - 2] - 2.0*F2AD[i - 1] + F2AD[i])*q*(q + 1.0)/2.0
         + (F2AD[i] - F2AD[i - 1])*q + F2AD[i];

      // Empirical correction in 1/n
      res += (CoAD[i]*(q + 1.0) - CoAD[i - 1]*q)/n;

      res = 1.0 - res;
      if (res >= 1.0)
         return 1.0;
      if (res <= 0.0)
         return 0.0;
   return res;
   }

   /**
    * Computes the inverse @f$x = F_n^{-1}(u)@f$ of the *Anderson–Darling*
    * distribution with parameter @f$n@f$.
    */
   public static double inverseF (int n, double u) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u must be in [0,1]");
      if (n == 1)
         return inverse_N_1 (u);
      if (u == 1.0)
         return Double.POSITIVE_INFINITY;
      if (u == 0.0)
         return 0.0;
      Function f = new Function (n,u);
      return RootFinder.brentDekker (0.0, 50.0, f, 1.0e-5);
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : n = " + n;
   }

}