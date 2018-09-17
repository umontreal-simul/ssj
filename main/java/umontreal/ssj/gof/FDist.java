/*
 * Class:        FDist
 * Description:  Empirical distributions
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
package umontreal.ssj.gof;

/**
 * This class provides methods to compute (or approximate) the distribution
 * functions of special types of goodness-of-fit test statistics.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class FDist {
   private FDist() {}

/**
 * Similar to  @ref umontreal.ssj.probdist.KolmogorovSmirnovPlusDist but for
 * the case where the distribution function @f$F@f$ has a jump of size
 * @f$a@f$ at a given point @f$x_0@f$, is zero at the left of @f$x_0@f$, and
 * is continuous at the right of @f$x_0@f$.
 *  The Kolmogorov-Smirnov statistic is defined in that case as
 * @anchor REF_gof_FDist_eq_KSPlusJumpOne
 * @f[
 *   D_N^+(a) = \sup_{a\le u\le1} \left(\hat{F}_N (F^{-1}(u)) - u\right) = \max_{\Rule{0.0pt}{7.0pt}{0.0pt} \lfloor1+aN \le j \le N} \left(j/N - F (V_{(j)})\right). \tag{KSPlusJumpOne}
 * @f]
 * where @f$V_{(1)},…,V_{(N)}@f$ are the observations sorted by increasing
 * order. The method returns an approximation of @f$P[D_N^+(a) \le x]@f$
 * computed via
 * @anchor REF_gof_FDist_DistKSJ1
 * @anchor REF_gof_FDist_DistKSJ2
 * @f{align}{
 *    P[D_N^+(a) \le x] 
 *    & 
 *   =
 *    1 - x \sum_{i=0}^{\lfloor N (1-a-x)\rfloor} \binom{N}{i} \left(\frac{i}{N} + x \right)^{i-1} \left(1 - \frac{i}{N} - x \right)^{N-i}. \tag{DistKSJ1} 
 *    \\  & 
 *   =
 *    x \sum_{j=0}^{\lfloor N (a+x) \rfloor} \binom{N}{j} \left(\frac{j}{N} - x \right)^j \left(1 - \frac{j}{N} + x \right)^{N-j-1}. \tag{DistKSJ2}
 * @f}
 * The current implementation uses formula (
 * {@link REF_gof_FDist_DistKSJ2 DistKSJ2} ) when @f$N (x+a) <
 * 6.5@f$ and @f$x+a < 0.5@f$, and uses (
 * {@link REF_gof_FDist_DistKSJ1 DistKSJ1} ) when @f$Nx
 * \ge6.5@f$ or @f$x+a \ge0.5@f$.
 *  Restriction: @f$0 < a < 1@f$.
 *  @param N            sample size
 *  @param a            size of the jump
 *  @param x            positive or negative Kolmogorov-Smirnov statistic
 *  @return the distribution function of the statistic evaluated at `x`
 */
public static double kolmogorovSmirnovPlusJumpOne (int N, double a,
                                                      double x) {
      final double EPSILONLR = 1.E-15;
      final double EPSILON = 1.0E-290;
      final double NXAPARAM = 6.5;   // frontier: alternating series
      double LogCom;
      double q, p1, q1;
      double Sum = 0.0;
      double term;
      double Njreal;
      double jreal;
      int Sign;
      int j;
      int jmax;

      if (N < 1)
        throw new IllegalArgumentException (
                             "Calling kolmogorovSmirnovPlusJumpOne "+
                             "with N < 1");
      if (a >= 1.0 || a <= 0.0)
        throw new IllegalArgumentException (
                             "Calling kolmogorovSmirnovPlusJumpOne "+
                             "with a outside (0, 1)");
      if (x <= 0.0)
         return 0.0;
      if (x + a >= 1.0)
         return 1.0;
      LogCom = Math.log ((double)N);

      //--------------------------------------------------------------------
      // the alternating series is stable and fast for N*(x + a) very small
      //--------------------------------------------------------------------
      if (N*(x + a) < NXAPARAM && a + x < 0.5) {
         jmax = (int)(N*(x + a));
         for (j = 1; j <= jmax; j++) {
            jreal = j;
            Njreal = N - j;
            q = jreal/N - x;
            if ((q < 0.0 && ((j & 1) != 0)) ||
               ((q > 1.0) && (((N - j - 1) & 1) != 0)))
               Sign = -1;
            else
               Sign = 1;

            // we must avoid log (0.0)
            q1 = Math.abs (q);
            p1 = Math.abs (1.0 - q);
            if (q1 > EPSILON && p1 > EPSILON) {
               term = LogCom + jreal*Math.log (q1) +
                     (Njreal - 1.0)*Math.log (p1);
               Sum += Sign*Math.exp (term);
            }
            LogCom += Math.log (Njreal/(jreal + 1.0));
         }
         // add the term j = 0
         Sum += Math.exp ((N - 1)*Math.log (1.0 + x));
         return Sum*x;
      }

      //---------------------------------------------
      // For N (x + a) >= NxaParam or (a + x) > 0.5,
      // use the non-alternating series.
      //---------------------------------------------

      // EpsilonLR because the distribution has a jump
      jmax = (int)(N*(1.0 - a - x - EPSILONLR));
      for (j = 1; j <= jmax; j++) {
         jreal = j;
         Njreal = N - jreal;
         q = jreal/N + x;
         if (1.0 - q > EPSILON) {
            term = LogCom + (jreal - 1.0)*Math.log (q) + Njreal*Math.log (1.0 - q);
            Sum += Math.exp (term);
         }
         LogCom += Math.log (Njreal/(jreal + 1.0));
      }
      Sum *= x;

      // add the term j = 0
      if (1.0 - x > EPSILON)
         Sum += Math.exp (N*Math.log (1.0 - x));
      return 1.0 - Sum;
   }

   /**
    * Returns @f$F (m)@f$, the distribution function of the scan statistic
    * with parameters @f$N@f$ and @f$d@f$, evaluated at @f$m@f$. For a
    * description of this statistic and its distribution, see
    * FBar.scan(int,double,int), which computes its complementary
    * distribution @f$\bar{F} (m) = 1 - F (m-1)@f$.
    *  @param N            sample size (@f$\ge2@f$)
    *  @param d            length of the test interval (@f$\in(0,1)@f$)
    *  @param m            scan statistic
    *  @return the distribution function of the statistic evaluated at `m`
    */
   public static double scan (int N, double d, int m) {
      return 1.0 - FBar.scan (N, d, m);
   }
}