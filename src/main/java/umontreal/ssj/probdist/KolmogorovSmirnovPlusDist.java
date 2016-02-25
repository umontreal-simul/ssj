/*
 * Class:        KolmogorovSmirnovPlusDist
 * Description:  Kolmogorov-Smirnov+ distribution
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
import umontreal.ssj.util.*;
import umontreal.ssj.functions.MathFunction;

/**
 * Extends the class  @ref ContinuousDistribution for the
 * *Kolmogorov–Smirnov*+ distribution (see @cite tDAR60a, @cite tDUR73a,
 * @cite tBRO07a&thinsp;). Given a sample of @f$n@f$ independent uniforms
 * @f$U_i@f$ over @f$[0,1]@f$, the *Kolmogorov–Smirnov*+ statistic
 * @f$D_n^+@f$ and the *Kolmogorov–Smirnov*@f$-@f$ statistic @f$D_n^-@f$, are
 * defined by
 * @anchor REF_probdist_KolmogorovSmirnovPlusDist_eq_DNp
 * @anchor REF_probdist_KolmogorovSmirnovPlusDist_eq_DNm
 * @f{align}{
 *    D_n^+ 
 *    & 
 *   =
 *    \max_{1\le j\le n} \left(j/n - U_{(j)}\right), \tag{DNp} 
 *    \\ 
 *   D_n^- 
 *    & 
 *   =
 *    \max_{1\le j\le n} \left(U_{(j)} - (j-1)/n\right), \tag{DNm}
 * @f}
 * where the @f$U_{(j)}@f$ are the @f$U_i@f$ sorted in increasing order. Both
 * statistics follows the same distribution function, i.e. @f$F_n(x) =
 * P[D_n^+ \le x] = P[D_n^- \le x]@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_edf
 */
public class KolmogorovSmirnovPlusDist extends ContinuousDistribution {
   protected int n;

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
    * Constructs an *Kolmogorov–Smirnov*+ distribution for a sample of
    * size @f$n@f$.
    */
   public KolmogorovSmirnovPlusDist (int n) {
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

/**
 * Computes the density of the *Kolmogorov–Smirnov*+ distribution with
 * parameter @f$n@f$.
 */
public static double density (int n, double x) {
      if (n <= 0)
        throw new IllegalArgumentException (
                             "Calling kolmogorovSmirnovPlus with n < 1");
      if (x <= 0.0 || x >= 1.0)
         return 0.0;
      if (n == 1)
         return 1.0;
      final double EPS = 1.0 / 100.0;
      final double D1 = dclem(n, x, EPS);
      final double D2 = dclem(n, x, 2.0 * EPS);
      final double RES = D1 + (D1 - D2) / 3.0;
      if (RES < 0.0)
         return 0.0;
      return RES;
   }

   /**
    * Computes the *Kolmogorov–Smirnov*+ distribution function
    * @f$F_n(x)@f$ with parameter @f$n@f$.
    * The distribution function can be approximated via the following
    * expressions:
    * @anchor REF_probdist_KolmogorovSmirnovPlusDist_DistKS1
    * @anchor REF_probdist_KolmogorovSmirnovPlusDist_DistKS2
    * @anchor REF_probdist_KolmogorovSmirnovPlusDist_DistKS3
    * @f{align}{
    *    P[D_n^+ \le x] 
    *    & 
    *   =
    *    1 - x \sum_{i=0}^{\lfloor n(1-x)\rfloor} \binom{n}{i} \left(\frac{i}{n} + x \right)^{i-1} \left(1 - \frac{i}{n} - x \right)^{n-i} \tag{DistKS1} 
    *    \\  & 
    *   =
    *    x \sum_{j=0}^{\lfloor nx \rfloor} \binom{n}{j} \left(\frac{j}{n} - x \right)^j \left(1 - \frac{j}{n} + x \right)^{n-j-1} \tag{DistKS2} 
    *    \\  & 
    *   \approx
    *    1 - e^{-2 n x^2} \left[1 - \frac{2x}{3} \left( 1 - x\left(1 - \frac{2 n x^2}{3}\right) \right.\right. \nonumber
    *    \\  &   
    *    \left.\left. - \frac{2}{3n} \left(\frac{1}{5} - \frac{19 n x^2}{15} + \frac{2n^2 x^4}{3}\right) \right) + O(n^{-2}) \right]. \tag{DistKS3}
    * @f}
    * Formula (
    * {@link REF_probdist_KolmogorovSmirnovPlusDist_DistKS1
    * DistKS1} ) and (
    * {@link REF_probdist_KolmogorovSmirnovPlusDist_DistKS2
    * DistKS2} ) can be found in @cite tDUR73a&thinsp;, equations (2.1.12)
    * and (2.1.16), while (
    * {@link REF_probdist_KolmogorovSmirnovPlusDist_DistKS3
    * DistKS3} ) can be found in @cite tDAR60a&thinsp;. Formula (
    * {@link REF_probdist_KolmogorovSmirnovPlusDist_DistKS2
    * DistKS2} ) becomes numerically unstable as @f$nx@f$ increases. The
    * approximation (
    * {@link REF_probdist_KolmogorovSmirnovPlusDist_DistKS3
    * DistKS3} ) is simpler to compute and excellent when @f$nx@f$ is
    * large.
    * The relative error on @f$F_n(x) = P[D_n^+ \le x]@f$ is always less
    * than @f$10^{-5}@f$.
    */
   public static double cdf (int n, double x) {
      if (n <= 0)
        throw new IllegalArgumentException (
                             "Calling kolmogorovSmirnovPlus with n < 1");
      if (x <= 0.0)
         return 0.0;
      if ((x >= 1.0) || (n*x*x >= 25.0))
         return 1.0;
      if (n == 1)
         return x;

      final double NXPARAM = 6.5;      // frontier: alternating series
      final int NPARAM = 4000;         // frontier: non-alternating series
      double q;
      double Sum = 0.0;
      double term;
      double Njreal;
      double jreal;
      double LogCom = Math.log ((double)n);
      int j;
      int jmax;

      //--------------------------------------------------------------
      // the alternating series is stable and fast for n*x very small
      //--------------------------------------------------------------

      if (n*x <= NXPARAM) {
         final double EPSILON = Num.DBL_MIN;
         jmax = (int)(n*x);
         int Sign = -1;
         for (j = 1; j <= jmax; j++) {
            jreal = j;
            Njreal = n - j;
            q = jreal/n - x;
            // we must avoid log (0.0) for j = jmax and n*x near an integer
            if (-q > EPSILON) {
               term = LogCom + jreal*Math.log(-q) + (Njreal - 1.0)
                       *Math.log1p (-q);
               Sum += Sign*Math.exp (term);
            }
            Sign = -Sign;
            LogCom += Math.log (Njreal/(j + 1));
         }
         // add the term j = 0
         Sum += Math.exp ((n - 1)*Math.log1p (x));
         return Sum*x;
      }

      //-----------------------------------------------------------
      // For nx > NxParam, we use the other exact series for small
      // n, and the asymptotic form for n larger than NPARAM
      //-----------------------------------------------------------

      if (n <= NPARAM) {
         jmax = (int)(n*(1.0 - x));
         if (1.0 - x - (double) jmax/n <= 0.0)
            --jmax;
         for (j = 1; j <= jmax; j++) {
            jreal = j;
            Njreal = n - j;
            q = jreal/n + x;
            term = LogCom+(jreal - 1.0)*Math.log(q) + Njreal*Math.log1p (-q);
            Sum += Math.exp (term);
            LogCom += Math.log (Njreal/(jreal + 1.0));
         }
         Sum *= x;

         // add the term j = 0; avoid log (0.0)
         if (1.0 > x)
            Sum += Math.exp (n*Math.log1p (-x));
         return 1.0 - Sum;
      }

      //--------------------------
      // Use an asymptotic formula
      //--------------------------

      term = 2.0/3.0;
      q = x*x*n;
      Sum = 1.0 - Math.exp (-2.0*q)*(1.0 - term*x*(1.0 - x*(1.0 - term*q)
            - term/n*(0.2 - 19.0/15.0*q + term*q*q)));
      return Sum;
   }


   private static double KSPlusbarAsymp (int n, double x) {
      /* Compute the probability of the complementary KSPlus distribution
         using an asymptotic formula */
      double t = (6.0*n*x + 1);
      double z = t*t/(18.0*n);
      double v = 1.0 - (2.0*z*z - 4.0*z - 1.0)/(18.0*n);
      if (v <= 0.0)
         return 0.0;
      v = v*Math.exp(-z);
      if (v >= 1.0)
         return 1.0;
      return v;
   }


//-------------------------------------------------------------------------

   static double KSPlusbarUpper (int n, double x) {
      /* Compute the probability of the complementary KS+ distribution in
         the upper tail using Smirnov's stable formula */

      if (n > 200000)
         return KSPlusbarAsymp (n, x);

      int jmax = (int) (n* (1.0 - x));
      // Avoid log(0) for j = jmax and q ~ 1.0
      if ((1.0 - x - (double)jmax/n) <= 0.0)
         jmax--;

      int jdiv;
      if (n > 3000)
         jdiv = 2;
      else
         jdiv = 3;
      int j = jmax / jdiv + 1;

      double LogCom = Num.lnFactorial(n) - Num.lnFactorial(j)
                      - Num.lnFactorial(n-j);
      final double LOGJM = LogCom;

      final double EPSILON = 1.0E-12;
      double q;
      double term;
      double t;
      double Sum = 0.0;

      while (j <= jmax) {
         q = (double)j / n + x;
         term = LogCom + (j - 1)*Math.log(q) + (n-j)*Math.log1p(-q);
         t = Math.exp (term);
         Sum += t;
         LogCom += Math.log ((double)(n - j)/(j + 1));
         if (t <= Sum*EPSILON)
            break;
         j++;
      }

      j = jmax / jdiv;
      LogCom = LOGJM + Math.log ((double)(j+1)/(n - j));

      while (j > 0) {
         q = (double)j / n + x;
         term = LogCom+(j - 1)*Math.log (q)+ (n - j)*Math.log1p (-q);
         t = Math.exp (term);
         Sum += t;
         LogCom += Math.log ((double)j/(n - j + 1));
         if (t <= Sum*EPSILON)
            break;
         j--;
      }

      Sum *= x;
      // add the term j = 0
      Sum += Math.exp (n*Math.log1p (-x));
      return Sum;
   }

/**
 * Computes the complementary distribution function @f$\bar{F}_n(x)@f$ with
 * parameter @f$n@f$.
 */
public static double barF (int n, double x) {
      if (n <= 0)
        throw new IllegalArgumentException (
                             "Calling kolmogorovSmirnovPlus with n < 1");
      if (x <= 0.0)
         return 1.0;
      if ((x >= 1.0) || (n*x*x >= 365.0))
         return 0.0;
      if (n == 1)
         return 1.0 - x;

      final double NXPARAM = 6.5;    // frontier: alternating series
      final int NPARAM = 4000;       // frontier: non-alternating series
      final int NASYMP = 200000;     // frontier: asymptotic

      // the alternating series is stable and fast for n*x very small
      if (n*x <= NXPARAM)
         return 1.0 - cdf(n, x);

      if (n >= NASYMP)
         return KSPlusbarAsymp (n, x);

      if ((n <= NPARAM) || (n*x*x > 1.0))
         return KSPlusbarUpper(n,x);

      return KSPlusbarAsymp (n, x);
      // return (1.0 - 2.0*x/3.0)*Math.exp(-2.0*n*x*x);
   }

   /**
    * Computes the inverse @f$x = F^{-1}(u)@f$ of the distribution with
    * parameter @f$n@f$.
    */
   public static double inverseF (int n, double u) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u must be in [0,1]");
      if (u == 1.0)
         return 1.0;
      if (u == 0.0)
         return 0.0;
      Function f = new Function (n,u);
      return RootFinder.brentDekker (0.0, 1.0, f, 1e-8);
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
      supportA = 0.0;
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