/*
 * Class:        CramerVonMisesDist
 * Description:  Cramér-von Mises distribution
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
 * Extends the class  @ref ContinuousDistribution for the Cramér-von Mises
 * distribution (see @cite tDUR73a, @cite tSTE70a, @cite tSTE86b&thinsp;).
 * Given a sample of @f$n@f$ independent uniforms @f$U_i@f$ over @f$[0,1]@f$,
 * the Cramér-von Mises statistic @f$W_n^2@f$ is defined by
 * @anchor REF_probdist_CramerVonMisesDist_eq_CraMis
 * @f[
 *   W_n^2 = \frac{1}{12n} + \sum_{j=1}^n \left(U_{(j)} - \frac{(j-0.5)}{n}\right)^2, \tag{CraMis}
 * @f]
 * where the @f$U_{(j)}@f$ are the @f$U_i@f$ sorted in increasing order. The
 * distribution function (the cumulative probabilities) is defined as
 * @f$F_n(x) = P[W_n^2 \le x]@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_edf
 */
public class CramerVonMisesDist extends ContinuousDistribution {
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
    * Constructs a <em>Cramér-von Mises</em> distribution for a sample of
    * size @f$n@f$.
    */
   public CramerVonMisesDist (int n) {
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

   public double getMean() {
      return CramerVonMisesDist.getMean (n);
   }

   public double getVariance() {
      return CramerVonMisesDist.getVariance (n);
   }

   public double getStandardDeviation() {
      return CramerVonMisesDist.getStandardDeviation (n);
   }

/**
 * Computes the density function for a <em>Cramér-von Mises</em> distribution
 * with parameter @f$n@f$.
 */
public static double density (int n, double x) {
      if (n <= 0)
        throw new IllegalArgumentException ("n <= 0");

      if (x <= 1.0/(12.0*n) || x >= n/3.0)
         return 0.0;

      if (n == 1)
         return 1.0 / Math.sqrt (x - 1.0/12.0);

      if (x <= 0.002 || x > 3.95)
         return 0.0;

      throw new UnsupportedOperationException("density not implemented.");
   }

   /**
    * Computes the Cramér-von Mises distribution function with parameter
    * @f$n@f$. Returns an approximation of @f$P[W_n^2 \le x]@f$, where
    * @f$W_n^2@f$ is the Cramér von Mises statistic (see @cite tSTE70a,
    * @cite tSTE86b, @cite tAND52a, @cite tKNO74a&thinsp;). The
    * approximation is based on the distribution function of @f$W^2 =
    * \lim_{n\to\infty} W_n^2@f$, which has the following series
    * expansion derived by Anderson and Darling @cite tAND52a&thinsp;:
    * @f[
    *   \qquad P(W^2 \le x)  =  \frac{1}{\pi\sqrt{x}} \sum_{j=0}^{\infty}(-1)^j \binom{-1/2}{j} \sqrt{4j+1}\;\; {exp}\left\{-\frac{(4j+1)^2}{16 x}\right\} K_{1/4}\left(\frac{(4j+1)^2}{16 x}\right),
    * @f]
    * where @f$K_{\nu}@f$ is the modified Bessel function of the second
    * kind. To correct for the deviation between @f$P(W_n^2\le x)@f$ and
    * @f$P(W^2\le x)@f$, we add a correction in @f$1/n@f$, obtained
    * empirically by simulation. For @f$n = 10@f$, 20, 40, the error is
    * less than 0.002, 0.001, and 0.0005, respectively, while for @f$n
    * \ge100@f$ it is less than 0.0005. For @f$n \to\infty@f$, we
    * estimate that the method returns at least 6 decimal digits of
    * precision. For @f$n = 1@f$, the method uses the exact distribution:
    * @f$P(W_1^2 \le x) = 2 \sqrt{x - 1/12}@f$ for @f$1/12 \le x
    * \le1/3@f$.
    */
   public static double cdf (int n, double x) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      if (n == 1) {
         if (x <= 1.0/12.0)
            return 0.0;
         if (x >= 1.0/3.0)
            return 1.0;
         return 2.0*Math.sqrt (x - 1.0/12.0);
      }

      if (x <= 1.0/(12.0*n))
         return 0.0;

      if (x <= (n + 3.0)/(12.0*n*n)) {
         double t = Num.lnFactorial(n) - Num.lnGamma (1.0 + 0.5*n) +
            0.5*n*Math.log (Math.PI*(x - 1.0/(12.0*n)));
         return Math.exp(t);
      }

      if (x <= 0.002)
         return 0.0;
      if (x > 3.95 || x >= n/3.0)
         return 1.0;

      final double EPSILON = Num.DBL_EPSILON;
      final int JMAX = 20;
      int j = 0;
      double Cor, Res, arg;
      double termX, termS, termJ;

      termX = 0.0625/x;            // 1 / (16x)
      Res = 0.0;

      final double A[] = {
         1.0,
         1.11803398875,
         1.125,
         1.12673477358,
         1.1274116945,
         1.12774323743,
         1.1279296875,
         1.12804477649,
         1.12812074678,
         1.12817350091
      };

      do {
         termJ = 4*j + 1;
         arg = termJ*termJ*termX;
         termS = A[j]*Math.exp (-arg)*Num.besselK025 (arg);
         Res += termS;
         ++j;
      } while (!(termS < EPSILON || j > JMAX));

      if (j > JMAX)
         System.err.println ("cramerVonMises: iterations have not converged");
      Res /= Math.PI*Math.sqrt (x);

      // Empirical correction in 1/n
      if (x < 0.0092)
         Cor = 0.0;
      else if (x < 0.03)
         Cor = -0.0121763 + x*(2.56672 - 132.571*x);
      else if (x < 0.06)
         Cor = 0.108688 + x*(-7.14677 + 58.0662*x);
      else if (x < 0.19)
         Cor = -0.0539444 + x*(-2.22024 + x*(25.0407 - 64.9233*x));
      else if (x < 0.5)
         Cor = -0.251455 + x*(2.46087 + x*(-8.92836 + x*(14.0988 -
                  x*(5.5204 + 4.61784*x))));
      else if (x <= 1.1)
         Cor = 0.0782122 + x*(-0.519924 + x*(1.75148 +
               x*(-2.72035 + x*(1.94487 - 0.524911*x))));
      else
         Cor = Math.exp (-0.244889 - 4.26506*x);

      Res += Cor/n;
      // This empirical correction is not very precise, so ...
      if (Res <= 1.0)
         return Res;
      else
         return 1.0;
   }

   /**
    * Computes the complementary distribution function @f$\bar{F}_n(x)@f$
    * with parameter @f$n@f$.
    */
   public static double barF (int n, double x) {
      return 1.0 - cdf(n,x);
   }

   /**
    * Computes @f$x = F_n^{-1}(u)@f$, where @f$F_n@f$ is the
    * <em>Cramér-von Mises</em> distribution with parameter @f$n@f$.
    */
   public static double inverseF (int n, double u) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u must be in [0,1]");

      if (u >= 1.0)
         return n/3.0;
      if (u <= 0.0)
         return 1.0/(12.0*n);

      if (n == 1)
         return 1.0/12.0 + 0.25*u*u;

      Function f = new Function (n,u);
      return RootFinder.brentDekker (0.0, 10.0, f, 1e-6);
   }

   /**
    * Returns the mean of the distribution with parameter @f$n@f$.
    *  @return the mean
    */
   public static double getMean (int n) {
      return 1.0 / 6.0;
   }

   /**
    * Returns the variance of the distribution with parameter @f$n@f$.
    *  @return variance
    */
   public static double getVariance (int n) {
      return (4.0 * n - 3.0) / (180.0 * n ); 
   }

   /**
    * Returns the standard deviation of the distribution with parameter
    * @f$n@f$.
    *  @return the standard deviation
    */
   public static double getStandardDeviation (int n) {
     return Math.sqrt(getVariance(n));
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
      supportA = 1.0/(12.0*n);
      supportB = n/3.0;
   }

   /**
    * Return an array containing the parameter @f$n@f$ of this object.
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