/*
 * Class:        WatsonUDist
 * Description:  Watson U  distribution 
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
 * Extends the class  @ref ContinuousDistribution for the *Watson U*
 * distribution (see @cite tDUR73a, @cite tSTE70a, @cite tSTE86b&thinsp;).
 * Given a sample of @f$n@f$ independent uniforms @f$u_i@f$ over @f$[0,1]@f$,
 * the *Watson* statistic @f$U_n^2@f$ is defined by
 * @anchor REF_probdist_WatsonUDist_eq_WatsonU
 * @f{align*}{
 *    W_n^2 
 *    & 
 *   =
 *    \frac{1}{12n} + \sum_{j=1}^n \left\{u_{(j)} - \frac{(j- 1/2)}{n}\right\}^2, 
 *    \\ 
 *   U_n^2 
 *    & 
 *   =
 *    W_n^2 - n\left(\bar{u}_n - 1/2\right)^2. \tag{WatsonU}
 * @f}
 * where the @f$u_{(j)}@f$ are the @f$u_i@f$ sorted in increasing order, and
 * @f$\bar{u}_n@f$ is the average of the observations @f$u_i@f$. The
 * distribution function (the cumulative probabilities) is defined as
 * @f$F_n(x) = P[U_n^2 \le x]@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_edf
 */
public class WatsonUDist extends ContinuousDistribution {
   private static final double XSEPARE = 0.15;
   private static final double PI = Math.PI;
   private static final int JMAX = 10;
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

   private static double cdfn (int n, double x) {
      // The 1/n correction for the cdf, for x < XSEPARE
      double terme;
      double v = Math.exp (-0.125/x);
      double somme = 0;
      int j = 0;

      do {
         double a = (2*j + 1)*(2*j + 1);
         terme = Math.pow (v, (double)(2*j + 1)*(2*j + 1));
         double der = terme*(a - 4.0*x)/(8.0*x*x);
         somme += (5.0*x - 1.0/12.0) * der / 12.0;
         der = terme* (a*a - 24.0*a*x + 48.0*x*x)/ (64.0*x*x*x*x);
         somme += x*x*der/6.0;
         ++j;
      } while (!(terme <= Math.abs(somme) * Num.DBL_EPSILON || j > JMAX));
      if (j > JMAX)
         System.err.println (x + ": watsonU:  somme 1/n has not converged");

      v = -2.0*somme/(n*Math.sqrt (2.0*PI*x));
      return v;
   }

   /**
    * Constructs a *Watson U* distribution for a sample of size @f$n@f$.
    */
   public WatsonUDist (int n) {
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
      return WatsonUDist.getMean (n);
   }

   public double getVariance() {
      return WatsonUDist.getVariance (n);
   }

   public double getStandardDeviation() {
      return WatsonUDist.getStandardDeviation (n);
   }

/**
 * Computes the density of the *Watson U* distribution with parameter
 * @f$n@f$.
 */
public static double density (int n, double x) {
      if (n < 2)
         throw new IllegalArgumentException ("n < 2");

      if (x <= 1.0/(12.0*n) || x >= n/12.0 || x >= XBIG)
         return 0.0;

      final double EPS = 1.0 / 100.0;
      return (cdf(n, x + EPS) - cdf(n, x - EPS)) / (2.0 * EPS);

/*
// This is the asymptotic density for n -> infinity
      int j;
      double signe;
      double v;
      double terme;
      double somme;

      if (x > XSEPARE) {
         // this series converges rapidly for x > 0.15
         v = Math.exp (-(x*2.0*PI*PI));
         signe = 1.0;
         somme = 0.0;
         j = 1;
         do {
            terme = j*j*Math.pow (v, (double)j*j);
            somme += signe*terme;
            signe = -signe;
            ++j;
         } while (!(terme < Num.DBL_EPSILON || j > JMAX));
         if (j > JMAX)
            System.err.println ("watsonU:  sum1 has not converged");
         return 4.0*PI*PI*somme;
      }

      // this series converges rapidly for x <= 0.15
      v = Math.exp (-0.125/x);
      somme = v;
      double somme2 = v;
      j = 2;
      do {
         terme = Math.pow (v, (double)(2*j - 1)*(2*j - 1));
         somme += terme;
         somme2 += terme * (2*j - 1)*(2*j - 1);
         ++j;
      } while (!(terme <= somme2 * Num.DBL_EPSILON || j > JMAX));
      if (j > JMAX)
         System.err.println ("watsonU:  sum2 has not converged");

      final double RACINE = Math.sqrt (2.0*PI*x);
      return -somme/(x*RACINE) + 2*somme2 * 0.125/ ((x*x) * RACINE);
*/
   }

   /**
    * Computes the Watson @f$U@f$ distribution function, i.e. returns
    * @f$P[U_n^2 \le x]@f$, where @f$U_n^2@f$ is the Watson statistic
    * defined in (
    * {@link REF_probdist_WatsonUDist_eq_WatsonU WatsonU} ).
    * We use the asymptotic distribution for @f$n \to\infty@f$, plus a
    * correction in @f$O(1/n)@f$, as given in @cite tCSO96a&thinsp;.
    */
   public static double cdf (int n, double x) {
      if (n < 2)
         throw new IllegalArgumentException ("n < 2");

      if (x <= 1.0/(12.0*n))
         return 0.0;
      if (x > 3.95 || x >= n/12.0)
         return 1.0;

      if (2 == n) {
         if (x <= 1.0/24.0)
            return 0.0;
         if (x >= 1.0/6.0)
            return 1.0;
         return 2.0*Math.sqrt(2.0*x - 1.0/12.0);
      }

      if (x > XSEPARE)
         return 1.0 - barF (n, x);

      // this series converges rapidly for x <= 0.15
      double terme;
      double v = Math.exp (-0.125/x);
      double somme = v;
      int j = 2;

      do {
         terme = Math.pow (v, (double)(2*j - 1)*(2*j - 1));
         somme += terme;
         ++j;
      } while (!(terme <= somme * Num.DBL_EPSILON || j > JMAX));
      if (j > JMAX)
         System.err.println (x + ": watsonU:  sum2 has not converged");

      v = 2.0*somme/Math.sqrt (2.0*PI*x);
      v += cdfn(n, x);
      if (v >= 1.0)
         return 1.0;
      if (v <= 0.0)
         return 0.0;
       return v;
   }

   /**
    * Computes the complementary distribution function @f$\bar{F}_n(x)@f$,
    * where @f$F_n@f$ is the *Watson* @f$U@f$ distribution with parameter
    * @f$n@f$.
    */
   public static double barF (int n, double x) {
      if (n < 2)
         throw new IllegalArgumentException ("n < 2");

      if (x <= 1.0/(12.0*n))
         return 1.0;
      if (x >= XBIG || x >= n/12.0)
         return 0.0;

      if (2 == n)
         return 1.0 - 2.0*Math.sqrt(2.0*x - 1.0/12.0);

      if (x > XSEPARE) {
         // this series converges rapidly for x > 0.15
         double terme, ter;
         double v = Math.exp (-2.0*PI*PI*x);
         double signe = 1.0;
         double somme = 0.0;
         double son = 0.0;
         int j = 1;

         do {
            terme = Math.pow (v, (double)j*j);
            somme += signe*terme;
            double h = 2*j*PI*x;
            ter = (5.0*x - h*h - 1.0/12.0)*j*j;
            son += signe*terme*ter;
            signe = -signe;
            ++j;
         } while (!(terme < Num.DBL_EPSILON || j > JMAX));
         if (j > JMAX)
            System.err.println (x + ": watsonU:  sum1 has not converged");
         v = 2.0*somme + PI*PI*son/(3.0*n);
         if (v <= 0.0) 
            return 0.0;
         if (v >= 1.0)
            return 1.0;
         return v;
      }

      return (1.0 - cdf(n, x));
   }

   /**
    * Computes @f$x = F_n^{-1}(u)@f$, where @f$F_n@f$ is the *Watson*
    * @f$U@f$ distribution with parameter @f$n@f$.
    */
   public static double inverseF (int n, double u) {
      if (n < 2)
         throw new IllegalArgumentException ("n < 2");
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u must be in [0,1]");
      if (u >= 1.0)
         return n/12.0;
      if (u <= 0.0)
         return 1.0/(12.0*n);

      if (2 == n)
         return 1.0/24.0 + u*u/8.0;

      Function f = new Function (n,u);
      return RootFinder.brentDekker (0.0, 2.0, f, 1e-7);
   }

   /**
    * Returns the mean of the *Watson* @f$U@f$ distribution with parameter
    * @f$n@f$.
    *  @return Returns the mean
    */
   public static double getMean (int n) {
      return 1.0/12.0;
   }

   /**
    * Returns the variance of the *Watson* @f$U@f$ distribution with
    * parameter @f$n@f$.
    *  @return the variance
    */
   public static double getVariance (int n) {
      return (n - 1)/(360.0*n);
   }

   /**
    * Returns the standard deviation of the *Watson* @f$U@f$ distribution
    * with parameter @f$n@f$.
    *  @return the standard deviation
    */
   public static double getStandardDeviation (int n) {
      return Math.sqrt (WatsonUDist.getVariance (n));
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
      if (n < 2)
         throw new IllegalArgumentException ("n < 2");
      this.n = n;
      supportA = 1.0/(12.0*n);
      supportB = n/12.0;
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