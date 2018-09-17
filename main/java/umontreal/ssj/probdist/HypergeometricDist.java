/*
 * Class:        HypergeometricDist
 * Description:  hypergeometric distribution
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
import umontreal.ssj.util.Num;

/**
 * Extends the class  @ref DiscreteDistributionInt for the *hypergeometric*
 * distribution @cite rGEN98a&thinsp; (page 101) with @f$k@f$ elements chosen
 * among @f$l@f$, @f$m@f$ being of one type, and @f$l-m@f$ of the other. The
 * parameters @f$m@f$, @f$k@f$ and @f$l@f$ are positive integers where
 * @f$1\le m\le l@f$ and @f$1\le k\le l@f$. Its mass function is given by
 * @anchor REF_probdist_HypergeometricDist_eq_fheperg
 * @f[
 *   p(x) = \frac{ \binom{m}{x} \binom{l - m}{k-x}}{\binom{l}{k}} \qquad\mbox{for } \max(0,k-l+m)\le x\le\min(k, m). \tag{fheperg}
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_discrete
 */
public class HypergeometricDist extends DiscreteDistributionInt {

   private int m;
   private int l;
   private int k;
   private double p0;

   /**
    * @name Constant
    * @{
    */

   /**
    * If the number of integers in the interval @f$[\max(0,k-l+m),
    * \min(k,m)]@f$ is larger than this constant, the tables will *not*
    * be precomputed by the constructor.
    */
   public static double MAXN = 100000;

   /**
    * @}
    */

   /**
    * Constructs an hypergeometric distribution with parameters @f$m@f$,
    * @f$l@f$ and @f$k@f$.
    */
   public HypergeometricDist (int m, int l, int k) {
      setParams (m, l, k);
   }


   public double prob (int x) {
      if (x < supportA || x > supportB)
         return 0.0;
      if (pdf == null || x < xmin || x > xmax)
         return prob (m, l, k, x);
      return pdf[x - xmin];
   }

   public double cdf (int x) {
      if (x < supportA)
         return 0.0;
      if (x >= supportB)
         return 1.0;
      if (cdf != null) {
         if (x >= xmax)
            return 1.0;
         if (x < xmin)
            return cdf (m, l, k, x);
         if (x <= xmed)
            return cdf[x - xmin];
         else
            // We keep the complementary distribution in the upper part of cdf
            return 1.0 - cdf[x + 1 - xmin];
      }
      else
         return cdf (m, l, k, x);
   }

   public double barF (int x) {
      if (x <= supportA)
         return 1.0;
      if (x > supportB)
         return 0.0;
      if (cdf != null) {
         if (x > xmax)
            return barF (m, l, k, x);
         if (x <= xmin)
            return 1.0;
         if (x > xmed)
           // We keep the complementary distribution in the upper part of cdf
           return cdf[x - xmin];
         else
            return 1.0 - cdf[x - 1 - xmin];
      }
      else
         return barF (m, l, k, x);
   }

   public int inverseFInt (double u) {
      if (u < 0 || u > 1)
         throw new IllegalArgumentException ("u is not in [0,1]");
      if (u <= 0.0)
         return Math.max (0, k - l + m);
      if (u >= 1.0)
         return Math.min (k, m);
      double p = p0;
     // Empirical correction, the original algorithm sets x=0.
     int x = Math.max (0, k - l + m);
     if (u <= p) return x;
     do {
        u = u - p;
        p = p*(m - x)*(k - x)/((x+1)*(l - m - k + 1.0 + x));
        x++;
     } while (u > p);
     return x;
   }

   public double getMean() {
      return HypergeometricDist.getMean (m, l, k);
   }

   public double getVariance() {
      return HypergeometricDist.getVariance (m, l, k);
   }

   public double getStandardDeviation() {
      return HypergeometricDist.getStandardDeviation (m, l, k);
   }

/**
 * Computes the hypergeometric probability @f$p(x)@f$ given by (
 * {@link REF_probdist_HypergeometricDist_eq_fheperg fheperg}
 * ).
 */
public static double prob (int m, int l, int k, int x) {
      final int SLIM = 70;
      final double MAXEXP = (Num.DBL_MAX_EXP - 1)*Num.LN2;// To avoid overflow
      if (l <= 0)
         throw new IllegalArgumentException ("l must be greater than 0");
      if (m <= 0 || m > l)
         throw new IllegalArgumentException ("m is invalid: 1 <= m < l");
      if (k <= 0 || k > l)
         throw new IllegalArgumentException ("k is invalid: 1 <= k < l");
      if (x < Math.max (0, k - l + m) || x > Math.min (k, m))
         return 0;

      if (l <= SLIM)
         return Num.combination (m, x)
            * Num.combination (l - m, k - x)/Num.combination (l, k);
      else {
         double res =
             Num.lnFactorial (m) + Num.lnFactorial (l-m) - Num.lnFactorial (l)
           - Num.lnFactorial (x) - Num.lnFactorial (k - x) + Num.lnFactorial (k)
           - Num.lnFactorial (m - x) - Num.lnFactorial (l - m - k + x)
           + Num.lnFactorial (l - k);
         if (res >= MAXEXP)
            throw new IllegalArgumentException ("term overflow");
         return Math.exp (res);
      }
   }

   /**
    * Computes the distribution function @f$F(x)@f$.
    */
   public static double cdf (int m, int l, int k, int x) {
      if (l <= 0)
         throw new IllegalArgumentException ("l must be greater than 0");
      if (m <= 0 || m > l)
         throw new IllegalArgumentException ("m is invalid: 1 <= m < l");
      if (k <= 0 || k > l)
         throw new IllegalArgumentException ("k is invalid: 1 <= k < l");
      int imin = Math.max (0, k - l + m);
      int imax = Math.min (k, m);
      if (x < imin)
         return 0.0;
      if (x >= imax)
         return 1.0;
      // Very inefficient
      double res = 0.0;
      for (int i = imin; i <= x; i++)
         res += prob (m, l, k, i);
      if (res >= 1.0)
         return 1.0;
      return res;
   }

   /**
    * Computes the complementary distribution function. *WARNING:* The
    * complementary distribution function is defined as @f$\bar{F}(x) =
    * P[X \ge x]@f$.
    */
   public static double barF (int m, int l, int k, int x) {
      if (l <= 0)
         throw new IllegalArgumentException ("l must be greater than 0");
      if (m <= 0 || m > l)
         throw new IllegalArgumentException ("m is invalid: 1 < =m < l");
      if (k <= 0 || k > l)
         throw new IllegalArgumentException ("k is invalid: 1 < =k < l");
      int imin = Math.max (0, k - l + m);
      int imax = Math.min (k, m);
      if (x <= imin)
         return 1.0;
      if (x > imax)
         return 0.0;
      // Very inefficient
      double res = 0.0;
      for (int i = imax; i >= x; i--)
         res += prob (m, l, k, i);
      if (res >= 1.0)
         return 1.0;
      return res;
   }

   /**
    * Computes @f$F^{-1}(u)@f$ for the hypergeometric distribution without
    * using precomputed tables. The inversion is computed using the
    * chop-down algorithm @cite sKAC85a&thinsp;.
    */
   public static int inverseF (int m, int l, int k, double u) {
      // algo hin dans Kachitvichyanukul
      if (u < 0 || u >= 1)
         throw new IllegalArgumentException ("u is not in [0,1]");
      if (u <= 0.0)
         return Math.max (0, k - l + m);
      if (u >= 1.0)
         return Math.min (k, m);
      double p = 0;
      if (k < l - m)
          p = Math.exp (Num.lnFactorial (l - m) + Num.lnFactorial (l - k)
                       -Num.lnFactorial (l) - Num.lnFactorial (l - m - k));
      else
          p = Math.exp (Num.lnFactorial (m) + Num.lnFactorial (k)
                      -Num.lnFactorial (k- l + m) - Num.lnFactorial (l));

     // Empirical correction, the original algorithm sets x=0.
     int x = Math.max (0, k - l + m);
     if (u <= p) return x;

     do {
         u = u - p;
         p = p*(m - x)*(k - x)/((x+1)*(l - m - k + 1.0 + x));
         x++;
     } while (u > p && p > 0);
     return x;
   }

   /**
    * Computes and returns the mean @f$E[X] = km/l@f$ of the
    * Hypergeometric distribution with parameters @f$m@f$, @f$l@f$ and
    * @f$k@f$.
    *  @return the mean of the hypergeometric distribution @f$E[X] = km /
    * l@f$
    */
   public static double getMean (int m, int l, int k) {
      if (l <= 0)
         throw new IllegalArgumentException ("l must be greater than 0");
      if (m <= 0 || m > l)
         throw new IllegalArgumentException ("m is invalid: 1<=m<l");
      if (k <= 0 || k > l)
         throw new IllegalArgumentException ("k is invalid: 1<=k<l");

      return ((double) k *  (double) m / (double) l);
   }

   /**
    * Computes and returns the variance
    *  @f$\mbox{Var}[X] = \frac{(km/l)(1 - m/l)(l - k)}{l - 1}@f$
    *  of the hypergeometric distribution with parameters @f$m@f$, @f$l@f$
    * and @f$k@f$.
    *  @return the variance of the Hypergeometric distribution
    * @f$\mbox{Var}[X] = (km / l)(1 - m / l)(l - k) / (l - 1)@f$
    */
   public static double getVariance (int m, int l, int k) {
      if (l <= 0)
         throw new IllegalArgumentException ("l must be greater than 0");
      if (m <= 0 || m > l)
         throw new IllegalArgumentException ("m is invalid: 1<=m<l");
      if (k <= 0 || k > l)
         throw new IllegalArgumentException ("k is invalid: 1<=k<l");

      return (((double) k * (double) m / (double) l) *
              ( 1 - ((double) m / (double) l)) * ((double) l - (double) k) /
              ((double) l - 1));
   }

   /**
    * Computes and returns the standard deviation of the hypergeometric
    * distribution with parameters @f$m@f$, @f$l@f$ and @f$k@f$.
    *  @return the standard deviation of the hypergeometric distribution
    */
   public static double getStandardDeviation (int m, int l, int k) {
      return Math.sqrt (HypergeometricDist.getVariance (m, l, k));
   }

   /**
    * Returns the @f$m@f$ associated with this object.
    */
   public int getM() {
      return m;
   }

   /**
    * Returns the @f$l@f$ associated with this object.
    */
   public int getL() {
      return l;
   }

   /**
    * Returns the @f$k@f$ associated with this object.
    */
   public int getK() {
      return k;
   }

   private void setHypergeometric() {
      int imin = Math.max (0, k - l + m);
      int imax = Math.min (k, m);
      supportA = imin;
      supportB = imax;
      int ns = imax - imin + 1;
      if (ns > MAXN) {
         pdf = null;
         cdf = null;
         return;
      }

      int offset = imin;
      imin = 0;
      imax -= offset;
      double[] P = new double[ns];
      double[] F = new double[ns];

      // Computes the mode (taken from UNURAN)
      int mode = (int)((k + 1.0)*(m + 1.0)/(l + 2.0));
      int imid = mode - offset;

      P[imid] = prob (m, l, k, mode);

      int i = imid;
      while (i > imin && Math.abs (P[i]) > EPSILON) {
         P[i-1] = P[i]*(i + offset)/(m-i-offset+1)
                  * (l - m - k + i + offset)/(k - i - offset + 1);
         i--;
      }
      imin = i;

      i = imid;
      while (i < imax && Math.abs (P[i]) > EPSILON) {
         P[i+1] = P[i]*(m - i - offset)/(i + offset + 1)
                  * (k - i - offset)/(l - m - k + i + offset + 1);
         i++;
      }
      imax = i;

      F[imin] = P[imin];
      i = imin;
      while (i < imax && F[i] < 0.5) {
         i++;
         F[i] = F[i-1] + P[i];
      }
      xmed = i;

      F[imax] = P[imax];
      i = imax - 1;
      while (i > xmed) {
         F[i] = P[i] + F[i + 1];
         i--;
      }

       xmin = imin + offset;
       xmax = imax + offset;
       xmed += offset;
       pdf  = new double[imax + 1 - imin];
       cdf  = new double[imax + 1 - imin];
       System.arraycopy (P, imin, pdf, 0, imax+1-imin);
       System.arraycopy (F, imin, cdf, 0, imax+1-imin);
   }

/**
 * Return a table containing the parameters of the current distribution. This
 * table is put in regular order: [@f$m@f$, @f$l@f$, @f$k@f$].
 */
public double[] getParams () {
      double[] retour = {m, l, k};
      return retour;
   }

   /**
    * Resets the parameters of this object to @f$m@f$, @f$l@f$ and
    * @f$k@f$.
    */
   public void setParams (int m, int l, int k) {
      if (l <= 0)
         throw new IllegalArgumentException ("l must be greater than 0");
      if (m <= 0 || m > l)
         throw new IllegalArgumentException ("m is invalid: 1 <= m < l");
      if (k <= 0 || k > l)
         throw new IllegalArgumentException ("k is invalid: 1 <= k < l");
      this.m = m;
      this.l = l;
      this.k = k;
      setHypergeometric();
      if (k < l - m)
          p0 = Math.exp (Num.lnFactorial (l - m) + Num.lnFactorial (l - k)
                               -Num.lnFactorial (l) - Num.lnFactorial (l - m - k));
      else
          p0 = Math.exp (Num.lnFactorial (m) + Num.lnFactorial (k)
                               -Num.lnFactorial (k- l + m) - Num.lnFactorial (l));
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : m = " + m + ", l = " + l + ", k = " + k;
   }

}