/*
 * Class:        BinomialDist
 * Description:  binomial distribution
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
 * Extends the class  @ref DiscreteDistributionInt for the *binomial*
 * distribution @cite sLAW00a&thinsp; (page 321) with parameters @f$n@f$ and
 * @f$p@f$, where @f$n@f$ is a positive integer and @f$0\le p\le1@f$. Its
 * mass function is given by
 * @anchor REF_probdist_BinomialDist_eq_fmass_binomial
 * @f[
 *   p(x) = \binom{n}{x} p^x (1-p)^{n-x} = \frac{n!}{x!(n-x)!}\; p^x (1-p)^{n-x} \qquad\mbox{for } x = 0, 1, 2,…n, \tag{fmass-binomial}
 * @f]
 * and its distribution function is
 * @anchor REF_probdist_BinomialDist_eq_Fbinomial
 * @f[
 *   F (x) = \sum_{j=0}^x \binom{n}{j}\; p^j (1-p)^{n-j} \qquad\mbox{for } x = 0, 1, 2, …n, \tag{Fbinomial}
 * @f]
 * where nCr@f$(n,x)@f$ is the number of possible combinations of @f$x@f$
 * elements chosen among a set of @f$n@f$ elements.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_discrete
 */
public class BinomialDist extends DiscreteDistributionInt {
   private int n;
   private double p;
   private double q;
   private static final double EPS2 = 100.0*EPSILON;

   private static class Function implements MathFunction {
      protected int m;
      protected int R;
      protected double mean;
      protected int f[];

      public Function (int m, double mean, int R, int f[]) {
         this.m = m;
         this.mean = mean;
         this.R = R;
         this.f = new int[f.length];
         System.arraycopy (f, 0, this.f, 0, f.length);
      }


      public double evaluate (double x) {
         if (x < R)
            return 1e100;

         double sum = 0.0;
         for (int j = 0; j < R; j++)
            sum += f[j] / (x - (double) j);

         return (sum + m * Math.log1p (-mean / x));
      }
   }

   /**
    * @name Constant
    * @{
    */

   /**
    * The value of the parameter @f$n@f$ above which the tables are *not*
    * precomputed by the constructor.
    */
   public static double MAXN = 100000;

   /**
    * @}
    */

   /**
    * Creates an object that contains the binomial terms (
    * {@link REF_probdist_BinomialDist_eq_fmass_binomial
    * fmass-binomial} ), for @f$0\le x\le n@f$, and the corresponding
    * cumulative function. These values are computed and stored in dynamic
    * arrays, unless @f$n@f$ exceeds `MAXN`.
    */
   public BinomialDist (int n, double p) {
      setBinomial (n, p);
   }

   public double prob (int x) {
      if (n == 0)
         return 1.0;

      if (x < 0 || x > n)
         return 0.0;

      if (p == 0.0) {
         if (x > 0)
            return 0.0;
         else
            return 1.0;
      }

      if (q == 0.0) {
         if (x < n)
            return 0.0;
         else
            return 1.0;
      }

      if (pdf == null)
         return prob (n, p, q, x);

      if (x > xmax || x < xmin)
         return prob (n, p, q, x);

      return pdf[x - xmin];
   }


   public double cdf (int x) {
      if (n == 0)
         return 1.0;

      if (x < 0)
         return 0.0;

      if (x >= n)
         return 1.0;

      if (p == 0.0)
         return 1.0;

      if (p == 1.0)
         return 0.0;

      if (cdf != null) {
         if (x >= xmax)
            return 1.0;
         if (x < xmin) {
            final int RMAX = 20;
            int i;
            double term = prob(x);
            double Sum = term;
            final double z = (1.0 - p) / p;
            i = x;
            while (i > 0 && i >= x - RMAX) {
               term *= z * i / (n - i + 1);
               i--;
               Sum += term;
            }
            return Sum;
         }
         if (x <= xmed)
            return cdf[x - xmin];
         else
            // We keep the complementary distribution in the upper part of cdf
            return 1.0 - cdf[x + 1 - xmin];
      } else
         return cdf (n, p, x);
   }


   public double barF (int x) {
      if (n == 0)
         return 1.0;

      if (x < 1)
         return 1.0;

      if (x > n)
         return 0.0;

      if (p == 0.0)
         return 0.0;

      if (p == 1.0)
         return 1.0;

      if (cdf != null) {
         if (x > xmax) {
            // Add IMAX dominant terms to get a few decimals in the tail
            final double q = 1.0 - p;
            double z, sum, term;
            int i;
            sum = term = prob(x);
            z = p / q;
            i = x;
            final int IMAX = 20;
            while (i < n && i < x + IMAX) {
               term = term * z * (n - i) / (i + 1);
               sum += term;
               i++;
            }
            return sum;
            // return fdist_Beta (x, n - x + 1, 10, p);
         }

         if (x <= xmin)
            return 1.0;
         if (x > xmed)
            // We keep the complementary distribution in the upper part of cdf
            return cdf[x - xmin];
         else
            return 1.0 - cdf[x - 1 - xmin];
      } else
         return 1.0 - cdf (n, p, x - 1);
   }


   public int inverseFInt (double u) {
      if ((cdf == null) || (u <= EPS2))
         return inverseF (n, p, u);
      else
         return super.inverseFInt (u);
   }

   public double getMean() {
      return BinomialDist.getMean (n, p);
   }

   public double getVariance() {
      return BinomialDist.getVariance (n, p);
   }

   public double getStandardDeviation() {
      return BinomialDist.getStandardDeviation (n, p);
   }

/**
 * Computes and returns the binomial probability @f$p(x)@f$ in eq. (
 * {@link REF_probdist_BinomialDist_eq_fmass_binomial
 * fmass-binomial} ).
 */
public static double prob (int n, double p, int x) {
      return prob (n, p, 1.0 - p, x);
   }

   /**
    * A generalization of the previous method. Computes and returns the
    * binomial term
    * @anchor REF_probdist_BinomialDist_eq_fmass_binom
    * @f[
    *   f(x) = \binom{n}{x} p^x q^{n-x} = \frac{n!}{x!(n-x)!}\; p^x q^{n-x}, \tag{fmass-binom}
    * @f]
    * where @f$p@f$ and @f$q@f$ are arbitrary real numbers (@f$q@f$ is not
    * necessarily equal to @f$1-p@f$). In the case where @f$0 \le p
    * \le1@f$ and @f$q = 1-p@f$, the returned value is a probability term
    * for the binomial distribution.
    */
   public static double prob (int n, double p, double q, int x) {
      final int SLIM = 50;          // To avoid overflow
      final double MAXEXP = (Num.DBL_MAX_EXP - 1)*Num.LN2;// To avoid overflow
      final double MINEXP = (Num.DBL_MIN_EXP - 1)*Num.LN2;// To avoid underflow
      int signe = 1;
      double Res;

      if (n < 0)
        throw new IllegalArgumentException ("n < 0");
      if (n == 0)
         return 1.0;
      if (x < 0 || x > n)
         return 0.0;

      // Combination (n, x) are symmetric between x and n-x
      if (x > n/2) {
         x = n - x;
         Res = p;
         p = q;
         q = Res;
      }

      if (p < 0.0) {
         p = -p;
         if ((x & 1) != 0)
            signe *= -1;             // odd x
      }
      if (q < 0.0) {
         q = -q;
         if (((n - x) & 1) != 0)
            signe *= -1;             // odd n - x
      }

      if (n <= SLIM) {
         Res = Math.pow (p, (double)x)*Num.combination (n, x)*Math.pow (q,
            (double)(n - x));
         return signe*Res;
      }
      else {
         /* This could be calculated with more precision as there is some
            cancellation because of subtraction of the large LnFactorial: the
            last few digits can be lost. But we need the function lgammal in
            long double precision. Another possibility would be to use an
            asymptotic expansion for the binomial coefficient. */

         Res = x*Math.log (p) + (n - x)*Math.log (q) + Num.lnFactorial (n)
            - Num.lnFactorial (n - x) - Num.lnFactorial (x);
         if (Res >= MAXEXP)
           throw new IllegalArgumentException ("term overflow");

         if (Res < MINEXP)
            return 0.0;

         return signe*Math.exp (Res);
      }
   }

   /**
    * Computes @f$F(x)@f$, the distribution function of a binomial random
    * variable with parameters @f$n@f$ and @f$p@f$, evaluated at @f$x@f$.
    * If @f$n \le10000@f$, the non-negligible terms of the sum are added
    * explicitly. If @f$n > 10000@f$ and @f$np (1-p) > 100@f$, the
    * Camp-Paulson normal approximation @cite tCAM51a,
    * @cite tMOL70a&thinsp; is used, otherwise, a Poisson approximation
    * due to Bol’shev @cite tBOL64a, @cite tMOL70a&thinsp; is used.
    */
   public static double cdf (int n, double p, int x) {
      final int NLIM1 = 10000;
      final double VARLIM = 50.0;
//      final double EPS = DiscreteDistributionInt.EPSILON * EPS_EXTRA;
      double y, z, q = 1.0 - p;
      double sum, term, termmid;
      int i, mid;
      boolean flag = false;

      if (p < 0.0 | p > 1.0)
        throw new IllegalArgumentException ("p not in [0,1]");
      if (n < 0)
        throw new IllegalArgumentException ("n < 0");

      if (n == 0)
         return 1.0;
      if (x < 0)
         return 0.0;
      if (x >= n)
         return 1.0;
      if (p <= 0.0)
         return 1.0;
      if (p >= 1.0)
         return 0.0;                 // For any x < n

      if (n < NLIM1) {               // Exact Binomial
         /* Sum RMAX terms to get a few decimals in the lower tail */
         final int RMAX = 20;
         mid = (int)((n + 1)*p);
         if (mid > x)
            mid = x;
         sum = term = termmid = prob (n, p, 1.0 - p, mid);

         z = q/p;
            i = mid;
            while (term >= EPSILON || i >= mid - RMAX) {
            term *= z * i / (n - i + 1);
            sum += term;
            i--;
            if (i == 0) break;
         }

         z = p/q;
         term = termmid;
         for (i = mid; i < x; i++) {
            term *= z*(n - i)/(i + 1);
            if (term < EPSILON)
               break;
            sum += term;
         }
         if (sum >= 1.0) return 1.0;
         return sum;

      } else {
         if (p > 0.5 || ((p == 0.5) && (x > n/2))) {
            // use F (p, n, x) = 1 - F (q, n, n-x-1)
            p = q;
            q = 1.0 - p;
            flag = true;
            x = n - x - 1;
         }
         if (n*p*q > VARLIM) {   // Normal approximation
            /* Uses the Camp-Paulson approximation based on the F-distribution.
               Its maximum absolute error is smaller than 0.007 / sqrt (npq).
               Ref: W. Molenaar; Approximations to the Poisson, Binomial,....
               QA273.6 M64, p. 93 (1970) */
            term = Math.pow ((x+1)*q/((n-x)*p), 1.0/3.0);
            y = term*(9.0 - 1.0/(x + 1)) - 9.0 + 1.0/(n - x);
            z = 3.0*Math.sqrt (term*term/(x + 1) + 1.0/(n - x));
            y /= z;
            if (flag)
               return NormalDist.barF01 (y);
            else
               return NormalDist.cdf01 (y);
         }
         else {                    // Poisson approximation
            /* Uses a Bol'shev approximation based on the Poisson distribution.
               Error is O (1/n^4) as n -> infinity. Ref: W. Molenaar;
             Approximations to the Poisson, Binomial,... QA273.6 M64, p. 107,
               Table 6.2, Formule lambda_9 (1970). */
            y = (2.0*n - x)*p/(2.0 - p);
            double t = 2.0*n - x;
            z = (2.0*y*y - x*y - x*(double)x - 2.0*x)/(6*t*t);
            z = y/(1.0 - z);
            if (flag)
               return PoissonDist.barF(z, x - 1);
            else
               return PoissonDist.cdf (z, x);
         }
      }
   }

   /**
    * Returns @f$\bar{F}(x) = P[X \ge x]@f$, the complementary
    * distribution function. The default implementation returns `1 -
    * cdf(x-1)`, which is not accurate when @f$F(x)@f$ is close to 1.
    */
   public static double barF (int n, double p, int x) {
      return 1.0 - cdf (n, p, x - 1);
   }

   /**
    * Computes @f$x = F^{-1}(u)@f$, the inverse of the binomial
    * distribution. If @f$n@f$ is larger than 10000, the linear search
    * starts from 0 and the  #cdf(int,double,int) static method is used to
    * compute @f$F(x)@f$ at different values of @f$x@f$, which is much
    * less efficient.
    */
   public static int inverseF (int n, double p, double u) {
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u not in [0,1]");
      if (u <= prob (n, p, 0))
         return 0;
      if ((u > 1.0 - prob (n, p, n)) || (u >= 1.0))
         return n;
      final int NLIM1 = 10000;
      int i, mid;
      final double q = 1.0 - p;
      double z, sum, term, termmid;

      if (n < NLIM1) {               // Exact Binomial
         i = (int)((n + 1)*p);
         if (i > n)
            i = n;
         term = prob (n, p, i);
         while ((term >= u) && (term > Double.MIN_NORMAL)) {
            i /= 2;
            term = prob (n, p, i);
         }

         z = q/p;
         if (term <= Double.MIN_NORMAL) {
            i *= 2;
            term = prob (n, p, i);
            while (term >= u && (term > Double.MIN_NORMAL)) {
               term *= z*i/(n - i + 1);
               i--;
            }
         }

         mid = i;
         term = prob (n, p, i);
         sum = termmid = term;

         z = q/p;
         for (i = mid; i > 0; i--) {
            term *= z*i/(n - i + 1);
            if (term < EPSILON)
               break;
            sum += term;
         }

         i = mid ;
         term = termmid;
         if (sum >= u) {
             while (sum >= u){
                 z = q/p;
                 sum -= term;
                 term *= z*i/(n - i + 1);
                 i--;
             }
             return i+1;
         }

        double prev = -1;
        while ((sum < u) && (sum > prev)) {
            z = p/q;
            term *=  z*(n - i)/(i + 1);
            prev = sum;
            sum += term;
            i++;
         }
         return i;

      } else{  // Normal or Poisson  approximation
          for(i = 0 ; i <= n ; i++){
              sum = cdf(n,p,i);
              if (sum >= u)
                  break;
          }
          return i;
      }
    }

   /**
    * Estimates the parameters @f$(n,p)@f$ of the binomial distribution
    * using the maximum likelihood method, from the @f$m@f$ observations
    * @f$x[i]@f$, @f$i = 0, 1,…, m-1@f$. The estimates are returned in a
    * two-element array, in regular order: [@f$n@f$, @f$p@f$].  The
    * maximum likelihood estimators are the values @f$(\hat{n},
    * \hat{p})@f$ that satisfy the equations:
    * @f{align*}{
    *    \hat{n}\hat{p} 
    *    & 
    *    = 
    *    \bar{x}_m
    *    \\ 
    *   \sum_{i=0}^{R-1} \frac{f_i}{\hat{n} - i} 
    *    & 
    *    = 
    *    - m \ln\left(1 - \frac{\bar{x}_m}{\hat{n}}\right)
    * @f}
    * where @f$\bar{x}_m@f$ is the average of @f$x[0],…,x[m-1]@f$,
    * @f$f_i@f$ is the number of @f$x[i]@f$’s that exceed @f$i@f$, and
    * @f$R =\max(x[0],…,x[m-1])@f$ is the largest observation
    * @cite tJOH69a&thinsp; (page 57).
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param m            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{n}@f$, @f$\hat{p}@f$]
    */
   public static double[] getMLE (int[] x, int m) {
      if (m <= 1)
         throw new UnsupportedOperationException(" m < 2");

      int i;
      int r = 0;
      double mean = 0.0;
      for (i = 0; i < m; i++) {
         mean += x[i];
         if (x[i] > r)
            r = x[i];
      }
      mean /= (double) m;

      double sum = 0.0;
      for (i = 0; i < m; i++)
         sum += (x[i] - mean)*(x[i] - mean);
      double var = sum / m;
      if (mean <= var)
         throw new UnsupportedOperationException("mean <= variance");

      int f[] = new int[r];
      for (int j = 0; j < r; j++) {
         f[j] = 0;
         for (i = 0; i < m; i++)
            if (x[i] > j) f[j]++;
      }

      double p = 1.0 - var/mean;
      double rup = (int) (5*mean/p);
      if (rup < 1) rup = 1;

      Function fct = new Function (m, mean, r, f);
      double parameters[] = new double[2];
      parameters[0] = (int) RootFinder.brentDekker (r-1, rup, fct, 1e-5);
      if (parameters[0] < r)
         parameters[0] = r;
      parameters[1] = mean / parameters[0];

      return parameters;
   }

   /**
    * Creates a new instance of a binomial distribution with both
    * parameters @f$n@f$ and @f$p@f$ estimated using the maximum
    * likelihood method, from the @f$m@f$ observations @f$x[i]@f$, @f$i =
    * 0, 1, …, m-1@f$.
    *  @param x            the list of observations to use to estimate the
    *                      parameters
    *  @param m            the number of observations to use to estimate
    *                      the parameters
    */
   public static BinomialDist getInstanceFromMLE (int[] x, int m) {
      double parameters[] = new double[2];
      parameters = getMLE (x, m);
      return new BinomialDist ((int) parameters[0], parameters[1]);
   }

   /**
    * Estimates the parameter @f$p@f$ of the binomial distribution with
    * given (fixed) parameter @f$n@f$, by the maximum likelihood method,
    * from the @f$m@f$ observations @f$x[i]@f$, @f$i = 0, 1,…, m-1@f$.
    * Returns the estimator in an array with a single element.
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param m            the number of observations used to evaluate
    *                      parameters
    *  @param n            the number of success
    *  @return returns the parameter [@f$\hat{p}@f$]
    */
   public static double[] getMLE (int[] x, int m, int n) {
      if (m <= 0)
         throw new IllegalArgumentException ("m <= 0");
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double parameters[] = new double[1];
      double mean = 0.0;
      for (int i = 0; i < m; i++)
         mean += x[i];
      mean /= (double) m;

      parameters[0] = mean / (double) n;
      return parameters;
   }

   /**
    * Creates a new instance of a binomial distribution with given (fixed)
    * parameter @f$n@f$, and with parameter @f$p@f$ estimated by the
    * maximum likelihood method based on the @f$m@f$ observations
    * @f$x[i]@f$, @f$i = 0, 1, …, m-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param m            the number of observations to use to evaluate
    *                      parameters
    *  @param n            the parameter n of the binomial
    */
   public static BinomialDist getInstanceFromMLE (int[] x, int m, int n) {
      double parameters[] = new double[1];
      parameters = getMLE (x, m, n);
      return new BinomialDist (n, parameters[0]);
   }

   /**
    * Computes the mean @f$E[X] = np@f$ of the binomial distribution with
    * parameters @f$n@f$ and @f$p@f$.
    *  @return the mean of the Binomial distribution @f$E[X] = np@f$
    */
   public static double getMean (int n, double p) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in range (0, 1)");

      return (n * p);
   }

   /**
    * Computes the variance @f$\mbox{Var}[X] = np(1 - p)@f$ of the
    * binomial distribution with parameters @f$n@f$ and @f$p@f$.
    *  @return the variance of the binomial distribution @f$\mbox{Var}[X]
    * = np(1 - p)@f$
    */
   public static double getVariance (int n, double p) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in range (0, 1)");

      return (n * p * (1 - p));
   }

   /**
    * Computes the standard deviation of the Binomial distribution with
    * parameters @f$n@f$ and @f$p@f$.
    *  @return the standard deviation of the binomial distribution
    */
   public static double getStandardDeviation (int n, double p) {
      return Math.sqrt (BinomialDist.getVariance (n, p));
   }

   private void setBinomial (int n, double p) {
     /*
      * Compute all probability terms of the binomial distribution; start near
      * the mean, and calculate probabilities on each side until they become
      * smaller than EPSILON, then stop there.
      * However, this is more general than the binomial probability distribu-
      * tion as this will compute the binomial terms when p + q != 1, and
      * even when p or q are negative. However in this case, the cumulative
      * terms will be meaningless.
      */

      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in range (0, 1)");
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      supportA = 0;
      supportB = n;
      double q = 1.0 - p;
      final double EPS = DiscreteDistributionInt.EPSILON * EPS_EXTRA;

      this.n = n;
      this.p = p;
      this.q = q;

      int i, mid;
      int imin, imax;
      double z = 0.0;
      double[] P;     // Binomial "probability" terms
      double[] F;     // Binomial cumulative "probabilities"

      // For n > MAXN, we shall not use pre-computed arrays
      if (n > MAXN) {
         pdf = null;
         cdf = null;
         return;
      }

      P = new double[1 + n];
      F = new double[1 + n];

      // the maximum term in absolute value
      mid = (int)((n + 1)*Math.abs (p)/(Math.abs (p) + Math.abs (q)));
      if (mid > n)
         mid = n;
      P[mid] = prob (n, p, q, mid);

      if (p != 0.0 || p != -0.0)
         z = q/p;
      i = mid;

      while (i > 0 && Math.abs (P[i]) > EPS) {
         P[i - 1] = P[i]*z*i/(n - i + 1);
         i--;
      }
      imin = i;

      if (q != 0.0 || q != -0.0)
         z = p/q;
      i = mid;

      while (i < n && Math.abs (P[i]) > EPS) {
         P[i + 1] = P[i]*z*(n - i)/(i + 1);
         i++;
      }
      imax = i;

   /* Here, we assume that we are dealing with a probability distribution.
      Compute the cumulative probabilities for F and keep them in the
      lower part of CDF.*/
      F[imin] = P[imin];
      i = imin;
      while (i < n && F[i] < 0.5) {
         i++;
         F[i] = F[i - 1] + P[i];
      }

      // This is the boundary between F (i <= xmed) and 1 - F (i > xmed) in
      // the array CDF
      xmed = i;

      // Compute the cumulative probabilities of the complementary
      // distribution and keep them in the upper part of the array
      F[imax] = P[imax];
      i = imax - 1;
      while (i > xmed) {
         F[i] = P[i] + F[i + 1];
         i--;
      }

       // Reset imin because we lose too much precision for a few terms near
       //   imin when we stop adding terms < epsilon.
       i = imin;
       while (i < xmed && F[i] < DiscreteDistributionInt.EPSILON)
          i++;
       xmin = imin = i;

       /* Same thing with imax */
       i = imax;
       while (i > xmed && F[i] < DiscreteDistributionInt.EPSILON)
          i--;
       xmax = imax = i;

       pdf  = new double[imax + 1 - imin];
       cdf  = new double[imax + 1 - imin];
       System.arraycopy (P, imin, pdf, 0, imax+1-imin);
       System.arraycopy (F, imin, cdf, 0, imax+1-imin);

   }

/**
 * Returns the parameter @f$n@f$ of this object.
 */
public int getN() {
      return n;
   }

   /**
    * Returns the parameter @f$p@f$ of this object.
    */
   public double getP() {
      return p;
   }

   /**
    * Returns a table that contains the parameters @f$(n,p)@f$ of the
    * current distribution, in regular order: [@f$n@f$, @f$p@f$].
    */
   public double[] getParams () {
      double[] retour = {n, p};
      return retour;
   }

   /**
    * Resets the parameters to these new values and recomputes everything
    * as in the constructor. From the performance viewpoint, it is
    * essentially the same as constructing a new  @ref BinomialDist
    * object.
    */
   public void setParams (int n, double p) {
      setBinomial (n, p);
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : n = " + n + ", p = " + p;
   }

}