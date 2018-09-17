/*
 * Class:        NegativeBinomialDist
 * Description:  negative binomial distribution
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
import optimization.*;

/**
 * Extends the class  @ref DiscreteDistributionInt for the *negative
 * binomial* distribution @cite sLAW00a&thinsp; (page 324) with real
 * parameters @f$n@f$ and @f$p@f$, where @f$n > 0@f$ and @f$0\le p\le1@f$.
 * Its mass function is
 * @anchor REF_probdist_NegativeBinomialDist_eq_fmass_negbin
 * @f[
 *   p(x) = \frac{\Gamma(n + x)}{\Gamma(n)\; x!} p^n (1 - p)^x, \qquad\mbox{for } x = 0, 1, 2, …\tag{fmass-negbin}
 * @f]
 * where @f$\Gamma(x)@f$ is the gamma function.
 *
 * If @f$n@f$ is an integer, @f$p(x)@f$ can be interpreted as the probability
 * of having @f$x@f$ failures before the @f$n@f$-th success in a sequence of
 * independent Bernoulli trials with probability of success @f$p@f$. This
 * special case is implemented as the Pascal distribution (see
 * @ref PascalDist ).
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_discrete
 */
public class NegativeBinomialDist extends DiscreteDistributionInt {
   protected double n;
   protected double p;
   private static final double EPS2 = 1000.0*EPSILON;

   private static class Func1 implements MathFunction {
      protected int m;
      protected int[] x;
      protected double p;

      public Func1 (double p, int[] x, int m) {
         this.p = p;
         this.m = m;
         this.x = x;
      }

      public double evaluate (double gam) {
         if (gam <= 0 ) return 1.0e100;

         double sum = 0.0;
         for (int j = 0; j < m; j++)
            sum += Num.digamma (gam + x[j]);
         return sum/m + Math.log (p) - Num.digamma (gam);
      }
   }


   private static class Function implements MathFunction {
      protected int m;
      protected int max;
      protected double mean;
      protected int[] Fj;

      public Function (int m, int max, double mean, int[] Fj) {
         this.m = m;
         this.max = max;
         this.mean = mean;
         this.Fj = new int[Fj.length];
         System.arraycopy(Fj, 0, this.Fj, 0, Fj.length);
      }

      public double evaluate (double s) {
    //     if (s <= 0 ) return 1.0e100;
         double sum = 0.0;
         double p = s / (s + mean);

         for (int j = 0; j < max; j++)
            sum += Fj[j] / (s + (double) j);

         return sum + m * Math.log (p);
      }
   }


   private static class FuncInv extends Function implements MathFunction {

      public FuncInv (int m, int max, double mean, int[] Fj) {
         super (m, max, mean, Fj);
      }

      public double evaluate (double nu) {
         double r = nu*mean;
         double sum = 0.;
         for (int j = 0; j < max; j++)
            sum += Fj[j] / (1.0 + nu* j);

         return (nu*sum - m * Math.log1p (r));
      }
   }


/************************
   // Class Optim seems to be useless
   private static class Optim implements Lmder_fcn
   {
      private double mean;
      private int N;
      private int max;
      private int [] Fj;

      public Optim (int N, int max, double mean, int[] Fj)
      {
         this.N = N;
         this.max = max;
         this.mean = mean;
         this.Fj = new int[max];
         System.arraycopy (Fj, 0, this.Fj, 0, max);
      }

      public void fcn (int m, int n, double[] x, double[] fvec, double[][] fjac,
                       int iflag[])
      {
         if (x[1] <= 0.0 || x[2] <= 0.0 || x[2] >= 1.0) {
             final double BIG = 1.0e100;
             fvec[1] = BIG;
             fvec[2] = BIG;
             fjac[1][1] = BIG;
             fjac[1][2] = 0.0;
             fjac[2][1] = 0.0;
             fjac[2][2] = BIG;
             return;
         }

         double trig;
         if (iflag[1] == 1)
         {
            double sum = 0.0;
            for (int j = 0; j < max; j++)
               sum += Fj[j] / (x[1] + j);
            fvec[1] = x[2] * mean  - x[1] * (1.0 - x[2]);
            fvec[2] =  N * Math.log (x[2]) + sum;

         } else if (iflag[1] == 2) {

            fjac[1][1] = x[2] - 1.0;
            fjac[1][2] = mean + x[1];
            double sum = 0.0;
            for (int j = 0; j < max; j++)
               sum += Fj[j] / ((x[1] + j)*(x[1] + j));
            fjac[2][1] = -sum;
            fjac[2][2] = N / x[2];
         }
      }
   }
****************************/

   /**
    * @name Constant
    * @{
    */

   /**
    * If the maximum term is greater than this constant, then the tables
    * will *not* be precomputed.
    */
   public static double MAXN = 100000;

/**
 * @}
 */

   protected NegativeBinomialDist() {}

/**
 * Creates an object that contains the probability terms (
 * {@link REF_probdist_NegativeBinomialDist_eq_fmass_negbin
 * fmass-negbin} ) and the distribution function for the negative binomial
 * distribution with parameters @f$n@f$ and @f$p@f$.
 */
public NegativeBinomialDist (double n, double p) {
      setParams (n, p);
   }

   public double prob (int x) {
      if (x < 0)
         return 0.0;

      if (p == 0.0)
         return 0.0;

      if (p == 1.0) {
         if (x > 0)
            return 0.0;
         else
            return 1.0;
      }

      if (pdf == null)
         return prob (n, p, x);

      if (x > xmax || x < xmin)
         return prob (n, p, x);

      return pdf[x - xmin];
   }

   public double cdf (int x) {
      if (x < 0)
         return 0.0;
      if (p >= 1.0)    // In fact, p == 1
         return 1.0;
      if (p <= 0.0)    // In fact, p == 0
         return 0.0;

      if (cdf != null) {
         if (x >= xmax)
            return 1.0;
         if (x < xmin)
            return cdf (n, p, x);
         if (x <= xmed)
            return cdf[x - xmin];
         else
            // We keep the complementary distribution in the upper part of cdf
            return 1.0 - cdf[x + 1 - xmin];

      }
      else
         return cdf (n, p, x);
   }

   public double barF (int x) {
      if (x < 1)
         return 1.0;
      if (p >= 1.0)   // In fact, p == 1
         return 0.0;
      if (p <= 0.0)   // In fact, p == 0
         return 1.0;

      if (cdf == null)
         //return BinomialDist.cdf (x - 1 + n, p, n - 1);
         return BetaDist.barF (n, x, p);

      if (x > xmax)
         //return BinomialDist.cdf (x - 1 + n, p, n - 1);
         return BetaDist.barF (n, x, p);

      if (x <= xmin)
         return 1.0;
      if (x > xmed)
         // We keep the complementary distribution in the upper part of cdf
         return cdf[x - xmin];
      else
         return 1.0 - cdf[x - 1 - xmin];
   }

   public int inverseFInt (double u) {
      if ((cdf == null) || (u <= EPS2))
         return inverseF (n, p, u);
      else
         return super.inverseFInt (u);
   }

   public double getMean() {
      return NegativeBinomialDist.getMean (n, p);
   }

   public double getVariance() {
      return NegativeBinomialDist.getVariance (n, p);
   }

   public double getStandardDeviation() {
      return NegativeBinomialDist.getStandardDeviation (n, p);
   }

/**
 * Computes the probability @f$p(x)@f$ defined in (
 * {@link REF_probdist_NegativeBinomialDist_eq_fmass_negbin
 * fmass-negbin} ).
 */
public static double prob (double n, double p, int x) {
      final int SLIM = 15;           // To avoid overflow
      final double MAXEXP = (Num.DBL_MAX_EXP - 1)*Num.LN2;// To avoid overflow
      final double MINEXP = (Num.DBL_MIN_EXP - 1)*Num.LN2;// To avoid underflow
      double y;

      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in [0, 1]");
      if (n <= 0.0)
         throw new IllegalArgumentException ("n <= 0.0");
      if (x < 0)
         return 0.0;
      if (p >= 1.0) {                // In fact, p == 1
         if (x == 0)
            return 1.0;
         else
            return 0.0;
      }
      if (p <= 0.0)                  // In fact, p == 0
         return 0.0;

      y = Num.lnGamma (n + x) - (Num.lnFactorial (x) + Num.lnGamma (n))
          + n * Math.log (p) + x * Math.log1p (-p) ;

      if (y >= MAXEXP)
         throw new IllegalArgumentException ("term overflow");
      return Math.exp (y);
   }

   /**
    * Computes the distribution function.
    */
   public static double cdf (double n, double p, int x) {
      final double EPSILON = DiscreteDistributionInt.EPSILON;
      final int LIM1 = 100000;
      double sum, term, termmode;
      int i, mode;
      final double q = 1.0 - p;

      if (p < 0.0 || p > 1.0)
        throw new IllegalArgumentException ("p not in [0, 1]");
      if (n <= 0.0)
        throw new IllegalArgumentException ("n <= 0.0");

      if (x < 0)
         return 0.0;
      if (p >= 1.0)                  // In fact, p == 1
         return 1.0;
      if (p <= 0.0)                  // In fact, p == 0
         return 0.0;

      // Compute the maximum term
      mode = 1 + (int) Math.floor ((n*q - 1.0)/p);
      if (mode < 0)
          mode = 0;
      else if (mode > x)
         mode = x;

      if (mode <= LIM1) {
         sum = term = termmode = prob (n, p, mode);
         for (i = mode; i > 0; i--) {
            term *= i/(q*(n + i - 1.0));
            if (term < EPSILON)
               break;
            sum += term;
         }

         term = termmode;
         for (i = mode; i < x; i++) {
            term *= q*(n + i)/(i + 1);
            if (term < EPSILON)
               break;
            sum += term;
         }
         if (sum <= 1.0)
            return sum;
         else
            return 1.0;
      }
      else
         //return 1.0 - BinomialDist.cdf (x + n, p, n - 1);
         return BetaDist.cdf (n, x + 1.0, p);
    }

   /**
    * Returns @f$\bar{F}(x) = P[X \ge x]@f$, the complementary
    * distribution function.
    */
   public static double barF (double n, double p, int x) {
      return 1.0 - cdf (n, p, x - 1);
   }

   /**
    * Computes the inverse function without precomputing tables.
    */
   public static int inverseF (double n, double p, double u) {
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u is not in [0,1]");
      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in [0, 1]");
      if (n <= 0.0)
         throw new IllegalArgumentException ("n <= 0");
      if (p >= 1.0)                  // In fact, p == 1
         return 0;
      if (p <= 0.0)                  // In fact, p == 0
         return 0;
      if (u <= prob (n, p, 0))
         return 0;
      if (u >= 1.0)
         return Integer.MAX_VALUE;

      double sum, term, termmode;
      final double q = 1.0 - p;

      // Compute the maximum term
      int mode = 1 + (int) Math.floor ((n * q - 1.0) / p);
      if (mode < 0)
         mode = 0;
      int i = mode;
      term = prob (n, p, i);
      while ((term >= u) && (term > Double.MIN_NORMAL)) {
         i /= 2;
         term = prob (n, p, i);
      }

      if (term <= Double.MIN_NORMAL) {
         i *= 2;
         term = prob (n, p, i);
         while (term >= u && (term > Double.MIN_NORMAL)) {
            term *= i / (q * (n + i - 1.0));
            i--;
         }
      }

      mode = i;
      sum = termmode = prob (n, p, i);

      for (i = mode; i > 0; i--) {
         term *= i / (q * (n + i - 1.0));
         if (term < EPSILON)
            break;
         sum += term;
      }

      term = termmode;
      i = mode;
      double prev = -1;
      if (sum < u) {
         // The CDF at the mode is less than u, so we add term to get >= u.
         while ((sum < u) && (sum > prev)){
            term *= q * (n + i) / (i + 1);
            prev = sum;
            sum += term;
            i++;
         }
      } else {
         // The computed CDF is too big so we substract from it.
         sum -= term;
         while (sum >= u) {
            term *= i / (q * (n + i - 1.0));
            i--;
            sum -= term;
         }
      }
      return i;
   }

   /**
    * Estimates the parameter @f$p@f$ of the negative binomial
    * distribution using the maximum likelihood method, from the @f$m@f$
    * observations @f$x[i]@f$, @f$i = 0, 1, …, m-1@f$. The parameter
    * @f$n@f$ is assumed known. The estimate @f$\hat{p}@f$ is returned in
    * element 0 of the returned array. The maximum likelihood estimator
    * @f$\hat{p}@f$ satisfies the equation @f$\hat{p} = n /(n +
    * \bar{x}_m)@f$, where @f$\bar{x}_m@f$ is the average of @f$x[0], …,
    * x[m-1]@f$.
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param m            the number of observations used to evaluate
    *                      parameters
    *  @param n            the first parameter of the negative binomial
    *  @return returns the parameters [@f$\hat{p}@f$]
    */
   public static double[] getMLE (int[] x, int m, double n) {
      if (m <= 0)
         throw new IllegalArgumentException ("m <= 0");
      double mean = 0.0;
      for (int i = 0; i < m; i++) {
         mean += x[i];
      }
      mean /= (double) m;
      double[] param = new double[1];
      param[0] = n / (n + mean);
      return param;
   }

   /**
    * Creates a new instance of a negative binomial distribution with
    * parameters @f$n@f$ given and @f$\hat{p}@f$ estimated using the
    * maximum likelihood method, from the @f$m@f$ observations @f$x[i]@f$,
    * @f$i = 0, 1, …, m-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param m            the number of observations to use to evaluate
    *                      parameters
    *  @param n            the first parameter of the negative binomial
    */
   public static NegativeBinomialDist getInstanceFromMLE (int[] x, int m,
                                                          double n) {
      double parameters[] = getMLE (x, m, n);
      return new NegativeBinomialDist (n, parameters[0]);
   }

   /**
    * Estimates the parameter @f$n@f$ of the negative binomial
    * distribution using the maximum likelihood method, from the @f$m@f$
    * observations @f$x[i]@f$, @f$i = 0, 1, …, m-1@f$. The parameter
    * @f$p@f$ is assumed known. The estimate @f$\hat{n}@f$ is returned in
    * element 0 of the returned array.  The maximum likelihood estimator
    * @f$\hat{p}@f$ satisfies the equation
    * @f[
    *   \frac{1}{m}\sum_{j=0}^{m-1} \psi(n +x_j) = \psi(n) - \ln(p)
    * @f]
    * where @f$\psi(x)@f$ is the digamma function, i.e. the logarithmic
    * derivative of the Gamma function @f$\psi(x) =
    * \Gamma^{\prime}(x)/\Gamma(x)@f$.
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param m            the number of observations used to evaluate
    *                      parameters
    *  @param p            the second parameter of the negative binomial
    *  @return returns the parameters [@f$\hat{n}@f$]
    */
   public static double[] getMLE1 (int[] x, int m, double p) {
      if (m <= 0)
         throw new IllegalArgumentException ("m <= 0");
      double mean = 0.0;
      for (int i = 0; i < m; i++)
         mean += x[i];
      mean /= m;

      double gam0 = mean*p/(1.0 - p);
      double[] param = new double[1];
      Func1 f = new Func1 (p, x, m);
      param[0] = RootFinder.brentDekker (gam0/100.0, 100.0*gam0, f, 1e-5);
      return param;
   }

   /**
    * Creates a new instance of a negative binomial distribution with
    * parameters @f$p@f$ given and @f$\hat{n}@f$ estimated using the
    * maximum likelihood method, from the @f$m@f$ observations @f$x[i]@f$,
    * @f$i = 0, 1, …, m-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param m            the number of observations to use to evaluate
    *                      parameters
    *  @param p            the second parameter of the negative binomial
    */
   public static NegativeBinomialDist getInstanceFromMLE1 (int[] x, int m,
                                                           double p) {
      double param[] = getMLE1 (x, m, p);
      return new NegativeBinomialDist (param[0], p);
   }

   /**
    * Estimates the parameter @f$(n, p)@f$ of the negative binomial
    * distribution using the maximum likelihood method, from the @f$m@f$
    * observations @f$x[i]@f$, @f$i = 0, 1, …, m-1@f$. The estimates are
    * returned in a two-element array, in regular order: [@f$n@f$,
    * @f$p@f$].  The maximum likelihood estimators are the values
    * @f$(\hat{n}@f$, @f$\hat{p})@f$ satisfying the equations
    * @f{align*}{
    *    \frac{\hat{n}(1 - \hat{p})}{\hat{p}} 
    *    & 
    *    = 
    *    \bar{x}_m
    *    \\ 
    *   \sum_{j=1}^{\infty} \frac{F_j}{(\hat{n} + j - 1)} 
    *    & 
    *    = 
    *    -m\ln(\hat{p})
    * @f}
    * where @f$\bar{x}_m@f$ is the average of @f$x[0],…,x[m-1]@f$, and
    * @f$F_j = \sum_{i=j}^{\infty} f_i@f$ = number of @f$x_i \ge j@f$
    * (see @cite tJOH69a&thinsp; (page 132)).
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param m            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{n}@f$, @f$\hat{p}@f$]
    */
   public static double[] getMLE (int[] x, int m) {
      if (m <= 0)
         throw new IllegalArgumentException ("m<= 0");

      int i, j;
      int max = Integer.MIN_VALUE;
      double sum = 0.0;
      for (i = 0; i < m; i++) {
         sum += x[i];
         if (x[i] > max)
            max = x[i];
      }
      double mean = sum / (double) m;

      double var = 0.0;
      for (i = 0; i < m; i++)
         var += (x[i] - mean) * (x[i] - mean);
      var /= (double) m;

      if (mean >= var) {
         throw new UnsupportedOperationException("mean >= variance");
      }
      double estimGamma = (mean * mean) / ( var - mean );

      int[] Fj = new int[max];
      for (j = 0; j < max; j++) {
         int prop = 0;
         for (i = 0; i < m; i++)
            if (x[i] > j)
               prop++;

         Fj[j] = prop;
      }

      double[] param = new double[3];
      Function f = new Function (m, max, mean, Fj);
      param[1] = RootFinder.brentDekker (estimGamma/100, estimGamma*100, f, 1.0e-5);
      param[2] = param[1] / (param[1] + mean);

/*    // Seems to be useless
      Optim system = new Optim (m, max, mean, Fj);
      double[] fvec = new double [3];
      double[][] fjac = new double[3][3];
      int[] iflag = new int[2];
      int[] info = new int[2];
      int[] ipvt = new int[3];

      Minpack_f77.lmder1_f77 (system, 2, 2, param, fvec, fjac, 1e-5, info, ipvt);
*/
      double parameters[] = new double[2];
      parameters[0] = param[1];
      parameters[1] = param[2];

      return parameters;
   }

   /**
    * Creates a new instance of a negative binomial distribution with
    * parameters @f$n@f$ and @f$p@f$ estimated using the maximum
    * likelihood method based on the @f$m@f$ observations @f$x[i]@f$, @f$i
    * = 0, 1, …, m-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param m            the number of observations used to evaluate
    *                      parameters
    */
   public static NegativeBinomialDist getInstanceFromMLE (int[] x, int m) {
      double parameters[] = getMLE (x, m);
      return new NegativeBinomialDist (parameters[0], parameters[1]);
   }

   /**
    * Estimates and returns the parameter @f$\nu= 1/\hat{n}@f$ of the
    * negative binomial distribution using the maximum likelihood method,
    * from the @f$m@f$ observations @f$x[i]@f$, @f$i = 0, 1, …, m-1@f$.
    * The maximum likelihood estimator is the value @f$\nu@f$ satisfying
    * the equation
    * @f[
    *   \sum_{j=1}^{\infty} \frac{\nu F_j}{1 + \nu(j - 1)} = m\ln(1 + \nu\bar{x}_m)
    * @f]
    * where @f$\bar{x}_m@f$ is the average of @f$x[0],…,x[m-1]@f$, and
    * @f$F_j = \sum_{i=j}^{\infty} f_i@f$ = number of @f$x_i \ge j@f$
    * (see @cite tJOH69a&thinsp; (page 132)).
    *  @param x            the list of observations used to evaluate
    *                      parameter
    *  @param m            the number of observations used to evaluate
    *                      parameter
    *  @return returns the parameter @f$\nu@f$
    */
   public static double getMLEninv (int[] x, int m) {
      if (m <= 0)
         throw new IllegalArgumentException ("m <= 0");

      int i, j;
      int max = Integer.MIN_VALUE;
      double mean = 0.0;
      for (i = 0; i < m; i++) {
         mean += x[i];
         if (x[i] > max)
            max = x[i];
      }
      mean /= (double) m;

      double var = 0.0;
      for (i = 0; i < m; i++)
         var += (x[i] - mean) * (x[i] - mean);
      var /= (double) m;

      if (mean >= var) {
         throw new UnsupportedOperationException("mean >= variance");
      }

      int[] Fj = new int[max];
      for (j = 0; j < max; j++) {
         int prop = 0;
         for (i = 0; i < m; i++)
            if (x[i] > j)
               prop++;

         Fj[j] = prop;
      }

      FuncInv f = new FuncInv (m, max, mean, Fj);
      double nu = RootFinder.brentDekker (1.0e-8, 1.0e8, f, 1.0e-5);
      return nu;
   }

   /**
    * Computes and returns the mean @f$E[X] = n(1 - p)/p@f$ of the
    * negative binomial distribution with parameters @f$n@f$ and @f$p@f$.
    *  @return the mean of the negative binomial distribution @f$E[X] =
    * n(1 - p) / p@f$
    */
   public static double getMean (double n, double p) {
      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in [0, 1]");
      if (n <= 0.0)
         throw new IllegalArgumentException ("n <= 0");

      return (n * (1 - p) / p);
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] = n(1 - p)/p^2@f$
    * of the negative binomial distribution with parameters @f$n@f$ and
    * @f$p@f$.
    *  @return the variance of the negative binomial distribution
    * @f$\mbox{Var}[X] = n(1 - p) / p^2@f$
    */
   public static double getVariance (double n, double p) {
      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in [0, 1]");
      if (n <= 0.0)
         throw new IllegalArgumentException ("n <= 0");

      return (n * (1 - p) / (p * p));
   }

   /**
    * Computes and returns the standard deviation of the negative binomial
    * distribution with parameters @f$n@f$ and @f$p@f$.
    *  @return the standard deviation of the negative binomial
    * distribution
    */
   public static double getStandardDeviation (double n, double p) {
      return Math.sqrt (NegativeBinomialDist.getVariance (n, p));
   }

   /**
    * Returns the parameter @f$n@f$ of this object.
    */
   @Deprecated
   public double getGamma() {
      return n;
   }

   /**
    * Returns the parameter @f$n@f$ of this object.
    */
   public double getN() {
      return n;
   }

   /**
    * Returns the parameter @f$p@f$ of this object.
    */
   public double getP() {
      return p;
   }

   /**
    * Sets the parameter @f$n@f$ and @f$p@f$ of this object.
    */
   public void setParams (double n, double p) {
      /**
      *  Compute all probability terms of the negative binomial distribution;
      *  start at the mode, and calculate probabilities on each side until they
      *  become smaller than EPSILON. Set all others to 0.
      */
      supportA = 0;
      int i, mode, Nmax;
      int imin, imax;
      double sum;
      double[] P;     // Negative Binomial mass probabilities
      double[] F;     // Negative Binomial cumulative

      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in [0, 1]");
      if (n <= 0.0)
         throw new IllegalArgumentException ("n <= 0");

      this.n  = n;
      this.p  = p;

      // Compute the mode (at the maximum term)
      mode = 1 + (int) Math.floor((n*(1.0 - p) - 1.0)/p);

      /**
       For mode > MAXN, we shall not use pre-computed arrays.
       mode < 0 should be impossible, unless overflow of long occur, in
       which case mode will be = LONG_MIN.
      */

      if (mode < 0.0 || mode > MAXN) {
         pdf = null;
         cdf = null;
         return;
      }

      /**
        In theory, the negative binomial distribution has an infinite range.
        But for i > Nmax, probabilities should be extremely small.
        Nmax = Mean + 16 * Standard deviation.
      */

      Nmax = (int)(n*(1.0 - p)/p + 16*Math.sqrt (n*(1.0 - p)/(p*p)));
      if (Nmax < 32)
         Nmax = 32;
      P = new double[1 + Nmax];

      double epsilon = EPSILON/prob (n, p, mode);

      // We shall normalize by explicitly summing all terms >= epsilon
      sum = P[mode] = 1.0;

      // Start from the maximum and compute terms > epsilon on each side.
      i = mode;
      while (i > 0 && P[i] >= epsilon) {
         P[i - 1] = P[i]*i/((1.0 - p)*(n + i - 1));
         i--;
         sum += P[i];
      }
      imin = i;

      i = mode;
      while (P[i] >= epsilon) {
         P[i + 1] = P[i]*(1.0 - p)*(n + i)/(i + 1);
         i++;
         sum += P[i];
         if (i == Nmax - 1) {
            Nmax *= 2;
            double[] nT = new double[1 + Nmax];
            System.arraycopy (P, 0, nT, 0, P.length);
            P = nT;
         }
      }
      imax = i;

      // Renormalize the sum of probabilities to 1
      for (i = imin; i <= imax; i++)
         P[i] /= sum;

      // Compute the cumulative probabilities for F and keep them in the
      // lower part of CDF.
      F = new double[1 + Nmax];
      F[imin] = P[imin];
      i = imin;
      while (i < imax && F[i] < 0.5) {
         i++;
         F[i] = F[i - 1] + P[i];
      }

      // This is the boundary between F (i <= xmed) and 1 - F (i > xmed) in
      // the array CDF
      xmed = i;

      // Compute the cumulative probabilities of the complementary
      // distribution 1 - F and keep them in the upper part of the array
      F[imax] = P[imax];
      i = imax - 1;
      do {
         F[i] = P[i] + F[i + 1];
         i--;
      } while (i > xmed);

     xmin = imin;
     xmax = imax;
     pdf = new double[imax + 1 - imin];
     cdf = new double[imax + 1 - imin];
     System.arraycopy (P, imin, pdf, 0, imax + 1 - imin);
     System.arraycopy (F, imin, cdf, 0, imax + 1 - imin);

   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$n@f$,
    * @f$p@f$].
    */
   public double[] getParams () {
      double[] retour = {n, p};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : n = " + n + ", p = " + p;
   }

}