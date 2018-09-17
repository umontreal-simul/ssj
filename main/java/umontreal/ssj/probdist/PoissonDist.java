/*
 * Class:        PoissonDist
 * Description:  Poisson distribution
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

/**
 * Extends the class  @ref DiscreteDistributionInt for the *Poisson*
 * distribution @cite sLAW00a&thinsp; (page 325) with mean
 * @f$\lambda\ge0@f$. The mass function is
 * @anchor REF_probdist_PoissonDist_eq_fmass_Poisson
 * @f{align}{
 *    p(x) 
 *    & 
 *    =
 *    \frac{e^{-\lambda} \lambda^x}{x!}, \qquad\mbox{for } x=0,1,…\tag{fmass-Poisson}
 * @f}
 * and the distribution function is
 * @anchor REF_probdist_PoissonDist_eq_FPoisson
 * @f{align}{
 *    F(x) 
 *    & 
 *    =
 *    e^{-\lambda} \sum_{j=0}^x\; \frac{\lambda^j}{j!}, \qquad\mbox{for } x=0,1,…. \tag{FPoisson}
 * @f}
 * If one has to compute @f$p(x)@f$ and/or @f$F(x)@f$ for several values of
 * @f$x@f$ with the same @f$\lambda@f$, where @f$\lambda@f$ is not too
 * large, then it is more efficient to instantiate an object and use the
 * non-static methods, since the functions will then be computed once and
 * kept in arrays.
 *
 * For the static methods that compute @f$F(x)@f$ and @f$\bar{F}(x)@f$, we
 * exploit the relationship @f$F(x) = 1 - G_{x+1}(\lambda)@f$, where
 * @f$G_{x+1}@f$ is the *gamma* distribution function with parameters
 * @f$(\alpha,\lambda) = (x+1, 1)@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_discrete
 */
public class PoissonDist extends DiscreteDistributionInt {

   private double lambda;

   /**
    * @name Constant
    * @{
    */

   /**
    * The value of the parameter @f$\lambda@f$ above which the tables are
    * *not* precomputed by the constructor.
    */
   public static double MAXLAMBDA = 100000;

   /**
    * @}
    */

   /**
    * Creates an object that contains the probability and distribution
    * functions, for the Poisson distribution with parameter `lambda`,
    * which are computed and stored in dynamic arrays inside that object.
    */
   public PoissonDist (double lambda) {
      setLambda (lambda);
   }


   public double prob (int x) {
      if (x < 0)
         return 0.0;
      if (pdf == null)
         return prob (lambda, x);
      if (x > xmax || x < xmin)
         return prob (lambda, x);
      return pdf[x - xmin];
   }

   public double cdf (int x) {
      double Sum = 0.0;
      // int j;

      if (x < 0)
         return 0.0;
      if (lambda == 0.0)
         return 1.0;

      /* For large lambda, we use the Chi2 distribution according to the exact
         relation, with 2x + 2 degrees of freedom

         cdf (lambda, x) = 1 - chiSquare (2x + 2, 2*lambda)

         which equals also 1 - gamma (x + 1, lambda) */
      if (cdf == null)
         return GammaDist.barF (x + 1.0, 15, lambda);

      if (x >= xmax)
         return 1.0;

      if (x < xmin) {
         // Sum a few terms to get a few decimals far in the lower tail. One
         // could also call GammaDist.barF instead.
         final int RMAX = 20;
         int i;
         double term = prob(lambda, x);
         Sum = term;
         i = x;
         while (i > 0 && i >= x - RMAX) {
            term = term * i / lambda;
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
   }


   public double barF (int x) {
      /*
       * poisson (lambda, x) = 1 - cdf (lambda, x - 1)
       */

      if (x <= 0)
         return 1.0;

      /* For large lambda,  we use the Chi2 distribution according to the exact
         relation, with 2x + 2 degrees of freedom

         cdf (lambda, x) = 1 - ChiSquare.cdf (2x + 2, 2*lambda)
         cdf (lambda, x) = 1 - GammaDist.cdf (x + 1, lambda)
       */

      if (cdf == null)
         return GammaDist.cdf ((double)x, 15, lambda);

      if (x > xmax)
//         return GammaDist.cdf ((double)x, 15, lambda);
         return PoissonDist.barF(lambda, x);
      if (x <= xmin)
         return 1.0;
      if (x > xmed)
         // We keep the complementary distribution in the upper part of cdf
         return cdf[x - xmin];
      else
         return 1.0 - cdf[x - 1 - xmin];
   }


   public int inverseFInt (double u) {
      if ((cdf == null) || (u <= EPSILON))
         return inverseF (lambda, u);
      return super.inverseFInt (u);
   }

   public double getMean() {
      return PoissonDist.getMean (lambda);
   }

   public double getVariance() {
      return PoissonDist.getVariance (lambda);
   }

   public double getStandardDeviation() {
      return PoissonDist.getStandardDeviation (lambda);
   }

/**
 * Computes and returns the Poisson probability @f$p(x)@f$ for @f$\lambda=
 * @f$ `lambda`, as defined in (
 * {@link REF_probdist_PoissonDist_eq_fmass_Poisson
 * fmass-Poisson} ).  If @f$\lambda\ge20@f$, this (static) method uses the
 * logarithm of the gamma function, defined in (
 * {@link REF_probdist_GammaDist_eq_Gamma Gamma} ), to estimate
 * the density.
 */
public static double prob (double lambda, int x) {
      if (x < 0)
         return 0.0;

      if (lambda >= 100.0) {
         if ((double) x >= 10.0*lambda)
            return 0.0;
      } else if (lambda >= 3.0) {
         if ((double) x >= 100.0*lambda)
            return 0.0;
      } else {
         if ((double) x >= 200.0*Math.max(1.0, lambda))
            return 0.0;
      }

      final double LAMBDALIM = 20.0;
      double Res;
      if (lambda < LAMBDALIM && x <= 100)
         Res = Math.exp (-lambda)*Math.pow (lambda, x)/Num.factorial (x);
      else {
         double y = x*Math.log (lambda) - Num.lnGamma (x + 1.0) - lambda;
         Res = Math.exp (y);
      }
      return Res;
   }

   /**
    * Computes and returns the value of the Poisson distribution function
    * @f$F(x)@f$ for @f$\lambda= @f$ `lambda`, as defined in (
    * {@link REF_probdist_PoissonDist_eq_FPoisson FPoisson}
    * ).  To compute @f$F(x)@f$, all non-negligible terms of the sum are
    * added if @f$\lambda\le200@f$; otherwise, the relationship
    * @f$F_{\lambda}(x) = 1 - G_{x + 1}(\lambda)@f$ is used, where
    * @f$G_{x+1}@f$ is the gamma distribution function with parameter
    * @f$\alpha= x+1@f$ (see  @ref GammaDist ).
    */
   public static double cdf (double lambda, int x) {
   /*
    * On our machine, computing a value using gamma is faster than the
    * naive computation for lambdalim > 200.0, slower for lambdalim < 200.0
    */
      if (lambda < 0.0)
        throw new IllegalArgumentException ("lambda < 0");
      if (lambda == 0.0)
         return 1.0;
      if (x < 0)
         return 0.0;

      if (lambda >= 100.0) {
         if ((double) x >= 10.0*lambda)
            return 1.0;
      } else {
         if ((double) x >= 100.0*Math.max(1.0, lambda))
            return 1.0;
      }

      /* If lambda > LAMBDALIM, use the Chi2 distribution according to the
         exact relation, with 2x + 2 degrees of freedom

         poisson (lambda, x) = 1 - chiSquare (2x + 2, 2*lambda)

         which also equals 1 - gamma (x + 1, lambda) */

      final double LAMBDALIM = 200.0;
      if (lambda > LAMBDALIM)
         return GammaDist.barF (x + 1.0, 15, lambda);

      if (x >= lambda)
         return 1 - PoissonDist.barF(lambda, x+1);

      // Naive computation: sum all prob. from i = x
      double sum = 1;
      double term = 1;
      for(int j = 1; j <= x; j++) {
         term *= lambda/j;
         sum += term;
      }
      return sum*Math.exp(-lambda);
   }

   /**
    * Computes and returns the value of the complementary Poisson
    * distribution function, for @f$\lambda= @f$ `lambda`. *WARNING:* The
    * complementary distribution function is defined as @f$\bar{F}(x) =
    * P[X \ge x]@f$.
    */
   public static double barF (double lambda, int x) {
      if (lambda < 0)
         throw new IllegalArgumentException ("lambda < 0");
      if (x <= 0)
         return 1.0;

      if (lambda >= 100.0) {
         if ((double) x >= 10.0*lambda)
            return 0.0;
      } else {
         if ((double) x >= 100 + 100.0*Math.max(1.0, lambda))
            return 0.0;
      }

      /* If lambda > LAMBDALIM, we use the Chi2 distribution according to the
         exact relation, with 2x + 2 degrees of freedom

         cdf (lambda, x) = 1 - ChiSquare.cdf (2x + 2, 2*lambda)

         which also equals   1 - GammaDist.cdf (x + 1, lambda) */

      final double LAMBDALIM = 200.0;
      if (lambda > LAMBDALIM)
         return GammaDist.cdf ((double)x, 15, lambda);

      if (x <= lambda)
         return 1.0 - PoissonDist.cdf(lambda, x - 1);

      // Naive computation: sum all prob. from i = x to i = oo
      double term, sum;
      final int IMAX = 20;

      // Sum at least IMAX prob. terms from i = s to i = oo
      sum = term = PoissonDist.prob(lambda, x);
      int i = x + 1;
      while (term > EPSILON || i <= x + IMAX) {
         term *= lambda/i;
         sum += term;
         i++;
      }
      return sum;
   }

   /**
    * Performs a linear search to get the inverse function without
    * precomputed tables.
    */
   public static int inverseF (double lambda, double u) {
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u is not in range [0,1]");
      if (lambda < 0.0)
         throw new IllegalArgumentException ("lambda < 0");
      if (u >= 1.0)
         return Integer.MAX_VALUE;
      if (u <= prob (lambda, 0))
         return 0;
      int i;

      final double LAMBDALIM = 700.0;
      if (lambda < LAMBDALIM) {
         double sumprev = -1.0;
         double term = Math.exp(-lambda);
         double sum = term;
         i = 0;
         while (sum < u && sum > sumprev) {
            i++;
            term *= lambda / i;
            sumprev = sum;
            sum += term;
         }
         return i;

      } else {
         i = (int)lambda;
         double term = PoissonDist.prob(lambda, i);
         while ((term >= u) && (term > Double.MIN_NORMAL)) {
            i /= 2;
            term = PoissonDist.prob (lambda, i);
         }
         if (term <= Double.MIN_NORMAL) {
            i *= 2;
            term = PoissonDist.prob (lambda, i);
            while (term >= u && (term > Double.MIN_NORMAL)) {
               term *= i / lambda;
               i--;
            }
         }
         int mid = i;
         double sum = term;
         double termid = term;

         while (term >= EPSILON*u && i > 0) {
            term *= i / lambda;
            sum += term;
            i--;
         }

         term = termid;
         i = mid;
        double prev = -1;
        if (sum < u) {
            while ((sum < u) && (sum > prev)) {
               i++;
               term *= lambda / i;
               prev = sum;
               sum += term;
            }
         } else {
            // The computed CDF is too big so we substract from it.
            sum -= term;
            while (sum >= u) {
               term *= i / lambda;
               i--;
               sum -= term;
            }
         }
      }
      return i;
   }

   /**
    * Estimates the parameter @f$\lambda@f$ of the Poisson distribution
    * using the maximum likelihood method, from the @f$n@f$ observations
    * @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$. The maximum likelihood estimator
    * @f$\hat{\lambda}@f$ satisfy the equation @f$\hat{\lambda} =
    * \bar{x}_n@f$, where @f$\bar{x}_n@f$ is the average of @f$x[0], …,
    * x[n-1]@f$ (see @cite sLAW00a&thinsp; (page 326)).
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameter [@f$\hat{\lambda}@f$]
    */
   public static double[] getMLE (int[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double parameters[];
      parameters = new double[1];
      double sum = 0.0;
      for (int i = 0; i < n; i++) {
         sum += x[i];
      }

      parameters[0] = (double) sum / (double) n;
      return parameters;
   }

   /**
    * Creates a new instance of a Poisson distribution with parameter
    * @f$\lambda@f$ estimated using the maximum likelihood method based
    * on the @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static PoissonDist getInstanceFromMLE (int[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new PoissonDist (parameters[0]);
   }

   /**
    * Computes and returns the mean @f$E[X] = \lambda@f$ of the Poisson
    * distribution with parameter @f$\lambda@f$.
    *  @return the mean of the Poisson distribution @f$E[X] = \lambda@f$
    */
   public static double getMean (double lambda) {
      if (lambda < 0.0)
       throw new IllegalArgumentException ("lambda < 0");

      return lambda;
   }

   /**
    * Computes and returns the variance @f$= \lambda@f$ of the Poisson
    * distribution with parameter @f$\lambda@f$.
    *  @return the variance of the Poisson distribution @f$\mbox{Var}[X] =
    * \lambda@f$
    */
   public static double getVariance (double lambda) {
      if (lambda < 0.0)
       throw new IllegalArgumentException ("lambda < 0");

      return lambda;
   }

   /**
    * Computes and returns the standard deviation of the Poisson
    * distribution with parameter @f$\lambda@f$.
    *  @return the standard deviation of the Poisson distribution
    */
   public static double getStandardDeviation (double lambda) {
      if (lambda < 0.0)
       throw new IllegalArgumentException ("lambda < 0");

      return Math.sqrt (lambda);
   }

   /**
    * Returns the @f$\lambda@f$ associated with this object.
    */
   public double getLambda() {
      return lambda;
   }

   /**
    * Sets the @f$\lambda@f$ associated with this object.
    */
   public void setLambda (double lambda) {
      supportA = 0;
      if (lambda < 0.0)
         throw new IllegalArgumentException ("lambda < 0");
      this.lambda = lambda;

      // For lambda > MAXLAMBDAPOISSON, we do not use pre-computed arrays
      if (lambda > MAXLAMBDA) {
         pdf = null;
         cdf = null;
         return;
      }

      double epsilon;
      int i, mid, Nmax;
      int imin, imax;
      double sum;
      double[] P;    // Poisson probability terms
      double[] F;    // Poisson cumulative probabilities

      // In theory, the Poisson distribution has an infinite range. But
      // for i > Nmax, probabilities should be extremely small.
      Nmax = (int)(lambda + 16*(2 + Math.sqrt (lambda)));
      P = new double[1 + Nmax];

      mid = (int)lambda;
      epsilon = EPSILON * EPS_EXTRA/prob (lambda, mid);
      // For large lambda, mass will lose a few digits of precision
      // We shall normalize by explicitly summing all terms >= epsilon
      sum = P[mid] = 1.0;

      // Start from the maximum and compute terms > epsilon on each side.
      i = mid;
      while (i > 0 && P[i] > epsilon) {
         P[i - 1] = P[i]*i/lambda;
         i--;
         sum += P[i];
      }
      xmin = imin = i;

      i = mid;
      while (P[i] > epsilon) {
         P[i + 1] = P[i]*lambda/(i + 1);
         i++;
         sum += P[i];
         if (i >= Nmax - 1) {
            Nmax *= 2;
            double[] nT = new double[1 + Nmax];
            System.arraycopy (P, 0, nT, 0, P.length);
            P = nT;
         }
      }
      xmax = imax = i;
      F = new double[1 + Nmax];

      // Renormalize the sum of probabilities to 1
      for (i = imin; i <= imax; i++)
         P[i] /= sum;

      // Compute the cumulative probabilities until F >= 0.5, and keep them in
      // the lower part of array, i.e. F[s] contains all P[i] for i <= s
      F[imin] = P[imin];
      i = imin;
      while (i < imax && F[i] < 0.5) {
         i++;
         F[i] = P[i] + F[i - 1];
      }
      // This is the boundary between F and 1 - F in the CDF
      xmed = i;

      // Compute the cumulative probabilities of the complementary distribution
      // and keep them in the upper part of the array. i.e. F[s] contains all
      // P[i] for i >= s
      F[imax] = P[imax];
      i = imax - 1;
      do {
         F[i] = P[i] + F[i + 1];
         i--;
      } while (i > xmed);

       /* Reset imin because we lose too much precision for a few terms near
      imin when we stop adding terms < epsilon. */
      i = imin;
      while (i < xmed && F[i] < EPSILON)
         i++;
      xmin = imin = i;

      /* Same thing with imax */
      i = imax;
      while (i > xmed && F[i] < EPSILON)
         i--;
      xmax = imax = i;

      pdf = new double[imax + 1 - imin];
      cdf = new double[imax + 1 - imin];
      System.arraycopy (P, imin, pdf, 0, imax-imin+1);
      System.arraycopy (F, imin, cdf, 0, imax-imin+1);
   }

   /**
    * Return a table containing the parameter of the current distribution.
    */
   public double[] getParams () {
      double[] retour = {lambda};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + ": lambda = " + lambda;
   }

}