/*
 * Class:        HypoExponentialDist
 * Description:  Hypo-exponential distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        January 2011
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
import java.util.Formatter;
import java.util.Locale;
import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;
import umontreal.ssj.util.*;
import umontreal.ssj.functions.MathFunction;

/**
 * This class implements the *hypoexponential* distribution, also called the
 * generalized Erlang distribution. Let the @f$X_j@f$, @f$j=1,…,k@f$, be
 * @f$k@f$ independent exponential random variables with different rates
 * @f$\lambda_j@f$, i.e. assume that @f$\lambda_j \neq\lambda_i@f$ for
 * @f$i \neq j@f$. Then the sum @f$\sum_{j=1}^kX_j@f$ is called a
 * *hypoexponential* random variable.
 *
 * Let the @f$k\times k@f$ upper triangular bidiagonal matrix
 * @anchor REF_probdist_HypoExponentialDist_eq_tail_hypomatrix
 * @f[
 *   \tag{tail-hypomatrix} \mathbf{A}= \begin{pmatrix}
 *    -\lambda_1 
 *    & 
 *    \lambda_1 
 *    & 
 *    0 
 *    & 
 *    …
 *    & 
 *    0 
 *    \\ 
 *   0 
 *    & 
 *    -\lambda_2
 *    & 
 *    \lambda_2 
 *    & 
 *    …
 *    & 
 *    0 
 *    \\ 
 *   \vdots
 *    & 
 *    \vdots
 *    & 
 *    \ddots
 *    & 
 *    &nbsp;&nbsp;\ddots
 *    & 
 *    \vdots
 *    \\ 
 *   0 
 *    & 
 *    …
 *    & 
 *    0 
 *    & 
 *    -\lambda_{k-1} 
 *    & 
 *    \lambda_{k-1} 
 *    \\ 
 *   0 
 *    & 
 *    …
 *    & 
 *    0 
 *    & 
 *    0 
 *    & 
 *    -\lambda_k 
 *   \end{pmatrix}
 * @f]
 * with @f$\lambda_j@f$ the rates of the @f$k@f$ exponential random
 * variables; then the cumulative complementary probability of the
 * hypoexponential distribution is given by @cite pNEU81a,
 * @cite pLAT99a&thinsp;
 * @anchor REF_probdist_HypoExponentialDist_eq_tail_hypoexp
 * @f[
 *   \tag{tail-hypoexp} \bar{F}(x) = 
 *    \mathbb P 
 *     \left[X_1 + \cdots+ X_k > x \right] = \sum_{j=1}^k \left(e^{\mathbf{A}x}\right)_{1j},
 * @f]
 * i.e., it is the sum of the elements of the first row of matrix
 * @f$e^{\mathbf{A}x}@f$. The density of the hypoexponential distribution is
 * @anchor REF_probdist_HypoExponentialDist_eq_fhypoexp
 * @f[
 *   f(x) = \left(-e^{\mathbf{A}x}\mathbf{A}\right)_{1k} = \lambda_k \left(e^{\mathbf{A}x}\right)_{1k}, \tag{fhypoexp}
 * @f]
 * i.e., it is element @f$(1,k)@f$ of matrix
 * @f$-e^{\mathbf{A}x}\mathbf{A}@f$. The distribution function is as usual
 * @f$F(x) = 1 - \bar{F}(x)@f$.
 *
 * See the class  @ref HypoExponentialDistQuick for alternative formulae for
 * the probabilities.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class HypoExponentialDist extends ContinuousDistribution {
   protected double[] m_lambda;

   protected static void testLambda (double[] lambda) {
      int m = lambda.length;
      for (int j = 0; j < m; ++j) {
         if (lambda[j] <= 0)
            throw new IllegalArgumentException ("lambda_j <= 0");
      }
   }


   // Builds the bidiagonal matrix A out of the lambda
   private static DoubleMatrix2D buildMatrix (double[] lambda, double x) {
      int m = lambda.length;
      testLambda (lambda);
      DoubleFactory2D F2 = DoubleFactory2D.dense;
      DoubleMatrix2D A = F2.make(m, m);
      for (int j = 0; j < m-1; j++) {
         A.setQuick(j, j, -lambda[j]*x);
         A.setQuick(j, j + 1, lambda[j]*x);
      }
      A.setQuick(m-1, m-1, -lambda[m-1]*x);
      return A;
   }


   private static class myFunc implements MathFunction {
      // For inverseF
      private double[] m_lam;
      private double m_u;

      public myFunc (double[] lam, double u) {
         m_lam = lam;
         m_u = u;
      }

      public double evaluate (double x) {
         return m_u - HypoExponentialDist.cdf(m_lam, x);
      }
   }

   /**
    * Constructs a `HypoExponentialDist` object, with rates @f$\lambda_i
    * = @f$ <tt>lambda[</tt>@f$i-1@f$<tt>]</tt>, @f$i = 1,…,k@f$.
    *  @param lambda       rates of the hypoexponential distribution
    */
   public HypoExponentialDist (double[] lambda) {
      supportA = 0.0;
      setLambda (lambda);
  }


   public double density (double x) {
      return density (m_lambda, x);
   }

   public double cdf (double x) {
      return cdf (m_lambda, x);
   }

   public double barF (double x) {
      return barF (m_lambda, x);
   }

   public double inverseF (double u) {
      return inverseF (m_lambda, u);
   }

   public double getMean() {
      return getMean (m_lambda);
   }

   public double getVariance() {
      return getVariance (m_lambda);
   }

   public double getStandardDeviation() {
      return getStandardDeviation (m_lambda);
   }

/**
 * Computes the density function @f$f(x)@f$, with @f$\lambda_i = @f$
 * <tt>lambda[</tt>@f$i-1@f$<tt>]</tt>, @f$i = 1,…,k@f$.
 *  @param lambda       rates of the hypoexponential distribution
 *  @param x            value at which the density is evaluated
 *  @return density at @f$x@f$
 */
public static double density (double[] lambda, double x) {
      if (x < 0)
         return 0;
      DoubleMatrix2D Ax = buildMatrix (lambda, x);
      DoubleMatrix2D T = DMatrix.expBidiagonal(Ax);
      int m = lambda.length;
      return lambda[m-1]*T.getQuick(0, m-1);
   }

   /**
    * Computes the distribution function @f$F(x)@f$, with @f$\lambda_i =
    * @f$ <tt>lambda[</tt>@f$i-1@f$<tt>]</tt>, @f$i = 1,…,k@f$.
    *  @param lambda       rates of the hypoexponential distribution
    *  @param x            value at which the distribution is evaluated
    *  @return distribution at @f$x@f$
    */
   public static double cdf (double[] lambda, double x) {
      if (x <= 0.0)
         return 0.0;
      if (x >= Double.MAX_VALUE)
         return 1.0;
      double mean = HypoExponentialDist.getMean (lambda);
      double std = HypoExponentialDist.getStandardDeviation(lambda);
      double LOW = mean - 1.5*std;
      if (x > LOW) {
         double p = 1.0 - HypoExponentialDist.barF(lambda, x);
         double LIMIT = 1.0e-3;
         if (p > LIMIT)
            return p;
      }

      DoubleMatrix2D T = buildMatrix (lambda, x);
      DoubleFactory1D fac1 = DoubleFactory1D.dense;
      int m = lambda.length;
      DoubleMatrix1D C = fac1.make(m, 1.0);
      DoubleMatrix1D B = DMatrix.expmiBidiagonal(T, C);
      return Math.abs(B.getQuick(0));
   }

   /**
    * Computes the distribution function @f$F(x)@f$, with @f$\lambda_i =
    * @f$ <tt>lambda[</tt>@f$i-1@f$<tt>]</tt>, @f$i = 1,…,k@f$. Returns
    * @f$1 - @f$<tt>barF(lambda, x)</tt>, which is much faster than `cdf`
    * but loses precision in the lower tail.
    *  @param lambda       rates of the hypoexponential distribution
    *  @param x            value at which the distribution is evaluated
    *  @return distribution at @f$x@f$
    */
   public static double cdf2 (double[] lambda, double x) {
      if (x <= 0.0)
         return 0.0;
      if (x >= Double.MAX_VALUE)
         return 1.0;
      return (1.0 - HypoExponentialDist.barF(lambda, x));
   }

   /**
    * Computes the complementary distribution @f$\bar{F}(x)@f$, with
    * @f$\lambda_i = @f$ <tt>lambda[</tt>@f$i-1@f$<tt>]</tt>, @f$i =
    * 1,…,k@f$.
    *  @param lambda       rates of the hypoexponential distribution
    *  @param x            value at which the complementary distribution
    *                      is evaluated
    *  @return complementary distribution at @f$x@f$
    */
   public static double barF (double[] lambda, double x) {
      if (x <= 0.0)
         return 1.0;
      if (x >= Double.MAX_VALUE)
         return 0.0;
      DoubleMatrix2D T = buildMatrix (lambda, x);
      DoubleFactory1D fac1 = DoubleFactory1D.dense;
      int m = lambda.length;
      DoubleMatrix1D C = fac1.make(m, 1.0);
      DoubleMatrix1D B = DMatrix.expBidiagonal(T, C);
      return B.getQuick(0);
   }

   /**
    * Computes the inverse distribution function @f$F^{-1}(u)@f$, with
    * @f$\lambda_i = @f$ <tt>lambda[</tt>@f$i-1@f$<tt>]</tt>, @f$i =
    * 1,…,k@f$.
    *  @param lambda       rates of the hypoexponential distribution
    *  @param u            value at which the inverse distribution is
    *                      evaluated
    *  @return inverse distribution at @f$u@f$
    */
   public static double inverseF (double[] lambda, double u) {
      if (u < 0.0 || u > 1.0)
          throw new IllegalArgumentException ("u not in [0,1]");
      if (u >= 1.0)
          return Double.POSITIVE_INFINITY;
      if (u <= 0.0)
          return 0.0;

      final double EPS = 1.0e-12;
      myFunc fonc = new myFunc (lambda, u);
      double x1 = getMean (lambda);
      double v = cdf (lambda, x1);
      if (u <= v)
         return RootFinder.brentDekker (0, x1, fonc, EPS);

      // u > v
      double x2 = 4.0*x1 + 1.0;
      v = cdf (lambda, x2);
      while (v < u) {
         x1 = x2;
         x2 = 4.0*x2;
         v = cdf (lambda, x2);
      }
      return RootFinder.brentDekker (x1, x2, fonc, EPS);
   }

   /**
    * Returns the mean, @f$E[X] = \sum_{i=1}^k 1/\lambda_i@f$, of the
    * hypoexponential distribution with rates @f$\lambda_i = @f$
    * <tt>lambda[</tt>@f$i-1@f$<tt>]</tt>, @f$i = 1,…,k@f$.
    *  @param lambda       rates of the hypoexponential distribution
    *  @return mean of the hypoexponential distribution
    */
   public static double getMean (double[] lambda) {
      testLambda (lambda);
      int k = lambda.length;
      double sum = 0;
      for (int j = 0; j < k; j++)
         sum += 1.0 / lambda[j];
      return sum;
   }

   /**
    * Returns the variance, @f$\mbox{Var}[X] = \sum_{i=1}^k
    * 1/\lambda_i^2@f$, of the hypoexponential distribution with rates
    * @f$\lambda_i = @f$ <tt>lambda[</tt>@f$i-1@f$<tt>]</tt>, @f$i =
    * 1,…,k@f$.
    *  @param lambda       rates of the hypoexponential distribution
    *  @return variance of the hypoexponential distribution
    */
   public static double getVariance (double[] lambda) {
      testLambda (lambda);
      int k = lambda.length;
      double sum = 0;
      for (int j = 0; j < k; j++)
         sum += 1.0 / (lambda[j]*lambda[j]);
      return sum;
   }

   /**
    * Returns the standard deviation of the hypoexponential distribution
    * with rates @f$\lambda_i = @f$ <tt>lambda[</tt>@f$i-1@f$<tt>]</tt>,
    * @f$i = 1,…,k@f$.
    *  @param lambda       rates of the hypoexponential distribution
    *  @return standard deviation of the hypoexponential distribution
    */
   public static double getStandardDeviation (double[] lambda) {
      double s = getVariance (lambda);
      return Math.sqrt(s);
   }

   /**
    * Returns the values @f$\lambda_i@f$ for this object.
    */
   public double[] getLambda() {
      return m_lambda;
   }

   /**
    * Sets the values @f$\lambda_i =
    * @f$<tt>lambda[</tt>@f$i-1@f$<tt>]</tt>, @f$i = 1,…,k@f$ for this
    * object.
    */
   public void setLambda (double[] lambda) {
      if (lambda == null)
         return;
      int k = lambda.length;
      m_lambda = new double[k];
      testLambda (lambda);
      System.arraycopy (lambda, 0, m_lambda, 0, k);
   }

   /**
    * Same as  #getLambda.
    */
   public double[] getParams() {
      return m_lambda;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      StringBuilder sb = new StringBuilder();
      Formatter formatter = new Formatter(sb, Locale.US);
      formatter.format(getClass().getSimpleName() + " : lambda = {" +
           PrintfFormat.NEWLINE);
      int k = m_lambda.length;
      for(int i = 0; i < k; i++) {
         formatter.format("   %g%n", m_lambda[i]);
      }
      formatter.format("}%n");
      return sb.toString();
   }

}