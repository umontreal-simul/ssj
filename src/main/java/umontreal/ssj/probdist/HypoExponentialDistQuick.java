/*
 * Class:        HypoExponentialDistQuick
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
import umontreal.ssj.util.*;
import umontreal.ssj.functions.MathFunction;

/**
 * This class is a subclass of  @ref HypoExponentialDist and also implements
 * the *hypoexponential* distribution. It uses different algorithms to
 * compute the probabilities. The formula (
 * {@link REF_probdist_HypoExponentialDist_eq_tail_hypoexp
 * tail-hypoexp} ) for the complementary distribution is mathematically
 * equivalent to (see @cite pROS07b&thinsp; (page 299) and
 * @cite pGER10a&thinsp; (Appendix B))
 * @anchor REF_probdist_HypoExponentialDistQuick_eq_convolution_hypo
 * @f[
 *   \bar{F}(x) = \mathbb P\left[X_1 + \cdots+ X_k > x \right] = \sum_{i=1}^k e^{-\lambda_i x} \prod_{\substack {j=1\\j\not i}}^k \frac{\lambda_j}{\lambda_j - \lambda_i}. \tag{convolution-hypo}
 * @f]
 * The expression (
 * {@link REF_probdist_HypoExponentialDistQuick_eq_convolution_hypo
 * convolution-hypo} ) is much faster to compute than the matrix exponential
 * formula (
 * {@link REF_probdist_HypoExponentialDist_eq_tail_hypoexp
 * tail-hypoexp} ), but it becomes numerically unstable when @f$k@f$ gets
 * large and/or the differences between the @f$\lambda_i@f$ are too small,
 * because it is an alternating sum with relatively large terms of similar
 * size. When the @f$\lambda_i@f$ are close, many of the factors
 * @f$\lambda_j - \lambda_i@f$ in (
 * {@link REF_probdist_HypoExponentialDistQuick_eq_convolution_hypo
 * convolution-hypo} ) are small, and the effect of this is amplified when
 * @f$k@f$ is large. This gives rise to large terms of opposite sign in the
 * sum and the formula becomes unstable due to subtractive cancellation. For
 * example, with the computations done in standard 64-bit floating-point
 * arithmetic, if the @f$\lambda_i@f$ are regularly spaced with differences
 * of @f$\lambda_{i+1} - \lambda_i = 0.1@f$ for all @f$i@f$, the formula (
 * {@link REF_probdist_HypoExponentialDistQuick_eq_convolution_hypo
 * convolution-hypo} ) breaks down already for @f$k \approx15@f$, while if
 * the differences @f$\lambda_{i+1} - \lambda_i = 3@f$, it gives a few
 * decimal digits of precision for @f$k@f$ up to @f$\approx300@f$.
 *
 * The formula (
 * {@link REF_probdist_HypoExponentialDist_eq_fhypoexp
 * fhypoexp} ) for the density is mathematically equivalent to the much
 * faster formula
 * @anchor REF_probdist_HypoExponentialDistQuick_eq_fhypoexp2
 * @f[
 *   f(x) = \sum_{i=1}^k\lambda_i e^{-\lambda_i x} \prod_{\substack {j=1\\j\not i}}^k \frac{\lambda_j}{\lambda_j - \lambda_i}, \tag{fhypoexp2}
 * @f]
 * which is also numerically unstable when @f$k@f$ gets large and/or the
 * differences between the @f$\lambda_i@f$ are too small.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class HypoExponentialDistQuick extends HypoExponentialDist {
   private double[] m_H;

   private static double[] computeH (double[] lambda) {
      int k = lambda.length;
      double[] H = new double[k];
      double tem;
      int j;
      for (int i = 0; i < k; i++) {
         tem = 1.0;
         for (j = 0; j < i; j++)
            tem *= lambda[j] / (lambda[j] - lambda[i]);
         for (j = i + 1; j < k; j++)
            tem *= lambda[j] / (lambda[j] - lambda[i]);
         H[i] = tem;
      }
      return H;
   }


   private static double m_density(double[] lambda, double[] H, double x) {
      if (x < 0)
         return 0;

      int k = lambda.length;
      double tem;
      double prob = 0;
      for (int j = 0; j < k; j++) {
         tem = Math.exp (-lambda[j] * x);
         if (tem > 0)
            prob += lambda[j] * H[j] * tem;
      }

      return prob;
   }


   private static double m_barF(double[] lambda, double[] H, double x) {
      if (x <= 0.)
         return 1.;

      int k = lambda.length;
      double tem;
      double prob = 0;            // probability
      for (int j = 0; j < k; j++) {
         tem = Math.exp (-lambda[j] * x);
         if (tem > 0)
            prob += H[j] * tem;
      }
      return prob;
	}


   private static double m_cdf(double[] lambda, double[] H, double x) {
      if (x <= 0.)
         return 0.;

      int k = lambda.length;
      double tem = Math.exp(-lambda[0] * x);
      if (tem <= 0)
         return 1.0 - HypoExponentialDistQuick.m_barF(lambda, H, x);

      double prob = 0;            // cumulative probability
      for (int j = 0; j < k; j++) {
         tem = Math.expm1 (-lambda[j] * x);
         prob += H[j] * tem;
      }
      return -prob;
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
         return m_u - HypoExponentialDistQuick.cdf(m_lam, x);
      }
   }

   /**
    * Constructs a `HypoExponentialDistQuick` object, with rates
    * @f$\lambda_i = @f$ <tt>lambda[</tt>@f$i-1@f$<tt>]</tt>, @f$i =
    * 1,…,k@f$.
    *  @param lambda       rates of the hypoexponential distribution
    */
   public HypoExponentialDistQuick (double[] lambda) {
      super(lambda);
   }

  /*  These public methods are necessary so that methods cdf,
   *  barF and inverseF used are those of the present
   *  class and not those of the mother class.
   */
   public double density (double x) {
     return m_density (m_lambda, m_H, x);
   }

   public double cdf (double x) {
      return m_cdf (m_lambda, m_H, x);
   }

   public double barF (double x) {
      return m_barF (m_lambda, m_H, x);
   }

   public double inverseF (double u) {
      return m_inverseF (m_lambda, m_H, u);
   }

/**
 * Computes the density function @f$f(x)@f$, with @f$\lambda_i = @f$
 * <tt>lambda[</tt>@f$i-1@f$<tt>]</tt>, @f$i = 1,…,k@f$.
 *  @param lambda       rates of the hypoexponential distribution
 *  @param x            value at which the density is evaluated
 *  @return density at @f$x@f$
 */
public static double density (double[] lambda, double x) {
      testLambda (lambda);
      double[] H = computeH (lambda);
      return m_density(lambda, H, x);
   }

   /**
    * Computes the distribution function @f$F(x)@f$, with @f$\lambda_i =
    * @f$ <tt>lambda[</tt>@f$i-1@f$<tt>]</tt>, @f$i = 1,…,k@f$.
    *  @param lambda       rates of the hypoexponential distribution
    *  @param x            value at which the distribution is evaluated
    *  @return value of distribution at @f$x@f$
    */
   public static double cdf (double[] lambda, double x) {
      testLambda (lambda);
      double[] H = computeH (lambda);
      return m_cdf(lambda, H, x);
   }

   /**
    * Computes the complementary distribution @f$\bar{F}(x)@f$, with
    * @f$\lambda_i = @f$ <tt>lambda[</tt>@f$i-1@f$<tt>]</tt>, @f$i =
    * 1,…,k@f$.
    *  @param lambda       rates of the hypoexponential distribution
    *  @param x            value at which the complementary distribution
    *                      is evaluated
    *  @return value of complementary distribution at @f$x@f$
    */
   public static double barF (double[] lambda, double x) {
      testLambda (lambda);
      double[] H = computeH (lambda);
      return m_barF(lambda, H, x);
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
      testLambda (lambda);
      double[] H = computeH (lambda);
      return m_inverseF(lambda, H, u);
   }

   private static double m_inverseF(double[] lambda, double[] H, double u) {
      if (u < 0.0 || u > 1.0)
          throw new IllegalArgumentException ("u not in [0,1]");
      if (u >= 1.0)
          return Double.POSITIVE_INFINITY;
      if (u <= 0.0)
          return 0.0;

      final double EPS = 1.0e-12;
      myFunc fonc = new myFunc (lambda, u);
      double x1 = getMean (lambda);
      double v = m_cdf(lambda, H, x1);
      if (u <= v)
         return RootFinder.brentDekker (0, x1, fonc, EPS);

      // u > v
      double x2 = 4.0*x1 + 1.0;
      v = m_cdf(lambda, H, x2);
      while (v < u) {
         x1 = x2;
         x2 = 4.0*x2;
         v = m_cdf(lambda, H, x2);
      }
      return RootFinder.brentDekker (x1, x2, fonc, EPS);
   }


   public void setLambda (double[] lambda) {
      if (lambda == null)
         return;
      super.setLambda (lambda);
      m_H = computeH (lambda);
   }

   public String toString() {
      StringBuilder sb = new StringBuilder();
      Formatter formatter = new Formatter(sb, Locale.US);
      formatter.format(getClass().getSimpleName() + " : lambda = {" +
           PrintfFormat.NEWLINE);
      int k = m_lambda.length;
      for(int i = 0; i < k; i++) {
         formatter.format("   %f%n", m_lambda[i]);
      }
      formatter.format("}%n");
      return sb.toString();
   }
}