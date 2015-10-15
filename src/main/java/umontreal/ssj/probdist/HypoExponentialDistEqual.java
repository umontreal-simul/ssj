/*
 * Class:        HypoExponentialDistEqual
 * Description:  Hypo-exponential distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        February 2014

 * SSJ is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, either version 3 of the License, or
 * any later version.

 * SSJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * A copy of the GNU General Public License is available at
   <a href="http://www.gnu.org/licenses">GPL licence site</a>.
 */
package umontreal.ssj.probdist;
import java.util.Formatter;
import java.util.Locale;
import umontreal.ssj.util.*;
import umontreal.ssj.functions.MathFunction;

/**
 * This class implements the *hypoexponential* distribution for the case of
 * equidistant @f$\lambda_i = (n+1-i)h@f$. We have @f$\lambda_{i+1} -
 * \lambda_i = h@f$, with @f$h@f$ a constant, and @f$n \ge k@f$ are
 * integers.
 *
 * The formula (
 * {@link REF_probdist_HypoExponentialDistQuick_eq_convolution_hypo
 * convolution-hypo} ) becomes
 * @anchor REF_probdist_HypoExponentialDistEqual_eq_conv_hypo_equal
 * @f[
 *   \bar{F}(x) = \mathbb P\left[X_1 + \cdots+ X_k > x \right] = \sum_{i=1}^k e^{-(n+1-i)h x} \prod_{\substack {j=1\\j\not i}}^k \frac{n+1-j}{i - j}. \tag{conv-hypo-equal}
 * @f]
 * The formula (
 * {@link REF_probdist_HypoExponentialDistQuick_eq_fhypoexp2
 * fhypoexp2} ) for the density becomes
 * @anchor REF_probdist_HypoExponentialDistEqual_eq_fhypoexp3
 * @f[
 *   f(x) = \sum_{i=1}^k (n+1-i)h e^{-(n+1-i)h x} \prod_{\substack {j=1\\j\not i}}^k \frac{n+1-j}{i - j}. \tag{fhypoexp3}
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class HypoExponentialDistEqual extends HypoExponentialDist {
   private double m_h;
   private int m_n;
   private int m_k;

   /**
    * Constructor for equidistant rates. The rates are @f$\lambda_i =
    * (n+1-i)h@f$, for @f$i = 1,…,k@f$.
    *  @param n            largest rate is @f$nh@f$
    *  @param k            number of rates
    *  @param h            difference between adjacent rates
    */
   public HypoExponentialDistEqual (int n, int k, double h) {
      super (null);
      setParams (n, k, h);
   }

   public double density (double x) {
      return density (m_n, m_k, m_h, x);
   }

   public double cdf (double x) {
      return cdf (m_n, m_k, m_h, x);
   }

   public double barF (double x) {
      return barF (m_n, m_k, m_h, x);
   }

   public double inverseF (double u) {
      return inverseF (m_n, m_k, m_h, u);
   }

/**
 * Computes the density function @f$f(x)@f$, with the same arguments as in
 * the constructor.
 *  @param n            max possible number of @f$\lambda_i@f$
 *  @param k            effective number of @f$\lambda_i@f$
 *  @param h            step between two successive @f$\lambda_i@f$
 *  @param x            value at which the distribution is evaluated
 *  @return density at @f$x@f$
 */
public static double density (int n, int k, double h, double x) {
      if (x < 0)
         return 0;
      double r = -Math.expm1(-h*x);
      double v = BetaDist.density(k, n - k + 1, r);
      return h*v*Math.exp(-h*x);
   }

   /**
    * Computes the distribution function @f$F(x)@f$, with arguments as in
    * the constructor.
    *  @param n            max possible number of @f$\lambda_i@f$
    *  @param k            effective number of @f$\lambda_i@f$
    *  @param h            step between two successive @f$\lambda_i@f$
    *  @param x            value at which the distribution is evaluated
    *  @return value of distribution at @f$x@f$
    */
   public static double cdf (int n, int k, double h, double x) {
      if (x <= 0)
         return 0;
      double r = -Math.expm1(-h*x);
      double u = BetaDist.cdf(k, n - k + 1, r);
      return u;
   }

   /**
    * Computes the complementary distribution @f$\bar{F}(x)@f$, as in
    * formula (
    * {@link REF_probdist_HypoExponentialDistEqual_eq_conv_hypo_equal
    * conv-hypo-equal} ).
    *  @param n            max possible number of @f$\lambda_i@f$
    *  @param k            effective number of @f$\lambda_i@f$
    *  @param h            step between two successive @f$\lambda_i@f$
    *  @param x            value at which the complementary distribution
    *                      is evaluated
    *  @return value of complementary distribution at @f$x@f$
    */
   public static double barF (int n, int k, double h, double x) {
      if (x <= 0)
         return 1.0;
      double r = Math.exp(-h*x);
      double v = BetaDist.cdf(n - k + 1, k, r);
      return v;
   }

   /**
    * Computes the inverse distribution @f$x=F^{-1}(u)@f$, with arguments
    * as in the constructor.
    *  @param n            max possible number of @f$\lambda_i@f$
    *  @param k            effective number of @f$\lambda_i@f$
    *  @param h            step between two successive @f$\lambda_i@f$
    *  @param u            value at which the inverse distribution is
    *                      evaluated
    *  @return inverse distribution at @f$u@f$
    */
   public static double inverseF (int n, int k, double h, double u) {
      if (u < 0.0 || u > 1.0)
          throw new IllegalArgumentException ("u not in [0,1]");
      if (u >= 1.0)
          return Double.POSITIVE_INFINITY;
      if (u <= 0.0)
          return 0.0;

      double z = BetaDist.inverseF(k, n - k + 1, u);
      return -Math.log1p(-z) / h;
   }

   /**
    * Returns the three parameters of this hypoexponential distribution as
    * array @f$(n, k, h)@f$.
    *  @return parameters of the hypoexponential distribution
    */
   public double[] getParams() {
      double[] par = new double[]{m_n, m_k, m_h};
      return par;
   }


   public void setParams (int n, int k, double h) {
      m_n = n;
      m_k = k;
      m_h = h;
      m_lambda = new double[k];
      for(int i = 0; i < k; i++) {
         m_lambda[i] = (n - i)*h;
      }
   }


   public String toString() {
      StringBuilder sb = new StringBuilder();
      Formatter formatter = new Formatter(sb, Locale.US);
      formatter.format(getClass().getSimpleName() + " : params = {" +
           PrintfFormat.NEWLINE);
      formatter.format("   %d, %d, %f", m_n, m_k, m_h);
      formatter.format("}%n");
      return sb.toString();
   }
}