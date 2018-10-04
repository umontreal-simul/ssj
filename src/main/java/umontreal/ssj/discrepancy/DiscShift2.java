/*
 * Class:        DiscShift2
 * Description:  computes a discrepancy for the randomly shifted points of a set
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        January 2009

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
package umontreal.ssj.discrepancy;
   import umontreal.ssj.hups.PointSet;
import umontreal.ssj.util.Num;

/**
 * This class computes the **discrepancy** in @cite vHIC02a&thinsp; (eq. 15)
 * for the *randomly shifted* points of a set @f$\mathcal{P}@f$. It is given
 * by
 * @anchor REF_discrepancy_DiscShift2_shift2
 * @f[
 *   [\mathcal{D}(\mathcal{P})]^2 = -1 + \frac{1}{n^2} \sum_{i=1}^n \sum_{j=1}^n \prod_{r=1}^s \left[1 + \frac{\gamma_r^2}{2} B_2(\{x_{ir} - x_{jr}\}) - \frac{\gamma_r^4}{12}B_4(\{x_{ir} - x_{jr}\}) \right], \tag{shift2}
 * @f]
 * where @f$n@f$ is the number of points of @f$\mathcal{P}@f$, @f$s@f$ is the
 * dimension of the points, @f$x_{ir}@f$ is the @f$r@f$-th coordinate of
 * point @f$i@f$, and the @f$\gamma_r@f$ are arbitrary positive weights. The
 * @f$B_{\alpha}(x)@f$ are the Bernoulli polynomials @cite mABR70a&thinsp;
 * (chap. 23) of degree @f$\alpha@f$ (see
 * umontreal.ssj.util.Num.bernoulliPoly in class <tt>util/Num</tt>). The
 * notation @f$\{x\}@f$ means the fractional part of @f$x@f$, defined as
 * @f$\{x\} = x \bmod1@f$. In one dimension, the formula simplifies to
 * @anchor REF_discrepancy_DiscShift2_shift2dim1
 * @f[
 *   [\mathcal{D}(\mathcal{P})]^2 = \frac{1}{n^2} \sum_{i=1}^n \sum_{j=1}^n \left[\frac{\gamma^2}{2} B_2(\{x_i - x_j\}) - \frac{\gamma^4}{12}B_4(\{x_i - x_j\}) \right]. \tag{shift2dim1}
 * @f]
 * The discrepancy represents a worst-case error criterion for the
 * approximation of integrals, when the integrands have a certain degree of
 * smoothness and lie in a Hilbert space @f$\mathcal{H}@f$ with a reproducing
 * kernel @f$K@f$ given by
 * @f[
 *   K(\mathbf{x},\mathbf{y}) = \prod_{r=1}^s \left[ - \frac{\gamma_r^4}{4!} B_4(\{x_r-y_r\}) + \sum_{\alpha=0}^2 \frac{\gamma_r^{2\alpha}}{(\alpha!)^2} B_{\alpha}(x_r)B_{\alpha}(y_r) \right],
 * @f]
 * The norm of the vectors in @f$\mathcal{H}@f$ is defined by
 * @f[
 *   \| f\|^2 = \sum_{u \subseteq S}\sum_{v \subseteq u}\gamma_u^{-2}\gamma_v^{-2} \int_{[0,1]^v}d\mathbf{x}_v\left[\int_{[0,1]^{S-v}} \frac{\partial^{|u| + |v|}f}{\partial\mathbf{x}_u\partial\mathbf{x}_v} d\mathbf{x}_{S-v} \right]^2,
 * @f]
 * where @f$S= \{1, …, s\}@f$ is a set of coordinate indices, @f$u
 * \subseteq S@f$, and @f$\gamma_u = \prod_{r\in u} \gamma_r@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class DiscShift2 extends Discrepancy {

   protected static double[] setC (double gam) {
      // 1 dimension
      double[] C = new double[2];
      double v = gam*gam;
      C[0] = 0.5 * v;
      v = v*v;
      C[1] = v / 12.0;
      return C;
   }

   protected static void setC (double[] C1, double[] C2, double[] gam, int s) {
      for (int i = 0; i < s; i++) {
         double v = gam[i]*gam[i];
         C1[i] = 0.5 * v;
         v = v*v;
         C2[i] = v / 12.0;
      }
   }

   /**
    * Constructor with the @f$n@f$ points @f$P_i = @f$ `points[i]` in
    * dimension @f$s@f$, with all weights @f$\gamma_j =1@f$.
    * `points[i][j]` is the <tt>j</tt>-th coordinate of point `i`. Both
    * indices `i` and `j` starts at 0.
    */
   public DiscShift2 (double[][] points, int n, int s) {
      super (points, n, s);
   }

   /**
    * Constructor with the @f$n@f$ points @f$P_i = @f$ `points[i]` in
    * dimension @f$s@f$, with the weights @f$\gamma_j = @f$ `gamma[j-1]`,
    * @f$j = 1, …, s@f$. `points[i][j]` is the <tt>j</tt>-th coordinate of
    * point `i`. Both indices `i` and `j` starts at 0.
    */
   public DiscShift2 (double[][] points, int n, int s, double[] gamma) {
      super (points, n, s, gamma);
   }

   /**
    * The number of points is @f$n@f$, the dimension @f$s@f$, and the
    * @f$s@f$ weight factors are <tt>gamma[</tt>@f$j@f$<tt>]</tt>, @f$j =
    * 0, 1, …, (s-1)@f$. The @f$n@f$ points will be chosen later.
    */
   public DiscShift2 (int n, int s, double[] gamma) {
      super (n, s, gamma);
   }

   /**
    * Constructor with the point set `set`. All the points are copied in
    * an internal array.
    */
   public DiscShift2 (PointSet set) {
      super(set);
   }

   /**
    * Empty constructor. One *must* set the points, the dimension, and the
    * weight factors before calling any method.
    */
   public DiscShift2() {
   }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShift2_shift2 shift2} ) for
    * the first @f$n@f$ points of set `points` in dimension @f$s@f$. All
    * weights @f$\gamma_r = 1@f$.
    */
   public double compute (double[][] points, int n, int s) {
      setONES (s);
      return compute (points, n, s, ONES);
   }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShift2_shift2 shift2} ) for
    * the first @f$n@f$ points of set `points` in dimension @f$s@f$ and
    * with weight @f$\gamma_r = @f$ `gamma[r-1]`.
    */
   public double compute (double[][] points, int n, int s, double[] gamma) {
      double[] C1 = new double[s];
      double[] C2 = new double[s];
      setC (C1, C2, gamma, s);

      double pol1 = UNSIX;            // BernoulliPoly(2, 0);
      double pol2 = -UNTRENTE;        // BernoulliPoly(4, 0);
      double prod = 1.0;
      for (int r = 0; r < s; ++r)
         prod *= (1.0 + C1[r] * pol1 - C2[r] * pol2);
      double disc = prod / n;

      double sum = 0.0;
      for (int i = 0; i < n - 1; ++i) {
         for (int j = i + 1; j < n; ++j) {
            prod = 1.0;
            for (int r = 0; r < s; ++r) {
                double u = points[i][r] - points[j][r];
                if (u < 0.0)
                   u += 1.0;
                pol1 = u*(u - 1.0) + UNSIX;   // Bernoulli(2, u)
                pol2 = ((u - 2.0) * u + 1.0)*u*u - UNTRENTE;    // Bernoulli(4,h)
                prod *= 1.0 + C1[r] * pol1 - C2[r] * pol2;
            }
            sum += prod;
         }
      }

      disc += 2.0 * sum / ((long)n*n) - 1.0;
      if (disc < 0.0)
         return -1.0;
      return Math.sqrt(disc);
    }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShift2_shift2dim1
    * shift2dim1} ) for the first @f$n@f$ points of @f$T@f$ in 1
    * dimension, with weight @f$\gamma= 1@f$.
    */
   public double compute (double[] T, int n) {
      return compute (T, n, 1.0);
    }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShift2_shift2dim1
    * shift2dim1} ) for the first @f$n@f$ points of @f$T@f$ in 1
    * dimension, with weight @f$\gamma=@f$ `gamma`.
    */
   public double compute (double[] T, int n, double gamma) {
      double pol1 = UNSIX;            // BernoulliPoly(2, 0);
      double pol2 = -UNTRENTE;        // BernoulliPoly(4, 0);
      double[] C = setC(gamma);
      double C1 = C[0];
      double C2 = C[1];
      double disc = (C1 * pol1 - C2 * pol2) / n;

      double sum = 0.0;
      for (int i = 0; i < n - 1; ++i) {
         for (int j = i + 1; j < n; ++j) {
            double u = T[i] - T[j];
            if (u < 0.0)
                u += 1.0;
            pol1 = u*(u - 1.0) + UNSIX;   // Bernoulli(2, u)
            pol2 = ((u - 2.0) * u + 1.0)*u*u - UNTRENTE;    // Bernoulli(4,h)
            sum += (C1 * pol1 - C2 * pol2);
         }
      }

      disc += 2.0 * sum / ((long)n*n);
      if (disc < 0.0)
         return -1.0;
      return Math.sqrt(disc);
    }

}