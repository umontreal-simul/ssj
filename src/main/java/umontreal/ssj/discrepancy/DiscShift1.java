/*
 * Class:        DiscShift1
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
 * This class computes the **discrepancy** for randomly shifted points of a
 * set @f$\mathcal{P}@f$ @cite vHIC02a&thinsp; (eq. 15). It is given by
 * @anchor REF_discrepancy_DiscShift1_shift1
 * @f[
 *   [\mathcal{D}(\mathcal{P})]^2 = -1 + \frac{1}{n^2} \sum_{i=1}^n\sum_{j=1}^n \prod_{r=1}^s \left[1 + \gamma_r^2 B_2(\{x_{ir} - x_{jr}\})\right], \tag{shift1}
 * @f]
 * where @f$n@f$ is the number of points of @f$\mathcal{P}@f$, @f$s@f$ is the
 * dimension of the points, @f$x_{ir}@f$ is the @f$r@f$-th coordinate of
 * point @f$i@f$, and the @f$\gamma_r@f$ are arbitrary positive weights. The
 * @f$B_{\alpha}(x)@f$ are the Bernoulli polynomials @cite mABR70a&thinsp;
 * (chap. 23) of degree @f$\alpha@f$ (see
 * umontreal.ssj.util.Num.bernoulliPoly in class <tt>util/Num</tt>), and the
 * notation @f$\{x\}@f$ means the fractional part of @f$x@f$, defined here as
 * @f$\{x\} = x \bmod1@f$. In one dimension, the formula simplifies to
 * @anchor REF_discrepancy_DiscShift1_shift1dim1
 * @f[
 *   [\mathcal{D}(\mathcal{P})]^2 = \frac{1}{n^2} \sum_{i=1}^n\sum_{j=1}^n B_2(\{x_i - x_j\}), \tag{shift1dim1}
 * @f]
 * where @f$z_i@f$ is the point @f$i@f$.
 *
 * The discrepancy represents a worst-case error criterion for the
 * approximation of integrals, when the integrands have a certain degree of
 * smoothness and lie in a Hilbert space @f$\mathcal{H}@f$ with a reproducing
 * kernel @f$K@f$ given by
 * @f[
 *   K(\mathbf{x},\mathbf{y}) = \prod_{r=1}^s \left[ \frac{\gamma_r^2}{2} B_2(\{x_r-y_r\}) + \sum_{\alpha=0}^1 \gamma_r^{2\alpha} B_{\alpha}(x_r)B_{\alpha}(y_r) \right],
 * @f]
 * The norm of the vectors in @f$\mathcal{H}@f$ is defined by
 * @f[
 *   \| f\|^2 = \sum_{u \subseteq S} \gamma_u^{-2} \int_{[0,1]^u}d\mathbf{x}_u\left[\int_{[0,1]^{S-u}} \frac{\partial^{|u|}f}{\partial\mathbf{x}_u}d\mathbf{x}_{S-u} \right]^2,
 * @f]
 * where @f$S= \{1, …, s\}@f$ is a set of coordinate indices, @f$u
 * \subseteq S@f$, and @f$\gamma_u = \prod_{r\in u} \gamma_r@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class DiscShift1 extends Discrepancy {

   static protected void setC (double[] C1, double[] gam, int s) {
      for (int i = 0; i < s; i++)
         C1[i] = gam[i]*gam[i];
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in @f$s@f$
    * dimensions and with all weights @f$\gamma_r =1@f$. `points[i][j]`
    * is the <tt>j</tt>-th coordinate of point `i`. Indices `i` and `j`
    * start at 0.
    */
   public DiscShift1 (double[][] points, int n, int s) {
      super (points, n, s);
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in @f$s@f$
    * dimensions, and with the weights @f$\gamma_r = @f$ `gamma[r-1]`,
    * @f$r = 1, …, s@f$. `points[i][j]` is the <tt>j</tt>-th coordinate of
    * point `i`. Indices `i` and `j` start at 0.
    */
   public DiscShift1 (double[][] points, int n, int s, double[] gamma) {
      super (points, n, s, gamma);
   }

   /**
    * The number of points is @f$n@f$, the dimension @f$s@f$, and the
    * @f$s@f$ weight factors are <tt>gamma[</tt>@f$j@f$<tt>]</tt>, @f$j =
    * 0, 1, …, (s-1)@f$. The @f$n@f$ points will be chosen later.
    */
   public DiscShift1 (int n, int s, double[] gamma) {
      super (n, s, gamma);
   }

   /**
    * Constructor with the point set `set`. All the points are copied in
    * an internal array.
    */
   public DiscShift1 (PointSet set) {
      super(set);
   }

   /**
    * Empty constructor. One *must* set the points, the dimension, and the
    * weight factors before calling any method.
    */
   public DiscShift1() {
   }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShift1_shift1 shift1} ) for
    * the first @f$n@f$ points of set `points` in dimension @f$s@f$. All
    * weights @f$\gamma_r = 1@f$.
    */
   public double compute (double[][] points, int n, int s) {
      setONES (s);
      return compute (points, n, s, ONES);
   }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShift1_shift1 shift1} ) in
    * dimension @f$s@f$ with @f$\gamma_r = @f$ `gamma[r-1]`.
    */
   public double compute (double[][] points, int n, int s, double[] gamma) {
      double[] C1 = new double[s];      // gamma_r^2
      setC (C1, gamma, s);

      double pol1 = UNSIX;        // BernoulliPoly(2, 0);
      double prod = 1.0;
      for (int r = 0; r < s; ++r)
         prod *= (1.0 + C1[r] * pol1);
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
                prod *= 1.0 + C1[r] * pol1;
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
    * {@link REF_discrepancy_DiscShift1_shift1dim1
    * shift1dim1} ) for the 1-dimensional set of @f$n@f$ points @f$T@f$.
    */
   public double compute (double[] T, int n) {
      double pol1 = UNSIX;            // BernoulliPoly(2, 0);
      double disc = pol1 / n;

      double sum = 0.0;
      for (int i = 0; i < n - 1; ++i) {
         for (int j = i + 1; j < n; ++j) {
            double h = T[i] - T[j];
            if (h < 0.0)
               h += 1.0;
            pol1 = h*(h - 1.0) + UNSIX;   // Bernoulli(2, h)
            sum += pol1;
          }
      }

      disc += 2.0 * sum / ((long)n*n);
      if (disc < 0.0)
         return -1.0;
      return Math.sqrt(disc);
    }

}