/*
 * Class:        DiscShiftBaker1
 * Description:  computes a discrepancy for a randomly shifted, then baker
                 folded point set
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
 * This class computes the **discrepancy** for *randomly shifted, then baker
 * folded* points of a set @f$\mathcal{P}@f$. It is given by
 * @cite vHIC02a&thinsp; (eq. 15)
 * @anchor REF_discrepancy_DiscShiftBaker1_baker1
 * @f{align}{
 *    [\mathcal{D}(\mathcal{P})]^2 
 *    & 
 *   =
 *    -1 + \frac{1}{n^2} \sum_{i=1}^n \sum_{j=1}^n \prod_{r=1}^s \left(1 - \frac{4\gamma_r^2}{3} \left[B_4(\{x_{ir} - x_{jr}\}) - B_4(\{\{x_{ir} - x_{jr}\}-1/2\})\right]\right. \nonumber
 *    \\  &   
 *    \qquad\qquad{} - \frac{\gamma_r^4}{9} \left[7B_4(\{x_{ir} - x_{jr}\}) - 2B_4(\{\{x_{ir} - x_{jr}\}-1/2\})\right] \tag{baker1} 
 *    \\  &   
 *    \qquad\qquad{} \left. {} - \frac{16\gamma_r^4}{45} \left[B_6(\{x_{ir} - x_{jr}\}) -B_6(\{\{x_{ir} - x_{jr}\}-1/2\})\right] \right)\nonumber,
 * @f}
 * where @f$n@f$ is the number of points of @f$\mathcal{P}@f$, @f$s@f$ is the
 * dimension of the points, @f$x_{ir}@f$ is the @f$r@f$-th coordinate of
 * point @f$i@f$, and the @f$\gamma_r@f$ are arbitrary positive weights. The
 * @f$B_{\alpha}(x)@f$ are the Bernoulli polynomials @cite mABR70a&thinsp;
 * (chap. 23) of degree @f$\alpha@f$ (see
 * umontreal.ssj.util.Num.bernoulliPoly in class <tt>util/Num</tt>), and the
 * notation @f$\{x\}@f$ means the fractional part of @f$x@f$, defined here as
 * @f$\{x\} = x \bmod1@f$. In one dimension, the formula simplifies to
 * @anchor REF_discrepancy_DiscShiftBaker1_baker1dim1
 * @f{align}{
 *    [\mathcal{D}(\mathcal{P})]^2 
 *    & 
 *   =
 *    \frac{1}{n^2} \sum_{i=1}^n \sum_{j=1}^n \left[ - \frac{4\gamma^2}{3} \left[B_4(\{x_i - x_j\}) - B_4(\{\{x_i - x_j\}-1/2\})\right]\right. \nonumber
 *    \\  &   
 *    \qquad\qquad{} - \frac{\gamma^4}{9} \left[7B_4(\{x_i - x_j\}) - 2B_4(\{\{x_i - x_j\}-1/2\})\right] \tag{baker1dim1} 
 *    \\  &   
 *    \qquad\qquad{} \left. {} - \frac{16\gamma^4}{45} \left[B_6(\{x_i - x_j\}) -B_6(\{\{x_i - x_j\}-1/2\})\right] \right]\nonumber.
 * @f}
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
public class DiscShiftBaker1 extends Discrepancy {

   static protected double []  setC (double gam) { // dimension = 1
      double[] C = new double[3];
      double v = gam*gam;
      C[0] = v * 4.0 / 3.0;
      v = v*v;
      C[1] = v / 9.0;
      C[2] = v * 16.0 / 45.0;
      return C;
   }


   static protected void setC (double[] C1, double[] C2, double[] C3,
                               double[] gam, int s) {
      for (int i = 0; i < s; i++) {
         double v = gam[i]*gam[i];
         C1[i] = v * 4.0 / 3.0;
         v = v*v;
         C2[i] = v / 9.0;
         C3[i] = v * 16.0 / 45.0;
      }
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in @f$s@f$
    * dimensions, with all the weights @f$\gamma_r = 1@f$. `points[i][r]`
    * is the <tt>r</tt>-th coordinate of point `i`. Indices `i` and `r`
    * start at 0.
    */
   public DiscShiftBaker1 (double[][] points, int n, int s) {
      super (points, n, s);
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in @f$s@f$
    * dimensions, with weights @f$\gamma_r = @f$ `gamma[r-1]`.
    * `points[i][r]` is the <tt>r</tt>-th coordinate of point `i`. Indices
    * `i` and `r` start at 0.
    */
   public DiscShiftBaker1 (double[][] points, int n, int s, double[] gamma) {
      super (points, n, s, gamma);
   }

   /**
    * The number of points is @f$n@f$, the dimension @f$s@f$, and the
    * @f$s@f$ weight factors are <tt>gamma[</tt>@f$r@f$<tt>]</tt>, @f$r =
    * 0, 1, …, (s-1)@f$. The @f$n@f$ points will be chosen later.
    */
   public DiscShiftBaker1 (int n, int s, double[] gamma) {
      super (n, s, gamma);
   }

   /**
    * Constructor with the point set `set`. All the points are copied in
    * an internal array.
    */
   public DiscShiftBaker1 (PointSet set) {
      super(set);
   }

   /**
    * Empty constructor. One *must* set the points, the dimension, and the
    * weight factors before calling any method.
    */
   public DiscShiftBaker1() {
   }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShiftBaker1_baker1 baker1}
    * ) for the @f$s@f$-dimensional points of set `points`, containing
    * @f$n@f$ points. All weights @f$\gamma_r = 1@f$.
    */
   public double compute (double[][] points, int n, int s) {
      setONES (s);
      return compute (points, n, s, ONES);
   }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShiftBaker1_baker1 baker1}
    * ) for the first @f$n@f$ points of `points` in dimension @f$s@f$ and
    * with weight @f$\gamma_r = @f$ `gamma[r-1]`.
    */
   public double compute (double[][] points, int n, int s, double[] gamma) {
      double[] C1 = new double[s];
      double[] C2 = new double[s];
      double[] C3 = new double[s];
      setC (C1, C2, C3, gamma, s);

      double pol1 = -1.0/30.0;        // BernoulliPoly(4, 0);
      double pol2 =  7.0/240.0;       // BernoulliPoly(4, 0.5);
      double pol3 =  1.0/42.0;        // BernoulliPoly(6, 0);
      double pol4 = -31.0/1344.0;     // BernoulliPoly(6, 0.5);

      double prod = 1.0;
      double temp;
      for (int r = 0; r < s; ++r) {
         temp = C1[r] * (pol1 - pol2) +
                C2[r] * (7.0*pol1 - 2.0*pol2) +
                C3[r] * (pol3 - pol4);
         prod *= 1.0 - temp;
      }
      double disc = prod / n;

      double sum = 0.0;
      for (int i = 0; i < n - 1; ++i) {
         for (int j = i + 1; j < n; ++j) {
            prod = 1.0;
            for (int r = 0; r < s; ++r) {
                double u = points[i][r] - points[j][r];
                if (u < 0.0)
                   u += 1.0;
                pol1 = ((u - 2.0) * u + 1.0)*u*u - UNTRENTE; // Bernoulli(4,u)
                double v = u - 0.5;
                if (v < 0.0)
                   v += 1.0;
                pol2 = ((v - 2.0) * v + 1.0)*v*v - UNTRENTE; // Bernoulli(4,v)

                // Bernoulli (6, u) and Bernoulli (6, v)
                pol3 = (((u - 3.0) * u + 2.5) * u*u - 0.5) * u*u + QUARAN;
                pol4 = (((v - 3.0) * v + 2.5) * v*v - 0.5) * v*v + QUARAN;

                temp = C1[r] * (pol1 - pol2) +
                       C2[r] * (7.0*pol1 - 2.0*pol2) +
                       C3[r] * (pol3 - pol4);
                prod *= 1.0 - temp;
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
    * {@link REF_discrepancy_DiscShiftBaker1_baker1dim1
    * baker1dim1} ) for the first @f$n@f$ points of @f$T@f$ in 1
    * dimension, with weight @f$\gamma= 1@f$.
    */
   public double compute (double[] T, int n) {
      return compute (T, n, 1.0);
    }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShiftBaker1_baker1dim1
    * baker1dim1} ) for the first @f$n@f$ points of @f$T@f$ in 1
    * dimension, with weight @f$\gamma=@f$ `gamma`.
    */
   public double compute (double[] T, int n, double gamma) {
      double[] C = setC(gamma);
      double C1 = C[0];
      double C2 = C[1];
      double C3 = C[2];

      double pol1 = -1.0/30.0;        // BernoulliPoly(4, 0);
      double pol2 =  7.0/240.0;       // BernoulliPoly(4, 0.5);
      double temp = C1 * (pol1 - pol2) +
                    C2 * (7.0*pol1 - 2.0*pol2) +
                    C3 * (1.0/42.0 + 31.0/1344.0);
      double disc = -temp / n;

      double sum = 0.0;
      for (int i = 0; i < n - 1; ++i) {
         for (int j = i + 1; j < n; ++j) {
            double h = T[i] - T[j];
            if (h < 0.0)
               h += 1.0;
            pol1 = ((h - 2.0) * h + 1.0)*h*h - UNTRENTE;    // Bernoulli(4,h)
            double v = h - 0.5;
            if (v < 0.0)
               v += 1.0;
            pol2 = ((v - 2.0) * v + 1.0)*v*v - UNTRENTE;    // Bernoulli(4,v)
            temp = C1 * (pol1 - pol2) +
                   C2 * (7.0*pol1 - 2.0*pol2);
            // Bernoulli (6, h) and Bernoulli (6, v)
            pol1 = (((h - 3.0) * h + 2.5) * h*h - 0.5) * h*h + QUARAN;
            pol2 = (((v - 3.0) * v + 2.5) * v*v - 0.5) * v*v + QUARAN;
            temp += C3 * (pol1 - pol2);
            sum -= temp;
          }
      }

      disc += 2.0*sum / ((long)n*n);
      if (disc < 0.0)
         return -1.0;
      return Math.sqrt(disc);
    }

}