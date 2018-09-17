/*
 * Class:        DiscL2Hickernell
 * Description:  computes the Hickernell L_2 star discrepancy
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

/**
 * This class computes the Hickernell @f$\mathcal{L}_2@f$-<strong>star
 * discrepancy</strong> in @cite vHIC98a&thinsp; (eq. 5.1c) for a point set.
 * It is based on the reproducing kernel Hilbert space which consists of
 * functions with square-integrable first-order derivatives, and whose
 * reproducing kernel is given by
 * @f[
 *   K(\mathbf{x},\mathbf{y}) = \prod_{j=1}^s \eta_j(x_j,y_j),
 * @f]
 * where
 * @f[
 *   \eta_j(x,y) = 1 + \frac{|x-a_j| + |y - a_j| - |x-y|}{2}.
 * @f]
 * The inner product is given by
 * @f[
 *   \langle f, g \rangle= \sum_{\mathfrak u\subseteq S} \int_{[0, 1)^{|\mathfrak u|}} \frac{\partial^{|\mathfrak u|}f(\mathbf{x}_{\mathfrak u}, \boldsymbol a)}{\partial\mathbf{x_{\mathfrak u}}} \frac{\partial^{|\mathfrak u|}g(\mathbf{x}_{\mathfrak u}, \boldsymbol a)}{\partial\mathbf{x_{\mathfrak u}}} d\mathbf{x_{\mathfrak u}},
 * @f]
 * where @f$S@f$ is the set of coordinates indices @f$\{1, 2, …, s\}@f$,
 * @f$|\mathfrak u|@f$ denotes the cardinality of @f$\mathfrak u@f$,
 * @f$\mathbf{x_{\mathfrak u}}@f$ denotes the vector made of the components
 * of @f$ \mathbf{x} @f$ whose indices are in @f$\mathfrak u@f$, and
 * @f$\boldsymbol a \in[0, 1]^s@f$ is called the anchor.
 *
 * The worst-case error for this function space is the @f$\mathcal{L}_2@f$
 * version of the star discrepancy @f$\mathcal{D}_2^*@f$. Choosing
 * @f$\boldsymbol a =\mathbf{1}@f$, Hickernell obtained the formula
 * @anchor REF_discrepancy_DiscL2Hickernell_disc_hicks
 * @f[
 *   [\mathcal{D}_2^*(\mathcal{P})]^2 = \left(\frac{4}{3}\right)^s - \frac{2}{n} \sum_{i=1}^n \prod_{k=1}^s \left(\frac{3 - z_{ik}^2}{2}\right) + \frac{1}{n^2} \sum_{i=1}^n\sum_{j=1}^n \prod_{k=1}^s \Bigl(2 - \max(z_{ik}, z_{jk})\Bigr), \tag{disc.hicks}
 * @f]
 * where @f$n@f$ is the number of points of set @f$\mathcal{P}@f$, @f$s@f$ is
 * the dimension of the points, and @f$z_{ik}@f$ is the @f$k@f$-th coordinate
 * of point @f$i@f$.
 *
 * In 1 dimension, the formula is equivalent to
 * @anchor REF_discrepancy_DiscL2Hickernell_disc_hickD1
 * @f[
 *   [\mathcal{D}_2^*(\mathcal{P})]^2 = \frac{1}{3} + \frac{1}{n} \sum_{i=1}^n {z_i^2} - \frac{1}{n^2} \sum_{i=1}^n\sum_{j=1}^n \max(z_i, z_j), \tag{hickD1}
 * @f]
 * where @f$z_i@f$ is the coordinate of point @f$i@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class DiscL2Hickernell extends Discrepancy {

   public double compute (double[][] points, int n, int s, double[] gamma) {
      return compute (points, n, s);
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in dimension
    * @f$s@f$. `points[i][j]` is the @f$j@f$-th coordinate of point
    * @f$i@f$. Both @f$i@f$ and @f$j@f$ start at 0.
    */
   public DiscL2Hickernell (double[][] points, int n, int s) {
      super (points, n, s);
   }

   /**
    * Constructor with @f$n@f$ points in dimension @f$s@f$. The @f$n@f$
    * points will be chosen later.
    */
   public DiscL2Hickernell (int n, int s) {
      super (null, n, s);
   }

   /**
    * Constructor with the point set `set`. All the points are copied in
    * an internal array.
    */
   public DiscL2Hickernell (PointSet set) {
      super(set);
   }

   /**
    * Empty constructor. One *must* set the points and the dimension
    * before calling any method.
    */
   public DiscL2Hickernell() {
      super();
   }

   /**
    * Computes the Hickernell
    * @f$\mathcal{L}_2@f$-<strong>discrepancy</strong> (
    * {@link REF_discrepancy_DiscL2Hickernell_disc_hicks
    * disc.hicks} ) for the set of @f$n@f$ @f$s@f$-dimensional points
    * `points`.
    */
   public double compute (double[][] points, int n, int s) {
      final int N = n;
      final double RAC3 = 1.73205080756887729352;  // sqrt(3)
      double sum = 0.0;
      for (int i = 0; i < N; ++i) {
         double prod = 1.0;
         for (int j = 0; j < s; ++j) {
             double u = points[i][j];
             prod *= (RAC3 - u) * (RAC3 + u);
         }
         sum += prod;
      }
      double disc = -Math.pow(0.5, (double) (s-1))*sum/N;

      sum = 0.0;
      for (int i = 0; i < N; ++i) {
         double prod = 1.0;
         for (int j = 0; j < s; ++j)
             prod *= 2.0 - points[i][j];
         sum += prod;
      }

      double sum2 = 0.0;
      for (int i = 0; i < N - 1; ++i) {
         for (int j = i + 1; j < N; ++j) {
            double prod = 1.0;
            for (int k = 0; k < s; ++k)
                prod *= 2.0 - Math.max(points[i][k], points[j][k]);
            sum2 += prod;
         }
      }

      disc += (sum + 2.0*sum2)/((long)N*N);
      disc += Math.pow(4.0/3.0, s);
      if (disc < 0.0)
         return -1.0;
      return Math.sqrt(disc);
    }

   /**
    * Computes the Hickernell
    * @f$\mathcal{L}_2@f$-<strong>discrepancy</strong> (
    * {@link REF_discrepancy_DiscL2Hickernell_disc_hickD1
    * hickD1} ) for the set of @f$n@f$ 1-dimensional points `T`.
    */
   public double compute (double[] T, int n) {
      // In 1 dimension, L2*Hickernell = L2*
      double sum = 0.0;
      for (int i = 0; i < n; ++i)
         sum += T[i]*T[i];
      double disc = sum/n;

      sum = 0.0;
      for (int i = 0; i < n; ++i)
         sum += T[i];

      double sum2 = 0.0;
      for (int i = 0; i < n - 1; ++i) {
         for (int j = i + 1; j < n; ++j)
            sum2 += Math.max(T[i], T[j]);
      }

      disc -= (sum + 2.0*sum2)/((long)n*n);
      disc += 1.0/3.0;
      if (disc < 0.0)
         return -1.0;
      return Math.sqrt(disc);
    }

}