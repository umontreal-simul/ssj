/*
 * Class:        DiscL2Star
 * Description:  computes the traditional L_2 star discrepancy
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
 * This class computes the traditional @f$\mathcal{L}_2@f$-<strong>star
 * discrepancy</strong> @f$\mathcal{D}_2^*(\mathcal{P})@f$ for a set of
 * @f$n@f$ points @f$\mathcal{P}@f$ @cite rWAR72a, @cite rWAR95a,
 * @cite tHIC99a&thinsp;. This discrepancy is also known as the Cramér-von
 * Mises goodness-of-fit statistic. Let @f$t = [t_1, t_2, …, t_s]@f$ be the
 * coordinates of a point in @f$[0, 1]^s@f$ and define the volume of the
 * parallelepiped anchored at the origin as Vol@f$([0,t)) = t_1t_2 …t_s@f$.
 * Let the function @f$R_n(t, P)@f$ be defined as
 * @f[
 *   R_n(t, P) = \frac{1}{n} \sum_{j=1}^n I_{[0, t)^s}(z_j) - \mbox{Vol}([0,t)),
 * @f]
 * where @f$I_{[0, t)^s}@f$ is the indicator function of @f$[0, t)^s@f$.
 * @f$R_n@f$ can be described as the difference between the number of points
 * in @f$[0, t)^s@f$ and the volume of that subregion. The
 * @f$\mathcal{L}_2@f$-<strong>star discrepancy</strong> is defined as the
 * @f$\mathcal{L}_2@f$ norm of @f$R_n(t, P)@f$, that is
 * @anchor REF_discrepancy_DiscL2Star_eq_discstar0
 * @f[
 *   \mathcal{D}_2^*(\mathcal{P}) = \left(\int_{[0, 1]^s} d^st\; R_n^2(t, P)\right)^{1/2} \tag{discstar0}
 * @f]
 * It can be calculated explicitly to give
 * @anchor REF_discrepancy_DiscL2Star_eq_discstar
 * @f[
 *   [\mathcal{D}_2^*(\mathcal{P})]^2 = \left(\frac{1}{3}\right)^s - \frac{2}{n} \sum_{i=1}^n \prod_{k=1}^s \left(\frac{1 - z_{ik}^2}{2}\right) + \frac{1}{n^2} \sum_{i=1}^n\sum_{j=1}^n \prod_{k=1}^s \left[1 - \max(z_{ik}, z_{jk})\right], \tag{discstar}
 * @f]
 * where @f$n@f$ is the number of points of @f$\mathcal{P}@f$, @f$s@f$ is the
 * dimension, and @f$z_{ik}@f$ is the @f$k@f$-th coordinate of point @f$i@f$.
 *
 * In one dimension, formula (
 * {@link REF_discrepancy_DiscL2Star_eq_discstar discstar} ) is
 * equivalent to
 * @anchor REF_discrepancy_DiscL2Star_eq_cramermises
 * @f[
 *   [\mathcal{D}_2^*(\mathcal{P})]^2 = \frac{1}{12n^2} + \frac{1}{n} \sum_{i=1}^n \left(z_{(i)} - \frac{i - 1/2}{n}\right)^2, \tag{cramermises}
 * @f]
 * where @f$n@f$ is the number of points of @f$\mathcal{P}@f$, and the
 * @f$z_{(i)}@f$ are the points of @f$\mathcal{P}@f$ *sorted* in increasing
 * order. This formula is faster to compute than the general formula (
 * {@link REF_discrepancy_DiscL2Star_eq_discstar discstar} ).
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class DiscL2Star extends Discrepancy {

   public double compute (double[][] points, int n, int s, double[] gamma) {
      return compute (points, n, s);
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in dimension
    * @f$s@f$. `points[i][j]` is the @f$j@f$-th coordinate of point
    * @f$i@f$. Both @f$i@f$ and @f$j@f$ start at 0.
    */
   public DiscL2Star (double[][] points, int n, int s) {
      super (points, n, s);
   }

   /**
    * Constructor with @f$n@f$ points in dimension @f$s@f$. The @f$n@f$
    * points will be chosen later.
    */
   public DiscL2Star (int n, int s) {
      super (null, n, s);
   }

   /**
    * Constructor with the point set `set`. All the points are copied in
    * an internal array.
    */
   public DiscL2Star (PointSet set) {
      super(set);
   }

   /**
    * Empty constructor. One *must* set the points and the dimension
    * before calling any method.
    */
   public DiscL2Star() {
      super();
   }

   /**
    * Computes the traditional @f$\mathcal{L}_2@f$-<strong>star
    * discrepancy</strong> (
    * {@link REF_discrepancy_DiscL2Star_eq_discstar
    * discstar} ) for the first @f$n@f$ points of `points`, in dimension
    * @f$s@f$.
    */
   public double compute (double[][] points, int n, int s) {
      double sum = 0.0;
      for (int i = 0; i < n; ++i) {
         double prod = 1.0;
         for (int j = 0; j < s; ++j) {
             double u = points[i][j];
             prod *= (1.0 - u) * (1.0 + u);
         }
         sum += prod;
      }
      double disc = -Math.pow(0.5, (double) (s-1))*sum/n;

      sum = 0.0;
      for (int i = 0; i < n; ++i) {
         double prod = 1.0;
         for (int j = 0; j < s; ++j)
             prod *= 1.0 - points[i][j];
         sum += prod;
      }

      double sum2 = 0.0;
      for (int i = 0; i < n - 1; ++i) {
         for (int j = i + 1; j < n; ++j) {
            double prod = 1.0;
            for (int k = 0; k < s; ++k)
                prod *= 1.0 - Math.max(points[i][k], points[j][k]);
            sum2 += prod;
         }
      }

      disc += (sum + 2.0*sum2)/((long)n*n);
      disc += Math.pow(1.0/3.0, s);
      if (disc < 0.0)
         return 0.0;
      return Math.sqrt(disc);
    }

   /**
    * Computes the traditional @f$\mathcal{L}_2@f$-<strong>star
    * discrepancy</strong> for the set of @f$n@f$ 1-dimensional points
    * @f$T@f$, using formula (
    * {@link REF_discrepancy_DiscL2Star_eq_cramermises
    * cramermises} ) above. The points of @f$T@f$ *must be sorted* before
    * calling this method.
    */
   public double compute (double[] T, int n) {
      double W2 = 0.0;
      double v;
      for (int i = 0; i < n; i++) {
         v = T[i] - (i + 0.5) / n;
         W2 += v * v;
      }
      W2 += 1.0 / (12.0 * n);

      return Math.sqrt(W2 / n);
    }

}