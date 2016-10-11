/*
 * Class:        DiscL2Unanchored
 * Description:  computes the L_2 unanchored discrepancy
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
 * A discrepancy is said to be reflection-invariant if it has the same value
 * when the points are reflected through any plane @f$x_j= 1/2@f$, passing
 * through the center of the unit hypercube, i.e. when any one of the
 * coordinates, say @f$z_j@f$, is replaced by @f$1 - z_j@f$ for all the
 * points. The star discrepancy is not reflection-invariant because it is
 * anchored at the origin, but the unanchored discrepancy is. This
 * discrepancy counts the points in all boxes @f$[x, y) \in[0,1)^s@f$.
 *
 * This class computes the @f$\mathcal{L}_2@f$-<strong>unanchored
 * discrepancy</strong> for a set of points @f$\mathcal{P}@f$ @cite vMOR94a,
 * @cite tHIC99a&thinsp;, given by
 * @anchor REF_discrepancy_DiscL2Unanchored_disc_unan
 * @f[
 *   [\mathcal{D}(\mathcal{P})]^2 = \left(\frac{1}{12}\right)^s - \frac{2}{n} \sum_{i=1}^n \prod_{k=1}^s \left[\frac{z_{ik}(1 - z_{ik})}{2}\right] + \frac{1}{n^2} \sum_{i=1}^n\sum_{j=1}^n \prod_{k=1}^s \left[\min(z_{ik}, z_{jk}) - z_{ik} z_{jk}\right], \tag{disc.unan}
 * @f]
 * where @f$n@f$ is the number of points of @f$\mathcal{P}@f$, @f$s@f$ is the
 * dimension, and @f$z_{ik}@f$ is the @f$k@f$-th coordinate of point @f$i@f$.
 *
 * In one dimension, formula (
 * {@link REF_discrepancy_DiscL2Unanchored_disc_unan disc.unan}
 * ) is equivalent to
 * @anchor REF_discrepancy_DiscL2Unanchored_disc_unanD1
 * @f[
 *   [\mathcal{D}(\mathcal{P})]^2 = \frac{1}{12} - \frac{1}{n} \sum_{i=1}^n {z_i(1 - z_i)} + \frac{1}{n^2} \sum_{i=1}^n\sum_{j=1}^n (\min(z_i, z_j) - z_i z_j), \tag{unanD1}
 * @f]
 * where @f$z_i@f$ is the point @f$i@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class DiscL2Unanchored extends Discrepancy {

   public double compute (double[][] points, int n, int s, double[] gamma) {
      return compute (points, n, s);
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in @f$s@f$
    * dimensions. `points[i][j]` is the @f$j@f$-th coordinate of point
    * @f$i@f$. Both @f$i@f$ and @f$j@f$ start at 0.
    */
   public DiscL2Unanchored (double[][] points, int n, int s) {
      super (points, n, s);
   }

   /**
    * Constructor with @f$n@f$ points in dimension @f$s@f$. The @f$n@f$
    * points will be chosen later.
    */
   public DiscL2Unanchored (int n, int s) {
      super (null, n, s);
   }

   /**
    * Constructor with the point set `set`. All the points are copied in
    * an internal array.
    */
   public DiscL2Unanchored (PointSet set) {
      super(set);
   }

   /**
    * Empty constructor. One *must* set the points and the dimension
    * before calling any method.
    */
   public DiscL2Unanchored() {
      super();
   }

   /**
    * Computes the @f$\mathcal{L}_2@f$-<strong>unanchored
    * discrepancy</strong> (
    * {@link REF_discrepancy_DiscL2Unanchored_disc_unan
    * disc.unan} ) for the set of @f$n@f$ @f$s@f$-dimensional points
    * `points`.
    */
   public double compute (double[][] points, int n, int s) {
      double sum = 0.0;
      for (int i = 0; i < n; ++i) {
         double prod = 1.0;
         for (int k = 0; k < s; ++k)
             prod *= points[i][k]*(1.0 - points[i][k]);
         sum += prod;
      }
      double disc = sum/n * (1.0/n - Math.pow(0.5, (double) (s-1)));

      sum = 0.0;
      for (int i = 0; i < n - 1; ++i) {
         for (int j = i + 1; j < n; ++j) {
            double prod = 1.0;
            for (int k = 0; k < s; ++k)
                prod *= Math.min(points[i][k], points[j][k]) -
                         points[i][k] * points[j][k];
            sum += prod;
         }
      }

      disc += 2.0*sum / ((long)n*n);
      disc += Math.pow(1.0/12.0, s);
      if (disc < 0.0)
         return -1.0;
      return Math.sqrt(disc);
    }

   /**
    * Computes the @f$\mathcal{L}_2@f$-<strong>unanchored
    * discrepancy</strong> for the 1-dimensional set of @f$n@f$ points
    * @f$T@f$, using formula (
    * {@link REF_discrepancy_DiscL2Unanchored_disc_unanD1
    * unanD1} ).
    */
   public double compute (double[] T, int n) {
      double sum = 0.0;
      for (int i = 0; i < n; ++i)
         sum += T[i]*(1.0 - T[i]);
      double disc = -(1.0 - 1.0/n) * sum / n;

      double sum2 = 0.0;
      for (int i = 0; i < n - 1; ++i) {
         for (int j = i + 1; j < n; ++j)
            sum2 += Math.min(T[i], T[j]) - T[i] * T[j];
      }

      disc += 2.0*sum2 / ((long)n*n);
      disc += 1.0/12.0;
      if (disc < 0.0)
         return -1.0;
      return Math.sqrt(disc);
    }

}