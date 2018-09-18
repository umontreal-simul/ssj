/*
 * Class:        DiscL2Symmetric
 * Description:  computes the L_2 symmetric discrepancy
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
 * COMPLÉTER LA DOC ICI.
 *
 * This class computes the @f$\mathcal{L}_2@f$-<strong>symmetric
 * discrepancy</strong> for a set of points @cite vHIC98a,
 * @cite rLIA00a&thinsp;, given by
 * @anchor REF_discrepancy_DiscL2Symmetric_disc_sym
 * @f[
 *   [\mathcal{D}_s(\mathcal{P})]^2 = \left(\frac{4}{3}\right)^s - \frac{2}{n} \sum_{i=1}^n \prod_{k=1}^s \left(1 + 2z_{ik} - 2z_{ik}^2\right) + \frac{2^s}{n^2} \sum_{i=1}^n\sum_{j=1}^n \prod_{k=1}^s \left(1 - |z_{ik} - z_{jk}|\right), \tag{sym}
 * @f]
 * where @f$n@f$ is the number of points of @f$\mathcal{P}@f$, @f$s@f$ is the
 * dimension, and @f$z_{ik}@f$ is the @f$k@f$-th coordinate of point @f$i@f$.
 *
 * In one dimension, formula (
 * {@link REF_discrepancy_DiscL2Symmetric_disc_sym sym} ) is
 * equivalent to
 * @anchor REF_discrepancy_DiscL2Symmetric_disc_symD1
 * @f[
 *   [\mathcal{D}_s(\mathcal{P})]^2 = \frac{4}{3} - \frac{4}{n} \sum_{i=1}^n z_i(1 - z_i) - \frac{2}{n^2} \sum_{i=1}^n\sum_{j=1}^n|z_i - z_j|, \tag{symD1}
 * @f]
 * where @f$z_i@f$ is the coordinate of point @f$i@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class DiscL2Symmetric extends Discrepancy {

   public double compute (double[][] points, int n, int s, double[] gamma) {
      return compute (points, n, s);
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in @f$s@f$
    * dimensions. `points[i][j]` is the @f$j@f$-th coordinate of point
    * @f$i@f$. Both @f$i@f$ and @f$j@f$ start at 0.
    */
   public DiscL2Symmetric (double[][] points, int n, int s) {
      super (points, n, s);
   }

   /**
    * Constructor with @f$n@f$ points in dimension @f$s@f$. The @f$n@f$
    * points will be chosen later.
    */
   public DiscL2Symmetric (int n, int s) {
      super (null, n, s);
   }

   /**
    * Constructor with the point set `set`. All the points are copied in
    * an internal array.
    */
   public DiscL2Symmetric (PointSet set) {
      super(set);
   }

   /**
    * Empty constructor. One *must* set the points and the dimension
    * before calling any method.
    */
   public DiscL2Symmetric() {
      super();
   }

   /**
    * Computes the @f$\mathcal{L}_2@f$-<strong>symmetric
    * discrepancy</strong> for the set of @f$n@f$ @f$s@f$-dimensional
    * points `points`, using formula (
    * {@link REF_discrepancy_DiscL2Symmetric_disc_sym sym}
    * ).
    */
   public double compute (double[][] points, int n, int s) {
      double sum = 0.0;
      for (int i = 0; i < n; ++i) {
         double prod = 1.0;
         for (int j = 0; j < s; ++j) {
             double u = 0.5 - points[i][j];
             prod *= 1.5 - 2.0*u*u;
         }
         sum += prod;
      }
      double disc = -2.0*sum/n;

      sum = n;
      double sum2 = 0.0;
      for (int i = 0; i < n - 1; ++i) {
         for (int j = i + 1; j < n; ++j) {
            double prod = 1.0;
            for (int k = 0; k < s; ++k)
                prod *= 1.0 - Math.abs(points[i][k] - points[j][k]);
            sum2 += prod;
         }
      }

      disc += (sum + 2.0*sum2)* Math.pow(2.0, s)/((long)n*n);
      disc += Math.pow(4.0/3.0, s);
      if (disc < 0.0)
         return -1.0;
      return Math.sqrt(disc);
    }

   /**
    * Computes the @f$\mathcal{L}_2@f$-<strong>symmetric
    * discrepancy</strong> for the set of @f$n@f$ 1-dimensional points
    * @f$T@f$, using formula (
    * {@link REF_discrepancy_DiscL2Symmetric_disc_symD1
    * symD1} ).
    */
   public double compute (double[] T, int n) {
      double sum = 0.0;
      for (int i = 0; i < n; ++i)
         sum += T[i]*(1.0 - T[i]);
      double disc = -4.0*sum/n;

      sum = 0.0;
      for (int i = 0; i < n - 1; ++i) {
         for (int j = i + 1; j < n; ++j)
            sum += Math.abs(T[i] - T[j]);
      }

      disc -= 4.0*sum / ((long)n*n);
      disc += 4.0/3.0;
      if (disc < 0.0)
         return 0.0;
      return Math.sqrt(disc);
    }

}