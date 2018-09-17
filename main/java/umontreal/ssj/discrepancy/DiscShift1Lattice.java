/*
 * Class:        DiscShift1Lattice
 * Description:  computes a discrepancy for the randomly shifted points of a
                 lattice
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
   import umontreal.ssj.hups.Rank1Lattice;

/**
 * This class computes the same **discrepancy** for randomly shifted points
 * of a set @f$\mathcal{L}@f$ as given in eq. (
 * {@link REF_discrepancy_DiscShift1_shift1 shift1} ) for class
 * @ref DiscShift1, but for the special case when the points are the nodes of
 * an integration lattice @cite vHIC02a&thinsp; (eq. 16). It is given by
 * @anchor REF_discrepancy_DiscShift1Lattice_shift1lat
 * @f[
 *   [\mathcal{D}(\mathcal{L})]^2 = -1 + \frac{1}{n} \sum_{i=1}^n \prod_{r=1}^s \left[1 + \gamma_r^2 B_2(x_{ir})\right], \tag{shift1lat}
 * @f]
 * where @f$n@f$ is the number of points of @f$\mathcal{L}@f$, @f$s@f$ is the
 * dimension of the points, @f$x_{ir}@f$ is the @f$r@f$-th coordinate of
 * point @f$i@f$, and the @f$\gamma_r@f$ are arbitrary positive weights.
 * @f$B_2(x)@f$ is the Bernoulli polynomials @cite mABR70a&thinsp; (chap. 23)
 * of degree @f$2@f$. For a *1-dimensional* lattice, the discrepancy becomes
 * @anchor REF_discrepancy_DiscShift1Lattice_shift1dim1lat
 * @f[
 *   [\mathcal{D}(\mathcal{L})]^2 = \frac{1}{n} \sum_{i=1}^n B_2(x_i). \tag{shift1dim1lat}
 * @f]
 * Computing the discrepancy for a lattice is much faster than for a general
 * point set.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class DiscShift1Lattice extends DiscShift1 {

   /**
    * Constructor with the @f$n@f$ points `points[i]` in @f$s@f$
    * dimensions with all weights @f$\gamma_r =1@f$. Element
    * `points[i][j]` is the <tt>j</tt>-th coordinate of point `i`. Indices
    * `i` and `j` start at 0.
    */
   public DiscShift1Lattice (double[][] points, int n, int s) {
      super (points, n, s);
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in @f$s@f$
    * dimensions with the weights @f$\gamma_r = @f$ `gamma[r-1]`, @f$r =
    * 1, …, s@f$. `points[i][j]` is the <tt>j</tt>-th coordinate of point
    * `i`. Indices `i` and `j` start at 0.
    */
   public DiscShift1Lattice (double[][] points, int n, int s, double[] gamma) {
      super (points, n, s, gamma);
   }

   /**
    * The number of points is @f$n@f$, the dimension @f$s@f$, and the
    * @f$s@f$ weight factors are <tt>gamma[</tt>@f$j@f$<tt>]</tt>, @f$j =
    * 0, 1, …, (s-1)@f$. The @f$n@f$ points will be chosen later.
    */
   public DiscShift1Lattice (int n, int s, double[] gamma) {
      super (n, s, gamma);
   }

   /**
    * Constructor with the lattice `set`. All the points are copied in an
    * internal array.
    */
   public DiscShift1Lattice (Rank1Lattice set) {
      super(set);
   }

   /**
    * Empty constructor. The points and parameters *must* be defined
    * before calling methods of this class.
    */
   public DiscShift1Lattice() {
   }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShift1Lattice_shift1lat
    * shift1lat} ) for the @f$s@f$-dimensional points of lattice `points`,
    * containing @f$n@f$ points. All weights @f$\gamma_r = 1@f$.
    */
   public double compute (double[][] points, int n, int s) {
      setONES (s);
      return compute (points, n, s, ONES);
   }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShift1Lattice_shift1lat
    * shift1lat} ) for the @f$s@f$-dimensional points of lattice `points`,
    * containing @f$n@f$ points, with weights @f$\gamma_r = @f$
    * `gamma[r-1]`.
    */
   public double compute (double[][] points, int n, int s, double[] gamma) {
      double[] C1 = new double[s];      // gamma_r^2
      setC (C1, gamma, s);

      double sum = 0.0;
      for (int i = 0; i < n; ++i) {
         double prod = 1.0;
         for (int r = 0; r < s; ++r) {
             double u = points[i][r];
             double pol1 = u*(u - 1.0) + UNSIX;   // Bernoulli(2, u)
             prod *= 1.0 + C1[r] * pol1;
         }
         sum += prod;
      }

      double disc = sum / n - 1.0;
      if (disc < 0.0)
         return -1.0;
      return Math.sqrt(disc);
    }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShift1Lattice_shift1dim1lat
    * shift1dim1lat} ) for the 1-dimensional lattice of @f$n@f$ points
    * @f$T@f$.
    */
   public double compute (double[] T, int n) {
      double sum = 0.0;
      for (int i = 0; i < n; ++i) {
         double h = T[i];
         double pol1 = h*(h - 1.0) + UNSIX;   // Bernoulli(2, h)
         sum += pol1;
      }

      double disc = sum / n;
      if (disc < 0.0)
         return -1.0;
      return Math.sqrt(disc);
    }

}