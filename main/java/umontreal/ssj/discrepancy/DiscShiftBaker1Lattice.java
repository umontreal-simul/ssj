/*
 * Class:        DiscShiftBaker1Lattice
 * Description:  computes a discrepancy for the randomly shifted, then baker
                 folded points of an integration lattice
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
 * This class computes the same **discrepancy** in @cite vHIC02a&thinsp; (eq.
 * 16) for the randomly shifted points of a set @f$\mathcal{L}@f$ as given in
 * eq. ( {@link REF_discrepancy_DiscShiftBaker1_baker1 baker1}
 * ) for class  @ref DiscShiftBaker1, but for the special case when the
 * points are the nodes of an integration lattice. It is given by
 * @anchor REF_discrepancy_DiscShiftBaker1Lattice_shiftBaker1lat
 * @f{align}{
 *    [\mathcal{D}(\mathcal{L})]^2 
 *    & 
 *   =
 *    -1 + \frac{1}{n} \sum_{i=1}^n \prod_{r=1}^s \left[1 - \frac{4\gamma_r^2}{3} \left[B_4(x_{ir}) - B_4(\{x_{ir}-1/2\})\right]\right. - \nonumber
 *    \\  &   
 *    \frac{\gamma_r^4}{9} \left[7B_4(x_{ir}) - 2B_4(\{x_{ir}-1/2\})\right] \left. {} - \frac{16\gamma_r^4}{45} \left[B_6(x_{ir}) -B_6(\{x_{ir}-1/2\})\right] \right], \tag{shiftBaker1lat}
 * @f}
 * where @f$n@f$ is the number of points of @f$\mathcal{L}@f$, @f$s@f$ is the
 * dimension of the points, @f$x_{ir}@f$ is the @f$r@f$-th coordinate of
 * point @f$i@f$, and the @f$\gamma_r@f$ are arbitrary positive weights. The
 * @f$B_{\alpha}(x)@f$ are the Bernoulli polynomials @cite mABR70a&thinsp;
 * (chap. 23) of degree @f$\alpha@f$. The discrepancy is much faster to
 * calculate for a lattice than for a general point set. For a
 * *1-dimensional* lattice, the discrepancy becomes
 * @anchor REF_discrepancy_DiscShiftBaker1Lattice_shiftBaker1latdim1
 * @f{align}{
 *    [\mathcal{D}(\mathcal{L})]^2 
 *    & 
 *   =
 *    - \frac{1}{n} \sum_{i=1}^n \left[\frac{4\gamma^2}{3} \left[B_4(x_i) - B_4(\{x_i-1/2\})\right]\right. + \nonumber
 *    \\  &   
 *    \frac{\gamma^4}{9} \left[7B_4(x_i) - 2B_4(\{x_i-1/2\})\right] \left. {} + \frac{16\gamma^4}{45} \left[B_6(x_i) -B_6(\{x_i-1/2\})\right] \right], \tag{shiftBaker1latdim1}
 * @f}
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class DiscShiftBaker1Lattice extends DiscShiftBaker1 {
   static final double TRENTEUN24 = 31.0/24.0;
   static final double SEPT24 = 7.0/24.0;

   private static double computeFactor (double x, double C1, double C2, double C3) {
      // Compute 1 factor of the product with Cr = C[r], r = 1,2,3, and x = Point[i][r]
      // This is a simplified form of computeFactor1.

      double pol1;   // Bernoulli(4, u) - Bernoulli(4, v)
      double pol2;   // 7*Bernoulli(4, u) - 2*Bernoulli(4, v)
      double pol3;   // Bernoulli(6, u) - Bernoulli(6, v)
      double temp;

      if (x >= 0.5) {
          // v = x - 0.5;
          pol1 = -0.5625 + x*(3.0 - x*(4.5 - 2.0*x));
          pol2 = -TRENTEUN24 + x*(6.0 - x*(4.0 + x*(6.0 - 5.0*x)));
          temp = 1.0 + x*(-6.0 + 4.0*x);
          pol3 = 0.046875 * temp * temp * (4.0*x - 3.0);
     } else {
          // v = x + 0.5;
          pol1 = -0.0625 + x*x*(1.5 - 2.0*x);
          pol2 = -SEPT24 + x*x*(8.0 - x*(14.0 - 5.0*x));
          temp = 1.0 + x*(2.0 - 4.0*x);
          pol3 = -0.046875 * temp * temp * (4.0*x - 1.0);
      }

      temp = C1 * pol1 + C2 * pol2 + C3 * pol3;
      return temp;
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in @f$s@f$
    * dimensions, and with all weights @f$\gamma_r =1@f$. `points[i][r]`
    * is the <tt>r</tt>-th coordinate of point `i`. Indices `i` and `r`
    * start at 0.
    */
   public DiscShiftBaker1Lattice (double[][] points, int n, int s) {
      super (points, n, s);
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in @f$s@f$
    * dimensions, with weights @f$\gamma_r = @f$ `gamma[r-1]`.
    * `points[i][r]` is the <tt>r</tt>-th coordinate of point `i`. Indices
    * `i` and `r` start at 0.
    */
   public DiscShiftBaker1Lattice (double[][] points, int n, int s,
                                  double[] gamma) {
      super (points, n, s, gamma);
   }

   /**
    * The number of points is @f$n@f$, the dimension @f$s@f$, and the
    * @f$s@f$ weight factors are <tt>gamma[</tt>@f$r@f$<tt>]</tt>, @f$r =
    * 0, 1, …, (s-1)@f$. The @f$n@f$ points will be chosen later.
    */
   public DiscShiftBaker1Lattice (int n, int s, double[] gamma) {
      super (n, s, gamma);
   }

   /**
    * Constructor with the point set `set`. All the points are copied in
    * an internal array.
    */
   public DiscShiftBaker1Lattice (Rank1Lattice set) {
      super(set);
   }

   /**
    * Empty constructor. The points and parameters *must* be defined
    * before calling methods of this class.
    */
   public DiscShiftBaker1Lattice() {
   }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShiftBaker1Lattice_shiftBaker1lat
    * shiftBaker1lat} ) for the @f$s@f$-dimensional points of lattice
    * `points`, containing @f$n@f$ points. All weights @f$\gamma_r =
    * 1@f$.
    */
   public double compute (double[][] points, int n, int s) {
      setONES (s);
      return compute (points, n, s, ONES);
   }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShiftBaker1Lattice_shiftBaker1lat
    * shiftBaker1lat} ) for the @f$s@f$-dimensional points of lattice
    * `points`, containing @f$n@f$ points, with weights @f$\gamma_r = @f$
    * `gamma[r-1]`.
    */
   public double compute (double[][] points, int n, int s, double[] gamma) {
      double[] C1 = new double[s];
      double[] C2 = new double[s];
      double[] C3 = new double[s];
      setC (C1, C2, C3, gamma, s);
      double prod, fact;
      int r;

      double sum = 0.0;
      for (int i = 0; i < n; ++i) {
         prod = 1.0;
         for (r = 0; r < s; ++r) {
            fact = computeFactor (points[i][r], C1[r], C2[r], C3[r]);
            prod *= 1.0 - fact;
         }
         sum += prod;
      }
      double disc = sum/n - 1.0;
      if (disc < 0.0)   // Lost all precision: result meaningless
         return -1.0;
      return Math.sqrt(disc);
    }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShiftBaker1Lattice_shiftBaker1latdim1
    * shiftBaker1latdim1} ) with weight @f$\gamma=1@f$ for the
    * 1-dimensional lattice of @f$n@f$ points @f$T@f$.
    */
   public double compute (double[] T, int n) {
      return compute (T, n, 1.);
   }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShiftBaker1Lattice_shiftBaker1latdim1
    * shiftBaker1latdim1} ) for the 1-dimensional lattice of @f$n@f$
    * points @f$T@f$, with weight @f$\gamma=@f$ `gamma`.
    */
   public double compute (double[] T, int n, double gamma) {
      double[] C = setC(gamma);
      double C1 = C[0];
      double C2 = C[1];
      double C3 = C[2];

      double sum = 0.0;
      for (int i = 0; i < n; ++i) {
         double h = T[i];
         double pol1 = ((h - 2.0) * h + 1.0)*h*h - UNTRENTE;    // Bernoulli(4,h)
         double v = h - 0.5;
         if (v < 0.0)
             v += 1.0;
         double pol2 = ((v - 2.0) * v + 1.0)*v*v - UNTRENTE;    // Bernoulli(4,v)
         double temp = C1 * (pol1 - pol2) + C2 * (7.0*pol1 - 2.0*pol2);

         // Bernoulli (6, h) and Bernoulli (6, v)
         pol1 = (((h - 3.0) * h + 2.5) * h*h - 0.5) * h*h + QUARAN;
         pol2 = (((v - 3.0) * v + 2.5) * v*v - 0.5) * v*v + QUARAN;
         temp += C3 * (pol1 - pol2);
         sum -= temp;
      }

      double disc = sum/n;
      if (disc < 0.0)
         return -1.0;
      return Math.sqrt(disc);
   }

}