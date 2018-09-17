/*
 * Class:        BigDiscShiftBaker1
 * Description:  computes a discrepancy for a randomly shifted, then baker
                 folded point set using multi-precision real numbers
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
import java.math.*;

/**
 * This class computes the same discrepancy as in
 * @ref umontreal.ssj.discrepancy.DiscShiftBaker1 [see eq. (
 * {@link REF_discrepancy_DiscShiftBaker1_baker1 baker1} )],
 * but uses multi-precision real numbers.
 *
 * **No method are implemented yet for this class.**
 *
 * <strong>Only the subclass
 * @ref umontreal.ssj.discrepancy.BigDiscShiftBaker1Lattice has implemented
 * methods to compute the discrepancy.</strong>
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class BigDiscShiftBaker1 extends BigDiscrepancy {
   static final BigDecimal BIG4s3 =
        new BigDecimal("1.3333333333333333333333333333333333");     // 4/3
   static final BigDecimal BIG1s9 =
        new BigDecimal("0.1111111111111111111111111111111111");     // 1/9
   static final BigDecimal BIG16s45 =
        new BigDecimal("0.35555555555555555555555555555555556");    // 16/45


   protected void setCBig (double[] gam, int s) {
      for (int i = 0; i < s; i++) {
         BigDecimal v = new BigDecimal(gam[i]);
         v = v.pow(2);                      // v = gam[i]*gam[i]
         C1Big[i] = v.multiply(BIG4s3);     // C1[i] = v * 4 / 3
         v = v.pow(2);                      // v = v*v
         C2Big[i] = v.multiply(BIG1s9);     // C2[i] = v / 9.0;
         C3Big[i] = v.multiply(BIG16s45);   // C3[i] = v * 16.0 / 45.0;
      }
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in @f$s@f$
    * dimensions, with all the weights @f$\gamma_r = 1@f$. `points[i][r]`
    * is the <tt>r</tt>-th coordinate of point `i`. Indices `i` and `r`
    * start at 0.
    */
   public BigDiscShiftBaker1 (double[][] points, int n, int s) {
      super (points, n, s);
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in @f$s@f$
    * dimensions, with weights @f$\gamma_r = @f$ `gamma[r-1]`.
    * `points[i][r]` is the <tt>r</tt>-th coordinate of point `i`. Indices
    * `i` and `r` start at 0.
    */
   public BigDiscShiftBaker1 (double[][] points, int n, int s,
                              double[] gamma) {
      super (points, n, s, gamma);
   }

   /**
    * Constructor for a lattice of @f$n@f$ points in @f$s@f$ dimensions,
    * with weights @f$\gamma_r = @f$ `gamma[r-1]`, @f$r = 1,2, …, s@f$.
    * The @f$n@f$ points will be chosen later.
    */
   public BigDiscShiftBaker1 (int n, int s, double[] gamma) {
      super (n, s, gamma);
   }

   /**
    * Constructor with the point set `set`. All the points are copied in
    * an internal array.
    */
   public BigDiscShiftBaker1 (PointSet set) {
      super(set);
   }

   /**
    * Empty constructor. One *must* set the points, the dimension, and the
    * weight factors before calling any method.
    */
   public BigDiscShiftBaker1() {
   }

   /**
    * NOT IMPLEMENTED.
    */
   public double compute (double[][] points, int n, int s) {
      setONES (s);
      return compute (points, n, s, ONES);
   }

   /**
    * NOT IMPLEMENTED.
    */
   public double compute (double[][] points, int n, int s, double[] gamma) {
      throw new UnsupportedOperationException ("method NOT IMPLEMENTED");
  }

}