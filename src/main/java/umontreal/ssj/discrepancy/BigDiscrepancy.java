/*
 * Class:        BigDiscrepancy
 * Description:  Base class of all discrepancy classes programmed with
                 multi-precision floating-point numbers
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
   import java.math.*;

/**
 * This *abstract* class is the base class of all discrepancy classes
 * programmed with floating-point numbers with multi-precision. For @f$n@f$
 * large, computing the discrepancy suffers from subtractive cancellation and
 * loses all precision if one uses `double` numbers. Using multi-precision
 * numbers allow us to compute the discrepancy for larger @f$n@f$, but the
 * computation is very slow.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public abstract class BigDiscrepancy extends Discrepancy {

   protected BigDecimal[] C1Big;    // functions of gamma[i]
   protected BigDecimal[] C2Big;    // functions of gamma[i]
   protected BigDecimal[] C3Big;    // functions of gamma[i]
   protected BigDecimal[] UBig;     // i/n
   protected BigDecimal[][] FactorBig;

   protected void reserveCBig (int s) {
      // Reserve memory for gamma factors
      C1Big = new BigDecimal[s];
      C2Big = new BigDecimal[s];
      C3Big = new BigDecimal[s];

      for (int i = 0; i < s; i++) {
         C1Big[i] = new BigDecimal(0);
         C2Big[i] = new BigDecimal(0);
         C3Big[i] = new BigDecimal(0);
      }
   }


   protected void reserveFactorBig (int n, int s) {
      // Reserve memory for factors
      FactorBig = new BigDecimal[n][s];
      for (int i = 0; i < n; i++) {
         FactorBig[i] = new BigDecimal[s];
         for (int j = 0; j < s; j++)
            FactorBig[i][j] = new BigDecimal(0);
      }
   }


   protected void setUBig (int n) {
      // Precompute all U[i] = i/n, with i = 0, 1, 2, ..., n-1
      UBig = new BigDecimal[n];
      BigDecimal Ninv = new BigDecimal(1);     // = 1/n
      Ninv = Ninv.divide(new BigDecimal(n), MathContext.DECIMAL128);

      for (int i = 0; i < n; i++) {
         UBig[i] = new BigDecimal(i);
         UBig[i] = UBig[i].multiply(Ninv);
      }
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in @f$s@f$
    * dimensions. `points[i][j]` is the @f$j@f$-th coordinate of point
    * @f$i@f$. Both @f$i@f$ and @f$j@f$ start at 0.
    */
   public BigDiscrepancy (double[][] points, int n, int s) {
      super(points, n, s);
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in @f$s@f$
    * dimensions with weight factors `gamma`. `points[i][j]` is the
    * @f$j@f$-th coordinate of point @f$i@f$. Both @f$i@f$ and @f$j@f$
    * start at 0.
    */
   public BigDiscrepancy (double[][] points, int n, int s, double[] gamma) {
      super(points, n, s, gamma);
   }

   /**
    * The number of points is @f$n@f$, the dimension @f$s@f$, and the
    * @f$s@f$ weight factors are <tt>gamma[</tt>@f$j@f$<tt>]</tt>, @f$j =
    * 0, 1, …, (s-1)@f$. The @f$n@f$ points will be chosen later.
    */
   public BigDiscrepancy (int n, int s, double[] gamma) {
      super(n, s, gamma);
   }

   /**
    * Constructor with the point set `set`. All the points are copied in
    * an internal array.
    */
   public BigDiscrepancy (PointSet set) {
      super(set);
   }

   /**
    * Empty constructor. The points and parameters *must* be defined
    * before calling methods of this or derived classes.
    */
   public BigDiscrepancy() {
   }

}