/*
 * Class:        BigDiscShiftBaker1Lattice
 * Description:  computes a discrepancy for a randomly shifted, then baker
                 folded points of a lattice using multi-precision real numbers
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
package umontreal.ssj.discrepancy; import java.math.*;

/**
 * This class computes the same discrepancy as in
 * @ref umontreal.ssj.discrepancy.DiscShiftBaker1Lattice [see eq. (
 * {@link REF_discrepancy_DiscShiftBaker1Lattice_shiftBaker1lat
 * shiftBaker1lat} )], but uses multi-precision real numbers.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class BigDiscShiftBaker1Lattice extends BigDiscShiftBaker1 {
   static final double TRENTEUN24 = 31.0/24.0;
   static final double SEPT24 = 7.0/24.0;

   static final BigDecimal BIGp5625 = new BigDecimal(0.5625);
   static final BigDecimal BIGp046875 = new BigDecimal(0.046875);
   static final BigDecimal BIGp0625 = new BigDecimal(0.0625);
   static final BigDecimal BIGp5 = new BigDecimal(0.5);
   static final BigDecimal BIG1p5 = new BigDecimal(1.5);
   static final BigDecimal BIG2 = new BigDecimal(2.0);
   static final BigDecimal BIG3 = new BigDecimal(3.0);
   static final BigDecimal BIG4 = new BigDecimal(4.0);
   static final BigDecimal BIG4p5 = new BigDecimal(4.5);
   static final BigDecimal BIG5 = new BigDecimal(5.0);
   static final BigDecimal BIG6 = new BigDecimal(6.0);
   static final BigDecimal BIG8 = new BigDecimal(8.0);
   static final BigDecimal BIG14 = new BigDecimal(14.0);
   static final BigDecimal BIG31s24 =
      new BigDecimal("1.2916666666666666666666666666666667");
   static final BigDecimal BIG7s24 =
      new BigDecimal("0.2916666666666666666666666666666667");

   /**
    * Constructor for a lattice of @f$n@f$ points in at most @f$s@f$
    * dimensions, with weights @f$\gamma_r = @f$ `gamma[r-1]`, @f$r =
    * 1,2, …, s@f$.
    */
   public BigDiscShiftBaker1Lattice (int n, int s, double[] gamma) {
      super (n, s, gamma);
      setUBig (n);
      reserveCBig (s);
      setCBig (this.gamma, s);
      reserveFactorBig (n, s);
      setFactorBig(n, s);
   }

   /**
    * Constructor for a lattice of @f$n@f$ points in at most @f$s@f$
    * dimensions, with weights @f$\gamma_r = 1@f$.
    */
   public BigDiscShiftBaker1Lattice (int n, int s) {
      this (n, s, null);
   }


   protected void setFactorBig (int n, int s) {
      // Precompute the factors in the discrepancy for each i/n
      for (int i = 0; i < n; i++) {
         for (int j = 0; j < s; j++) {
            FactorBig[i][j] = computeFactor (UBig[i], C1Big[j], C2Big[j],
                                             C3Big[j]);
         }
      }
   }


   private static BigDecimal computeFactor (BigDecimal x, BigDecimal C1,
                                            BigDecimal C2, BigDecimal C3) {
      // This is a simplified form of computeFactor(double, ...).
      // Compute 1 factor of the product with Cr = C[r], r = 1,2,3, and
      //    x = Point[i][r]

      BigDecimal pol1 = new BigDecimal(0); // Bernoulli(4,u) - Bernoulli(4,v)
      BigDecimal pol2 = new BigDecimal(0); // 7*Bernoulli(4,u) - 2*Bernoulli(4,v)
      BigDecimal pol3 = new BigDecimal(0); // Bernoulli(6,u) - Bernoulli(6,v)
      BigDecimal temp = new BigDecimal(0);
      pol1 = x;
      pol2 = x;
      pol3 = x;
      temp = x;

      if (x.compareTo(BIGp5) >= 0) {                    //  case v = x - 0.5;
          // pol1 = -0.5625 + x*(3.0 - x*(4.5 - 2.0*x));
          pol1 = pol1.multiply(BIG2);
          pol1 = pol1.subtract(BIG4p5);
          pol1 = pol1.multiply(x);
          pol1 = pol1.add(BIG3);
          pol1 = pol1.multiply(x);
          pol1 = pol1.subtract(BIGp5625);

          // pol2 = -31/24 + x*(6.0 - x*(4.0 + x*(6.0 - 5.0*x)));
          pol2 = pol2.multiply(BIG5);
          pol2 = pol2.negate();
          pol2 = pol2.add(BIG6);
          pol2 = pol2.multiply(x);
          pol2 = pol2.add(BIG4);
          pol2 = pol2.multiply(x);
          pol2 = pol2.negate();
          pol2 = pol2.add(BIG6);
          pol2 = pol2.multiply(x);
          pol2 = pol2.subtract(BIG31s24);

          // temp = 1.0 + x*(-6.0 + 4.0*x);
          temp = temp.multiply(BIG4);
          temp = temp.subtract(BIG6);
          temp = temp.multiply(x);
          temp = temp.add(BigDecimal.ONE);

          // pol3 = 0.046875 * temp * temp * (4.0*x - 3.0);
          pol3 = pol3.multiply(BIG4);
          pol3 = pol3.subtract(BIG3);
          pol3 = pol3.multiply(temp);
          pol3 = pol3.multiply(temp);
          pol3 = pol3.multiply(BIGp046875);

     } else {                                           //  case  v = x + 0.5;
          // pol1 = -0.0625 + x*x*(1.5 - 2.0*x);
          pol1 = pol1.multiply(BIG2);
          pol1 = pol1.negate();
          pol1 = pol1.add(BIG1p5);
          pol1 = pol1.multiply(x);
          pol1 = pol1.multiply(x);
          pol1 = pol1.subtract(BIGp0625);

          // pol2 = -7/24 + x*x*(8.0 - x*(14.0 - 5.0*x));
          pol2 = pol2.multiply(BIG5);
          pol2 = pol2.subtract(BIG14);
          pol2 = pol2.multiply(x);
          pol2 = pol2.add(BIG8);
          pol2 = pol2.multiply(x);
          pol2 = pol2.multiply(x);
          pol2 = pol2.subtract(BIG7s24);

          // temp = 1.0 + x*(2.0 - 4.0*x);
          temp = temp.multiply(BIG4);
          temp = temp.negate();
          temp = temp.add(BIG2);
          temp = temp.multiply(x);
          temp = temp.add(BigDecimal.ONE);

          // pol3 = -0.046875 * temp * temp * (4.0*x - 1.0);
          pol3 = pol3.multiply(BIG4);
          pol3 = pol3.subtract(BigDecimal.ONE);
          pol3 = pol3.multiply(temp);
          pol3 = pol3.multiply(temp);
          pol3 = pol3.multiply(BIGp046875);
          pol3 = pol3.negate();
       }

     // sum = C1 * pol1 + C2 * pol2 + C3 * pol3;
     BigDecimal sum = new BigDecimal(0);
     pol1 = pol1.multiply(C1);
     sum = sum.add(pol1);
     pol2 = pol2.multiply(C2);
     sum = sum.add(pol2);
     pol3 = pol3.multiply(C3);
     sum = sum.add(pol3);
     return sum;
   }

/**
 * Computes the discrepancy (
 * {@link REF_discrepancy_DiscShiftBaker1Lattice_shiftBaker1lat
 * shiftBaker1lat} ) for a rank-1 lattice in dimension @f$s@f$. This lattice
 * is generated by @f$a_r =@f$ `a[r-1]`, @f$r=1,2, …, s@f$.
 */
public double compute (long[] a, int s) {
      BigDecimal prod = new BigDecimal(1);
      BigDecimal tem = new BigDecimal(0);
      BigDecimal sum = new BigDecimal(0);
      int j, r;
      long nl = numPoints;

      for (long i = 0; i < nl; ++i) {
         prod = BigDecimal.ONE;
         for (j = 0; j < s; ++j) {
            r = (int)((i * a[j]) % nl);
            tem = FactorBig[r][j];
            tem = tem.negate();
            tem = tem.add(BigDecimal.ONE);
            prod = prod.multiply(tem);
         }
         sum = sum.add(prod);
      }

      sum = sum.divide(new BigDecimal(nl), MathContext.DECIMAL128);
      sum = sum.subtract(BigDecimal.ONE);
      double disc = sum.doubleValue();
      if (disc < 0.0)   // Lost all precision
         return -1.0;
      return Math.sqrt(disc);
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