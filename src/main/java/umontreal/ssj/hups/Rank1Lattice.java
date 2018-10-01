/*
 * Class:        Rank1Lattice
 * Description:  Rank-1 lattice
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author
 * @since
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package umontreal.ssj.hups;
import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.rng.RandomStream;

/**
 * This class implements point sets specified by integration lattices of rank
 * 1. They are defined as follows @cite vSLO94a&thinsp;. One selects an
 * arbitrary positive integer @f$n@f$ and a @f$s@f$-dimensional integer
 * vector @f$(a_0,…,a_{s-1})@f$. [Usually, @f$a_0=1@f$ and @f$0 \le a_j <
 * n@f$ for each @f$j@f$; when the @f$a_j@f$ are outside the interval
 * @f$[0,n)@f$, then we replace @f$a_j@f$ by (@f$a_j \bmod n@f$) in all
 * calculations.] The points are defined by
 * @f[
 *   \mathbf{u}_i = (i/n)(a_0, a_1, …, a_{s-1}) \bmod1
 * @f]
 * for @f$i=0,…,n-1@f$. These @f$n@f$ points are distinct provided that
 * @f$n@f$ and the @f$a_j@f$’s have no common factor.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class Rank1Lattice extends PointSet {

   protected int[] genAs;          // Lattice generator:  a[i]
   protected double[] v;           // Lattice vector:  v[i] = a[i]/n
   protected double normFactor;    // 1/n.

   private void initN (int n) {
      numPoints = n;
      normFactor = 1.0 / (double) n;
      for (int j = 0; j < dim; j++) {
         int amod = (genAs[j] % n) + (genAs[j] < 0 ? n : 0);
         v[j] = normFactor * amod;
      }
   }

   /**
    * Instantiates a  @ref Rank1Lattice with @f$n@f$ points and lattice
    * vector @f$a@f$ of dimension @f$s@f$.
    *  @param n            there are n points
    *  @param a            the lattice vector
    *  @param s            dimension of the lattice vector a
    */
   public Rank1Lattice (int n, int[] a, int s) {
      dim = s;
      v = new double[s];
      genAs = new int[s];
      for (int j = 0; j < s; j++) {
         genAs[j] = a[j];
      }
      initN (n);
   }

   /**
    * Resets the number of points of the lattice to @f$n@f$. The dimension
    * @f$s@f$ and the @f$a_j@f$ are unchanged.
    */
   public void setNumPoints (int n) {
      initN(n);
   }

   /**
    * Returns the generator @f$a_j@f$ of the lattice. (The original ones
    * before they are reset to @f$a_j \bmod n@f$). Its components are
    * returned as <tt>a[</tt>@f$j@f$<tt>]</tt>, for @f$j = 0, 1, …,
    * (s-1)@f$.
    */
   public int[] getAs() {
      return genAs;
   }

   /**
    * Adds a random shift to all the points of the point set, using stream
    * `stream` to generate the random numbers. For each coordinate @f$j@f$
    * from `d1` to `d2-1`, the shift @f$d_j@f$ is generated uniformly over
    * @f$[0, 1)@f$ and added modulo @f$1@f$ to all the coordinates of all
    * the points.
    *  @param d1           lower dimension of shift
    *  @param d2           upper dimension of shift is d2 - 1
    *  @param stream       random number stream used to generate uniforms
    */
   public void addRandomShift (int d1, int d2, RandomStream stream) {
      if (null == stream)
         throw new IllegalArgumentException (
              PrintfFormat.NEWLINE +
                  "   Calling addRandomShift with null stream");
      if (0 == d2)
         d2 = Math.max (1, dim);
      if (shift == null) {
         shift = new double[d2];
         capacityShift = d2;
      } else if (d2 > capacityShift) {
         int d3 = Math.max (4, capacityShift);
         while (d2 > d3)
            d3 *= 2;
         double[] temp = new double[d3];
         capacityShift = d3;
         for (int i = 0; i < d1; i++)
            temp[i] = shift[i];
         shift = temp;
      }
      dimShift = d2;
      for (int i = d1; i < d2; i++)
         shift[i] = stream.nextDouble ();
      shiftStream = stream;
   }

   /**
    * Clears the random shift.
    */
   public void clearRandomShift() {
      super.clearRandomShift();
      shift = null;
   }


   public String toString() {
      StringBuffer sb = new StringBuffer ("Rank1Lattice:" +
                                           PrintfFormat.NEWLINE);
      sb.append (super.toString());
      return sb.toString();
   }


   public double getCoordinate (int i, int j) {
      double x = (v[j] * i) % 1.0;
      if (shift != null) {
         if (j >= dimShift)   // Extend the shift.
            addRandomShift (dimShift, j + 1, shiftStream);
         x += shift[j];
         if (x >= 1.0)
            x -= 1.0;
         if (x <= 0.0)
            x = EpsilonHalf;  // avoid x = 0
       }
      return x;
   }


   // Recursive method that computes a^e mod m.
   protected long modPower (long a, int e, int m) {
      // If parameters a and m == numPoints could be omitted, then
      // the routine would run much faster due to reduced stack usage.
      // Note that a can be larger than m, e.g. in lattice sequences !

      if (e == 0)
         return 1;
      else if (e == 1)
         return a % m;
      else if ((e & 1) == 1)
         return (a * modPower(a, e - 1, m)) % m;
      else {
         long p = modPower(a, e / 2, m);
         return (p * p) % m;
      }
   }

   protected double radicalInverse (int base, int i) {
      double digit, radical, inverse;
      digit = radical = 1.0 / (double) base;
      for (inverse = 0.0; i > 0; i /= base) {
         inverse += digit * (double) (i % base);
         digit *= radical;
      }
      return inverse;
   }

   public PointSetIterator iterator() {
      return new Rank1LatticeIterator();
   }

// ************************************************************************

   protected class Rank1LatticeIterator extends PointSet.DefaultPointSetIterator
   {
      public double nextCoordinate() {
         // I tried with long's and with double's. The double version is
         // 4.5 times faster than the long version.
         if (curPointIndex >= numPoints || curCoordIndex >= dim)
            outOfBounds();
//      return (curPointIndex * v[curCoordIndex++]) % 1.0;
         double x = (curPointIndex * v[curCoordIndex]) % 1.0;
         if (shift != null) {
             if (curCoordIndex >= dimShift)   // Extend the shift.
                addRandomShift (dimShift, curCoordIndex + 1, shiftStream);
             x += shift[curCoordIndex];
             if (x >= 1.0)
                x -= 1.0;
             if (x <= 0.0)
                x = EpsilonHalf;  // avoid x = 0
         }
         curCoordIndex++;
         return x;
      }
   }
}