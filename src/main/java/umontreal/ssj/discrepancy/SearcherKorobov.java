/*
 * Class:        SearcherKorobov
 * Description:  searches the best Korobov lattices with respect to a 
                 given discrepancy
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
   import umontreal.ssj.util.Num;
   import umontreal.ssj.hups.KorobovLattice;

/**
 * This class implements searches to find the best *Korobov* lattices with
 * respect to a given discrepancy. Given a positive integer @f$n@f$ and a
 * multiplier @f$a@f$, the points are defined by
 * @f[
 *   \mathbf{u}_i = (i/n)(1, a, a^2, …, a^{s-1}) \bmod1
 * @f]
 * for @f$i=0,…,n-1@f$.
 *
 * The discrepancy object in the constructor *must* fix the number of points
 * @f$n@f$, the maximal dimension @f$s@f$ of the lattice and possibly, the
 * weight factors @f$\gamma_j@f$. Then the search program will examine
 * different lattices with @f$n@f$, @f$s@f$ and @f$\gamma_j@f$ fixed in
 * order to find the best amongst those examined.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class SearcherKorobov extends Searcher {
   protected int bestA;           // best generator a for Korobov lattices


   protected void print (int y) {
      // For debugging
      System.out.printf ("  a = %d%n", y);
   }


   private void calcAs(int n, int s, int a) {
      long b = a;
      bestAs[0] = 1;
      for (int j = 1; j < s; ++j)
         bestAs[j] = (int)((b * bestAs[j - 1]) % n);
   }


   private double exhaust (int s, boolean relPrime) {
      int n = disc.getNumPoints();
      double err;
      gamma = disc.getGamma();
      bestVal = Double.MAX_VALUE;
      bestA = -1;

      for (int i = 2; i < n; i++) {
         if (relPrime) {
            // Consider only values i relatively prime to n
            if (power2F) {
               if ((i & 1) == 0)
                  continue;
            } else {
               if (Num.gcd(n, i) != 1)
                  continue;
            }
         }
      //   print (i);
         lat = new KorobovLattice(n, i, s);
         err = disc.compute (lat, gamma);
         if (err < bestVal) {
            bestVal = err;
            bestA = i;
         }
      }

      calcAs(n, s, bestA);
      return bestVal;
   }


   private double random (int s, int k, boolean relPrime) {
      int n = disc.getNumPoints();
      if (k >= n)
         return exhaust (s, relPrime);

      gamma = disc.getGamma();
      final int nm1 = n - 1;
      double err;
      bestVal = Double.MAX_VALUE;
      bestA = -1;
      int a;
      int i = 0;

      while (i < k) {
         if (power2F) {
            a = gen.nextInt (2, nm1);
            if (relPrime)
               a |= 1;
         } else {
            do {
               a = gen.nextInt (2, nm1);
            } while (relPrime && (Num.gcd(n, a) != 1));
         }
     //    print (a);
         lat = new KorobovLattice(n, a, s);
         i++;
         err = disc.compute (lat, gamma);
         if (err < bestVal) {
            bestVal = err;
            bestA = a;
         }
      }

      calcAs(n, s, bestA);
      return bestVal;
   }

   /**
    * The number of points @f$n@f$, the dimension @f$s@f$, and possibly
    * the @f$s@f$ weight factors @f$\gamma_j@f$ must be given in `disc`.
    * The @f$n@f$ points of the lattice will be generated in the search.
    * The flag `primeN` indicates whether @f$n@f$ is a prime number
    * (<tt>true</tt>) or not (<tt>false</tt>). This is used in the
    * `*Prime` methods.
    */
   public SearcherKorobov (Discrepancy disc, boolean primeN) {
      super(disc, primeN);
  }

   /**
    * Exhaustive search to find the best Korobov lattice (i.e. with the
    * smallest discrepancy) in dimension @f$s@f$. The search runs
    * exhaustively over all values @f$ a = 2, 3, …, (n-1)@f$, where
    * @f$a@f$ is the generator of the lattice. The method returns the best
    * value found for the discrepancy.
    */
   public double exhaust (int s) {
      return exhaust (s, false);
   }

   /**
    * Similar to  {@link #exhaust(int) exhaust(s)}, except that only
    * values of @f$a@f$ *relatively prime* to @f$n@f$ are considered.
    */
   public double exhaustPrime (int s) {
      if (primeN)
         return exhaust (s, false);   // all a are rel. prime to n
      else
         return exhaust (s, true);
   }

   /**
    * Random search to find the best Korobov lattice (with the smallest
    * discrepancy) in dimension @f$s@f$. @f$k@f$ random values @f$a@f$ are
    * examined as generator of the lattice. The @f$a@f$ takes values over
    * the integers @f$ a = 2, 3, …, (n-1)@f$. The method returns the best
    * value found for the discrepancy.
    */
   public double random (int s, int k) {
      return random(s, k, false);
   }

   /**
    * Similar to  {@link #random(int,int) random(s, k)}, except that only
    * values of @f$a@f$ *relatively prime* to @f$n@f$ are considered.
    */
   public double randomPrime (int s, int k) {
      if (primeN)
         return random(s, k, false);
      else
         return random(s, k, true);
   }

   /**
    * Returns the generator @f$a@f$ of the lattice which gave the best
    * value of the discrepancy in the last search.
    */
   public int getBestA() {
      return bestA;
   }

}