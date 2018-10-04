/*
 * Class:        Searcher
 * Description:  searches the best lattices of rank 1 with respect to a 
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
   import umontreal.ssj.rng.*;
   import umontreal.ssj.hups.*;
   import umontreal.ssj.util.Num;

/**
 * This class implements methods to search for the best lattices of rank 1,
 * defined as follows @cite vSLO94a&thinsp;. Given a positive integer @f$n@f$
 * and a @f$s@f$-dimensional integer vector @f$(a_0,…,a_{s-1})@f$, where @f$0
 * \le a_j < n@f$ for each @f$j@f$, the points are defined by
 * @f[
 *   \mathbf{u}_i = (i/n)(a_0, a_1, …, a_{s-1}) \bmod1
 * @f]
 * for @f$i=0,…,n-1@f$. Here we always choose @f$a_0=1@f$.
 *
 * The discrepancy object in the `Searcher` constructor *must* fix the number
 * of points @f$n@f$, the maximal dimension @f$s@f$ of the lattice and
 * possibly, the weight factors @f$\gamma_j@f$. Then the search program will
 * examine different lattices with @f$n@f$, @f$s@f$ and @f$\gamma_j@f$ fixed
 * in order to find the best amongst those examined. Searches may be
 * exhaustive or random.
 *
 * One may consider only a subset of the possible lattices according to some
 * criterion. Some of these restricted searches are implemented in classes
 * derived from `Searcher`.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class Searcher {

   protected static RandomStream gen = new LFSR113(); // RNG used in random searches
   protected Discrepancy disc;
   protected PointSet lat;
   protected double[] gamma;
   protected double bestVal;       // best value of discr in dim s
   protected int[] bestAs;         // best generator vector for lattices
   protected boolean primeN = false;   // true if n is prime
   protected boolean power2F = false;  // true if n is a power of 2



   protected void print (int[] y, int s) {
      // For debugging
      System.out.printf ("  a = [ ");
      for (int j = 0; j < s; ++j) {
         System.out.printf ("%d  ", y[j]);
      }
      System.out.printf ("]%n");
   }


   private void incr (int[] y, int n, int s) {
      // increment y[i]; if y[i] == n, set y[i] = 1 and increment y[i-1]; recursively.
      // Thus we consider all y[i], each component taking all values in [1, 2, ..., n-1]
      // If y[1] >= n, stop. The exhaustive search of all y[j] is done.

      for (int i = s - 1; i >= 1; --i) {
         ++y[i];
         if (y[i] < n)
            return;
         if (i > 1)
            y[i] = 1;
      }
   }


   private void incrPrime (int[] y, int n, int s) {
      // increment y[i] as in incr; but keeps only y[i] relatively prime to n.

      for (int i = s - 1; i >= 1; --i) {
         ++y[i];
         if (power2F) {
            ++y[i];
         } else {
            while (Num.gcd(n, y[i]) != 1)
               ++y[i];
         }
         if (y[i] < n)
            return ;
         if (i > 1)
            y[i] = 1;
      }
   }


   private double exhaust (int s, boolean relPrime) {
      // If relPrime is true, only the multipliers y_j relatively prime to n
      // are considered; otherwise, all multipliers are.

      int n = disc.getNumPoints();
      gamma = disc.getGamma();
      double err;
      int j;
      int[] y = new int[s];
      for (j = 0; j < s; j++)
         y[j] = 1;
      bestVal = Double.MAX_VALUE;
      bestAs[0] = 1;

      while (y[1] < n) {
   //       print (y, s);
         lat = new Rank1Lattice(n, y, s);
         err = disc.compute (lat, gamma);
         // System.out.printf ("   disc = %8.4f%n", err);
         if (err < bestVal) {
            bestVal = err;
            // keep the best y in bestAs
            for (j = 1; j < s; j++)
               bestAs[j] = y[j];
         }
         if (relPrime)
            incrPrime (y, n, s);
         else
            incr (y, n, s);
      }

      return bestVal;
   }


   private double random (int s, int k, boolean relPrime) {
      // If relPrime is true, only the multipliers y_j relatively prime to n
      // are considered; otherwise, all multipliers are.

      int n = disc.getNumPoints();
      gamma = disc.getGamma();
      final int nm1 = n - 1;
      double err;
      int j;
      bestVal = Double.MAX_VALUE;
      int[] y = new int[s];
      y[0] = bestAs[0] = 1;

      for (int i = 0; i < k; i++) {
         for (j = 1; j < s; j++) {
            if (power2F) {
               y[j] = gen.nextInt (1, nm1);
               if (relPrime)
                  y[j] |= 1;
           } else {
               do {
                  y[j] = gen.nextInt (1, nm1);
               } while (relPrime && (Num.gcd(n, y[j]) != 1));
            }
         }
   //      print (y, s);
         lat = new Rank1Lattice(n, y, s);
         err = disc.compute (lat, gamma);
         if (err < bestVal) {
            bestVal = err;
            for (j = 1; j < s; j++)
               bestAs[j] = y[j];
         }
      }

      return bestVal;
   }

   /**
    * The number of points @f$n@f$, the dimension @f$s@f$, and possibly
    * the @f$s@f$ weight factors @f$\gamma_j@f$ must be given in `disc`.
    * The @f$n@f$ points of each lattice will be generated in the search.
    * The flag `primeN` indicates whether @f$n@f$ is a prime number
    * (<tt>true</tt>) or not (<tt>false</tt>). This is used in the
    * `*Prime` methods.
    */
   public Searcher (Discrepancy disc, boolean primeN) {
      this.disc = disc;
      int s = disc.getDimension();
      int n = disc.getNumPoints();
      bestAs = new int[s];
      this.primeN = primeN;
      if ((n & (n-1)) == 0)
         power2F = true;
      else
         power2F = false;
   }

   /**
    * Exhaustive search to find the lattice with the best (the smallest)
    * discrepancy in dimension @f$s@f$. The search runs over all values of
    * the generator @f$a_j = 1, 2, …, (n-1)@f$ and over all dimensions up
    * to @f$s@f$ inclusively. The first component @f$a_0@f$ of the
    * generator is always set to 1. The method returns the best value of
    * the discrepancy found in dimension @f$s@f$.
    */
   public double exhaust (int s) {
      return exhaust (s, false);
   }

   /**
    * Similar to  {@link #exhaust(int) exhaust(s)}, except that only
    * values of @f$a_j@f$ *relatively prime* to @f$n@f$ are considered.
    */
   public double exhaustPrime (int s) {
      if (primeN)
         return exhaust (s, false);   // all a_j are rel. prime to n
      else
         return exhaust (s, true);
   }

   /**
    * Random search to find the lattice with the best (the smallest)
    * discrepancy in dimension @f$s@f$. At most @f$k@f$ random
    * (@f$s@f$-dimensional) generators @f$\mathbf{a}@f$ are examined. Each
    * random component @f$a_j@f$ takes values over the integers @f$1, 2,
    * …, (n-1)@f$. The first component @f$a_0@f$ is always set to 1. The
    * method returns the best value of the discrepancy found in dimension
    * @f$s@f$.
    */
   public double random (int s, int k) {
      return random(s, k, false);
   }

   /**
    * Similar to  {@link #random(int,int) random(s, k)}, except that only
    * values of @f$a_j@f$ *relatively prime* to @f$n@f$ are considered.
    */
   public double randomPrime (int s, int k) {
      if (primeN)
         return random(s, k, false);
      else
         return random(s, k, true);
   }

   /**
    * Returns the best value of the discrepancy found in the last search.
    */
   public double getBestVal() {
      return bestVal;
   }

   /**
    * Returns the generator of the lattice which gave the best value of
    * the discrepancy in the last search. The components of this generator
    * are returned as @f$a_j@f$, for @f$j = 0, 1, …, (s-1)@f$.
    */
   public int[] getBestAs() {
      return bestAs;
   }

   /**
    * Initializes the random number generator used in random searches with
    * the starting seed `seed`. If this method is not called, a default
    * seed will be used.
    */
   public void initGen (int seed) {
      if (seed == 0 || seed == 1)
         seed = 7654321;
      final int s = 12345;
      int[] state = { seed, s, s, s };
      LFSR113.setPackageSeed(state);
      gen = new LFSR113();
   }

}