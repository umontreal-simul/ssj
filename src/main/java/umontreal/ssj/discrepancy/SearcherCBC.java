/*
 * Class:        SearcherCBC
 * Description:  searches the best rank-1 lattices with respect to a given
                 discrepancy, using component-by-component (CBC) searches.
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
   import umontreal.ssj.util.Num;

/**
 * This class implements searches to find the best rank-1 lattices with
 * respect to a given discrepancy, using component-by-component (CBC)
 * searches, random or exhaustive for each component. That is, one searches
 * for the best lattice by varying only one component at a time for each
 * dimension. Once the best component has been found for a given dimension,
 * then this value is fixed and we pass to the search for the next component.
 *
 * The discrepancy object in the `SearcherCBC` constructor *must* fix the
 * number of points @f$n@f$, the maximal dimension @f$s@f$ of the lattice and
 * possibly, the weight factors @f$\gamma_j@f$. Then the search program will
 * examine different lattices with @f$n@f$, @f$s@f$ and @f$\gamma_j@f$ fixed
 * in order to find the best amongst those examined.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class SearcherCBC extends Searcher {
   protected double[] bestVals;    // best value of discr in dim 0, 1, ..., (s-1)


   private double exhaust (int s, boolean relPrime) {
      int n = disc.getNumPoints();
      gamma = disc.getGamma();
      int pos = -1;
      double err;
      double best = 0;
      int i, j;
      bestAs[0] = 1;
      bestVals[0] = -1;

      for (j = 1; j < s; j++) {
         best = Double.MAX_VALUE;
         pos = -1;
         for (i = 1; i < n; i++) {
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
            bestAs[j] = i;
        //    print (bestAs, j + 1);
            lat = new Rank1Lattice(n, bestAs, j + 1);
            err = disc.compute (lat, gamma);
            if (err < best) {
               best = err;
               pos = i;
            }
         }
         bestAs[j] = pos;
         bestVals[j] = best;
      }

      bestVal = best;
      return best;
   }


   private  double random (int s, int k, boolean relPrime) {
      int n = disc.getNumPoints();
      if (k >= n)
         return exhaust (s, relPrime);

      gamma = disc.getGamma();
      final int nm1 = n - 1;
      int pos = -1;
      double err;
      double best = 0;
      int i, j;
      bestAs[0] = 1;
      bestVals[0] = -1;

      for (j = 1; j < s; j++) {
         best = Double.MAX_VALUE;
         i = 0;
         while (i < k) {
            if (power2F) {
               bestAs[j] = gen.nextInt (1, nm1);
               if (relPrime)
                  bestAs[j] |= 1;
            } else {
               do {
                  bestAs[j] = gen.nextInt (1, nm1);
               } while (relPrime && (Num.gcd(n, bestAs[j]) != 1));
            }
            i++;
       //     print (bestAs, j + 1);
            lat = new Rank1Lattice(n, bestAs, j + 1);
            err = disc.compute (lat, gamma);
            if (err < best) {
               best = err;
               pos = bestAs[j];
            }
         }
         bestAs[j] = pos;
         bestVals[j] = best;
      }

      bestVal = best;
      return best;
   }

   /**
    * The number of points @f$n@f$, the dimension @f$s@f$, and possibly
    * the @f$s@f$ weight factors @f$\gamma_j@f$ must be given in `disc`.
    * The @f$n@f$ points of the lattice will be generated in the search.
    * The flag `primeN` indicates whether @f$n@f$ is a prime number
    * (<tt>true</tt>) or not (<tt>false</tt>). This is used in the
    * `*Prime` methods.
    */
   public SearcherCBC (Discrepancy disc, boolean primeN) {
      super(disc, primeN);
      int s = disc.getDimension();
      bestVals = new double[s];
  }

   /**
    * Exhaustive CBC search to find the lattice with the best (the
    * smallest) discrepancy in dimension @f$s@f$. The search runs over all
    * values @f$ a_j = 1, 2, …, (n-1)@f$ for a given dimension @f$j@f$.
    * Once the best lattice has been found in dimension @f$j@f$, that
    * coefficient @f$a_j@f$ is fixed and the search runs over all values
    * of @f$a_{j+1}@f$ for the next dimension. The first component
    * @f$a_0@f$ is always set to 1. The method returns the best value of
    * the discrepancy in dimension @f$s@f$.
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
    * Random CBC search to find the lattice with the best (the smallest)
    * discrepancy in dimension @f$s@f$. @f$k@f$ random values @f$a_j@f$
    * are examined for each dimension @f$j@f$. Each random component
    * @f$a_j@f$ takes values over the integers @f$1, 2, …, (n-1)@f$. The
    * first component @f$a_0@f$ is always set to 1. The method returns the
    * best value of the discrepancy in dimension @f$s@f$.
    */
   public double random (int s, int k) {
      return random (s, k, false);
   }

   /**
    * Similar to  {@link #random(int,int) random(s, k)}, except that only
    * values of @f$a_j@f$ *relatively prime* to @f$n@f$ are considered.
    */
   public double randomPrime (int s, int k) {
      if (primeN)
         return random(s, k, false); // all a_j are rel. prime to n
      else
         return random(s, k, true);
   }

   /**
    * Returns the best value of the discrepancy found in the last search,
    * in each dimension up to @f$s@f$. The values returned are @f$V_j@f$,
    * for @f$j = 0, 1, …, (s-1)@f$.
    */
   public double[] getBestVals() {
      return bestVals;
   }

}