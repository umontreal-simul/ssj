/*
 * Class:        KorobovLatticeSequence
 * Description:  
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since

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
package umontreal.ssj.hups;

/**
 * @remark **Pierre:** This class is not yet fully implemented
 *
 * This class implements Korobov lattice sequences, defined as follows. One
 * selects a *basis* @f$b@f$ and a (large) multiplier @f$a@f$. For each
 * integer @f$k\ge0@f$, we may consider the @f$n@f$-point Korobov lattice
 * with modulus @f$n = b^k@f$ and multiplier @f$\tilde{a} = a \bmod n@f$. Its
 * points have the form
 * @anchor REF_hups_KorobovLatticeSequence_eq_Korobov_seq1
 * @f[
 *   \mathbf{u}_i = (a^i (1, a, a^2, …) \bmod n) / n = (\tilde{a}^i (1, \tilde{a}, \tilde{a}^2, …) \bmod n) / n \tag{Korobov-seq1}
 * @f]
 * for @f$i=0,…,n-1@f$. For @f$k = 0,1,…@f$, we have an increasing sequence
 * of lattices contained in one another.
 *
 * These embedded lattices contain an infinite sequence of points that can be
 * enumerated as follows @cite vHIC01a&thinsp;:
 * @anchor REF_hups_KorobovLatticeSequence_eq_Korobov_seq2
 * @f[
 *   \mathbf{u}_i = \psi_b(i) \left(1, a, a^2, …\right) \bmod1. \tag{Korobov-seq2}
 * @f]
 * where @f$\psi_b(i)@f$ is the radical inverse function in base @f$b@f$,
 * defined in  @ref RadicalInverse. The first @f$n=b^k@f$ points in this
 * sequence are exactly the same as the @f$n@f$ points in (
 * {@link REF_hups_KorobovLatticeSequence_eq_Korobov_seq1
 * Korobov-seq1} ), for each @f$k\ge0@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class KorobovLatticeSequence extends KorobovLattice { 
   int base;         // Base for radical inversion
   int inverse;      // global variables for radical inverssion,
   int n;            // since bloody JAVA cannot pass references

   // Method modPower is inherited from Rank1Lattice.

   /**
    * Constructs a new lattice sequence with base `b` and `generator` @f$
    * = a@f$.
    *  @param b            number of points (modulus) is a power of b
    *  @param a            multiplier @f$a@f$ of this lattice sequence
    */
   public KorobovLatticeSequence (int b, int a) {
// Pas termine: ne fonctionne pas
      super (2, 3, 1);
      if (a < 1)
         throw new IllegalArgumentException
             ("KorobovLatticeSequence:   Multiplier a must be >= 1");
//      dim       = Integer.MAX_VALUE;
//      numPoints = Integer.MAX_VALUE;
      base = b;
throw new UnsupportedOperationException ("NOT FINISHED");
   }
 

   // A very inefficient way of generating the points!
   public double getCoordinate (int i, int j) {
      int n;
      int inverse;
      if (i == 0)
         return 0.0;
      else if (j == 0)
         return radicalInverse (base, i);
      else {
         // integerRadicalInverse (i);
         n = 1;
         for (inverse = 0; i > 0; i /= base) {
            inverse = inverse * base + (i % base);
            n *= base;
         }
         return (double) ((inverse * modPower (genA, j, n)) % n) / (double) n;
      }
   }

   // ... has been unrolled in getCoordinate.
   private void integerRadicalInverse (int i) {
      // Attention: returns results in variables n and inverse.
      n = 1;
      for (inverse = 0; i > 0; i /= base) {
         inverse = inverse * base + (i % base);
         n *= base;
      }
   }
 
}