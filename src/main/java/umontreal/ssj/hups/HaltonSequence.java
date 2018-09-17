/*
 * Class:        HaltonSequence
 * Description:  
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

/**
 * This class implements the sequence of Halton @cite rHAL60a&thinsp;, which
 * is essentially a modification of Hammersley nets for producing an infinite
 * sequence of points having low discrepancy. The @f$i@f$th point in @f$s@f$
 * dimensions is
 * @anchor REF_hups_HaltonSequence_eq_Halton_point2
 * @f[
 *   \mathbf{u}_i = (\psi_{b_1}(i),\psi_{b_2}(i),…, \psi_{b_s}(i)), \tag{Halton-point2}
 * @f]
 * for @f$i=0,1,2,…@f$, where @f$\psi_b@f$ is the radical inverse function
 * in base @f$b@f$, defined in class  @ref RadicalInverse, and where @f$2 =
 * b_1 < \cdots< b_s@f$ are the @f$s@f$ smallest prime numbers in increasing
 * order.
 *
 * A fast method is implemented to generate randomized Halton sequences
 * @cite rSTR95a, @cite vWAN99a&thinsp;, starting from an arbitrary point
 * @f$x_0@f$.
 *
 * The points can be "scrambled" by applying a permutation to the digits of
 * @f$i@f$ before computing each coordinate via (
 * {@link REF_hups_overview_eq_Halton_point Halton-point} ), in
 * the same way as for the class  @ref HammersleyPointSet, for all
 * coordinates @f$j\ge0@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class HaltonSequence extends PointSet { 
   private int[] base;           // Vector of prime bases.
   private int[][] permutation;  // Digits permutation, for each dimension.
   private boolean permuted;     // Permute digits?
   private RadicalInverse[] radinv; // Vector of RadicalInverse's.
   private int[] start;          // starting indices
   private final static int positiveBitMask = ~Integer.reverse(1);

   /**
    * Constructs a new Halton sequence in `dim` dimensions.
    *  @param dim          dimension
    */
   public HaltonSequence (int dim) {
      if (dim < 1)
         throw new IllegalArgumentException
            ("Halton sequence must have positive dimension dim");
      this.dim  = dim;
      numPoints = Integer.MAX_VALUE;
      base = RadicalInverse.getPrimes (dim);
      start = new int[dim];
      java.util.Arrays.fill(start, 0);
   }

   /**
    * Initializes the Halton sequence starting at point `x0`. For each
    * coordinate @f$j@f$, the sequence starts at index @f$i_j@f$ such that
    * <tt>x0[</tt>@f$j@f$<tt>]</tt> is the radical inverse of @f$i_j@f$.
    * The dimension of `x0` must be at least as large as the dimension of
    * this object.
    *  @param x0           starting point of the Halton sequence
    */
   public void setStart (double[] x0) {
      for (int i = 0; i < dim; i++)
         start[i] = RadicalInverse.radicalInverseInteger(base[i], x0[i]);
   }

   /**
    * Initializes the Halton sequence starting at point `x0`. The
    * dimension of `x0` must be at least as large as the dimension of this
    * object.
    *  @param x0           starting point of the Halton sequence
    */
   public void init (double[] x0) {
      radinv = new RadicalInverse[dim];
      for (int i = 0; i < dim; i++)
         radinv[i] = new RadicalInverse (base[i], x0[i]);
   }

   /**
    * Permutes the digits using permutations from @cite vFAU09a&thinsp;
    * for all coordinates. After the method is called, the coordinates
    * @f$u_{i,j}@f$ are generated via
    * @f[
    *   u_{i,j} = \sum_{r=0}^{k-1} \pi_j[a_r] b_j^{-r-1},
    * @f]
    * for @f$j=0,…,s-1@f$, where @f$\pi_j@f$ is the Faure-Lemieux (2008)
    * permutation of @f$\{0,…,b_j-1\}@f$.
    */
   public void addFaureLemieuxPermutations() {
      permutation = new int[dim][];
      for (int i = 0; i < dim; i++) {
         permutation[i] = new int[base[i]];
         RadicalInverse.getFaureLemieuxPermutation (i, permutation[i]);
      }
      permuted = true;
   }

   /**
    * Permutes the digits using Faure permutations for all coordinates.
    * After the method is called, the coordinates @f$u_{i,j}@f$ are
    * generated via
    * @f[
    *   u_{i,j} = \sum_{r=0}^{k-1} \pi_j[a_r] b_j^{-r-1},
    * @f]
    * for @f$j=0,…,s-1@f$, where @f$\pi_j@f$ is the Faure permutation of
    * @f$\{0,…,b_j-1\}@f$.
    */
   public void addFaurePermutations() {
      permutation = new int[dim][];
      for (int i = 0; i < dim; i++) {
         permutation[i] = new int[base[i]];
         RadicalInverse.getFaurePermutation (base[i], permutation[i]);
      }
      permuted = true;
   }

   /**
    * Erases the permutations: from now on, the digits will not be
    * permuted.
    */
   public void ErasePermutations() {
      permuted = false;
      permutation = null;
   }

    
   public int getNumPoints () {
      return Integer.MAX_VALUE;
   }

   public double getCoordinate (int i, int j) {
      if (radinv != null) {
         if (!permuted) {
            return radinv[j].nextRadicalInverse ();
         } else {
            throw new UnsupportedOperationException (
            "Fast radical inverse is not implemented in case of permutation");
         }
      } else {
         int k = start[j] + i;
         // if overflow, restart at first nonzero point
         // (Struckmeier restarts at zero)
         if (k < 0)
            k = (k & positiveBitMask) + 1;
         if (permuted)
            return RadicalInverse.permutedRadicalInverse 
            (base[j], permutation[j], k);
         else 
            return RadicalInverse.radicalInverse (base[j], k);
      }
   }
}