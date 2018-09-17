/*
 * Class:        HammersleyPointSet
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
 * This class implements *Hammersley point sets*, which are defined as
 * follows. Let @f$2 = b_1 < b_2 < \cdots@f$ denote the sequence of all
 * prime numbers by increasing order. The Hammersley point set with @f$n@f$
 * points in @f$s@f$ dimensions contains the points
 * @anchor REF_hups_HammersleyPointSet_eq_Hammersley_point2
 * @f[
 *   \mathbf{u}_i = (i/n,\psi_{b_1}(i),\psi_{b_2}(i),…, \psi_{b_{s-1}}(i)), \tag{Hammersley-point2}
 * @f]
 * for @f$i=0,…,n-1@f$, where @f$\psi_b@f$ is the radical inverse function
 * in base @f$b@f$, defined in  @ref RadicalInverse. This class is not a
 * subclass of  @ref DigitalNet, because the basis is not the same for all
 * coordinates. We do obtain a net in a generalized sense if @f$n = b_1^{k_1}
 * = b_2^{k_2} = \cdots= b_{s-1}^{k_{s-1}}@f$ for some integers
 * @f$k_1,…,k_{s-1}@f$.
 *
 * The points of a Hammersley point set can be "scrambled" by applying a
 * permutation to the digits of @f$i@f$ before computing each coordinate via
 * ( {@link REF_hups_overview_eq_Hammersley_point
 * Hammersley-point} ). If
 * @f[
 *   i = a_0 + a_1 b_j + …+ a_{k_j-1} b_j^{k_j-1},
 * @f]
 * and @f$\pi_j@f$ is a permutation of the digits @f$\{0,…,b_j-1\}@f$, then
 * @f[
 *   \psi_{b_j}(i) = \sum_{r=0}^{k_j-1} a_r b_j^{-r-1}
 * @f]
 * is replaced in (
 * {@link REF_hups_overview_eq_Hammersley_point
 * Hammersley-point} ) by
 * @f[
 *   u_{i,j}= \sum_{r=0}^{k_j-1} \pi_j[a_r] b_j^{-r-1}.
 * @f]
 * The permutations @f$\pi_j@f$ can be deterministic or random. One
 * (deterministic) possibility implemented here is to use the Faure
 * permutation of @f$\{0,…,b_j\}@f$ for @f$\pi_j@f$, for each coordinate
 * @f$j > 0@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class HammersleyPointSet extends PointSet {
   private int[] base;           // Vector of prime bases.
   private int[][] permutation;  // Digits permutation, for each dimension.
   private boolean permuted;     // Permute digits?

   /**
    * Constructs a new Hammersley point set with `n` points in `dim`
    * dimensions.
    *  @param n            number of points
    *  @param dim          dimension of the point set
    */
   public HammersleyPointSet (int n, int dim) {
      if (n < 0 || dim < 1)
         throw new IllegalArgumentException
            ("Hammersley point sets must have positive dimension and n >= 0");
      numPoints = n;
      this.dim  = dim;
      if (dim > 1)
         base = RadicalInverse.getPrimes (dim - 1);
   }

   /**
    * Permutes the digits using Faure permutations for all coordinates.
    * After the method is called, the coordinates @f$u_{i,j}@f$ are
    * generated via
    * @f[
    *   u_{i,j} = \sum_{r=0}^{k-1} \pi_j[a_r] b_j^{-r-1},
    * @f]
    * for @f$j=1,…,s-1@f$ and @f$u_{i,0}=i/n@f$, where @f$\pi_j@f$ is the
    * Faure permutation of @f$\{0,…,b_j-1\}@f$.
    */
   public void addFaurePermutations() {
      if (dim > 1) {
         permutation = new int[dim][];
         for (int i = 0; i < dim - 1; i++) {
            permutation[i] = new int[base[i]];
            RadicalInverse.getFaurePermutation (base[i], permutation[i]);
         }
      }
      permuted = true;
   }

   /**
    * Erases the Faure permutations: from now on, the digits will not be
    * Faure permuted.
    */
   public void ErasePermutations() {
      permuted = false;
      permutation = null;
   }


   public double getCoordinate (int i, int j) {
      if (j == 0)
         return (double) i / (double) numPoints;
      if (permuted)
         return RadicalInverse.permutedRadicalInverse 
                   (base[j - 1], permutation[j - 1], i);
      else 
         return RadicalInverse.radicalInverse (base[j - 1], i);
   }
}