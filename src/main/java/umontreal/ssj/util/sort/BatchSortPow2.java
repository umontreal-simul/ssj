/*
 * Class:        BatchSortPow2
 * Description:  performs a batch sort on an array, batch numbers must be powers of 2.
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2015  Pierre L'Ecuyer and Universite de Montreal
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
package umontreal.ssj.util.sort;
 import java.util.Comparator;
 import java.util.Arrays;

/**
 * This is a subclass of  @ref BatchSort for which the batch numbers
 * @f$n_j@f$ are always powers of 2. One has @f$n_j = 2^{e_j}@f$ for all
 * @f$j@f$ and @f$p = 2^e@f$, where the @f$e_j@f$ and @f$e@f$ are integers.
 *
 * The batch numbers can only be specified via the proportion exponents
 * @f$\alpha_0,\alpha_1,…,\alpha_{d-1}@f$. The @f$e_j@f$ are computed as a
 * function of the number @f$n@f$ of objects to be sorted. This number
 * @f$n@f$ is first rounded up the nearest power of 2, say @f$p = 2^e
 * \ge n@f$. (Often, it should already be a power of 2.) Then the @f$e_j@f$
 * are defined as @f$e_0 = \lceil e \alpha_0 \rceil@f$ and @f$e_j =
 * \lceil(e - e_0 - \cdots- e_{j-1}) (\alpha_0 + \cdots+ \alpha_j)
 * \rceil@f$ for @f$1\le j < d@f$. This gives @f$e_0 + \cdots+ e_{d-1} =
 * e@f$.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class BatchSortPow2 extends BatchSort {
   int[] ej;   // The current exponents e_j.

   /**
    * Constructs a  @ref BatchSortPow2 that will use proportion exponents
    * in `batchExponents`, which must contain non-negative numbers that
    * sum to 1.
    *  @param batchExponents proportion exponents to compute the batches
    *                        sizes
    */
   public BatchSortPow2 (double[] batchExponents) {
      super (batchExponents);
      ej = new int[this.dimension];
   }

   /**
    * For a number @f$n@f$ of objects and a predefined vector of
    * proportion exponents @f$\alpha_0,\alpha_1,…,\alpha_{d-1}@f$,
    * computes and sets the corresponding vector of batch numbers, for
    * which each batch number is a power of 2 and their product is @f$p =
    * 2^e \ge n@f$.
    *  @param n            number of objects
    */
   public void setBatchNumbers (int n) {
      if (batchExponents == null)
         throw new NullPointerException("batchExponents is null");
      // Round up to the nearest power of 2.  Will be 2^e.
      int e = 1;
			while (n > (1 << e)) e++; 
			nSaved = n;
			batchProduct = 1 << e; // Product of batch numbers = p = 2^e.
      for (int j = 0; j < dimension; ++j) {
         ej[j] = (int) Math.ceil (e * batchExponents[j]); // log_2 of num of batches.
         e -= ej[j];   // The power of 2 that remains for each batch.
				 batchNumbers[j] = 1 << ej[j];
      }
   }

   /**
    * Returns the current vector of integers @f$e_j@f$.
    */
   public int[] getBitNumbers () {
			return ej;
   }

}
