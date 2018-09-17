/*
 * Class:        BatchSort
 * Description:  performs a batch sort on the arrays
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
package umontreal.ssj.util.sort;
 import java.util.Comparator;
 import java.util.Arrays;

/**
 * This class implements a \ref MultiDimSortComparable that performs a batch
 * sort on multivariate arrays. It separates the objects in @f$n_0@f$ batches
 * of (approximately) equal size such that each object in a batch is smaller
 * or equal, in the first coordinate, to the objects in the next batches.
 * Then, each batch is separated in @f$n_1@f$ batches in the same way but
 * using the second coordinate. And so on.
 *
 * One way of constructing a @ref BatchSort object is to specify the batch
 * numbers (integers) @f$n_0, n_1, …, n_{d-1}@f$ in an array of size @f$d@f$
 * named `batchNumbers` passed to the constructor. The size @f$d@f$ must not
 * exceed the dimension of the objects to be sorted, and the product @f$p =
 * n_0 n_1 \cdots n_{d-1}@f$ must normally be equal to the number @f$n@f$ of
 * objects to be sorted, which is `iMax - iMin` when doing a partial sort,
 * and the total number of objects otherwise. This means that one should
 * normally use always the same @f$n@f$ after the BatchSort object has
 * been constructed in this way, with a fixed @f$p@f$. The sorting method
 * always uses batch sizes at each level as if the number of objects was
 * exactly @f$p@f$. The batch size at level @f$j@f$ is @f$n_{j+1}
 * \cdots n_{d-1}@f$. If there are less than @f$p@f$ objects, some batches
 * (the last ones) can have fewer objects or be empty. What is done is
 * conceptually equivalent to adding dummy objects at the end of the array so
 * their total number is exactly @f$p@f$. For example, if @f$n_0 = 3@f$,
 * @f$p=15@f$ and @f$n=13@f$, then the batch size will be @f$15/3 = 5@f$ but
 * the last batch will have only 3 objects. The value of @f$n_{d-1}@f$ is
 * actually never used by the batch sort; at the last level, each batch is
 * sorted by the last coordinate regardless of its size.
 *
 * A second way of constructing a BatchSort object is by specifying an
 * array of non-negative doubles @f$\alpha_0,\alpha_1,…,\alpha_{d-1}@f$,
 * called *proportion exponents*, such that
 * @f$\alpha_0+\cdots+\alpha_{d-1} = 1@f$. For each number @f$n@f$ of
 * objects to be sorted, the batch numbers will be recomputed as
 * (approximately) @f$n_j = n^{\alpha_j}@f$. With this construction, one can
 * easily handle very different (arbitrary) values of @f$n@f$ with the same
 * BatchSort object. The vector of batch numbers is recomputed each time
 * we sort with a new value of @f$n@f$, as follows: @f$n_0 =
 * \lceil n^{\alpha_0} \rceil@f$, @f$n_1 = \lceil(n/n_0)^{\alpha_1}
 * \rceil@f$, …, @f$n_{d-1} = \lceil n / (n_0 \cdots n_{d-2}) \rceil@f$.
 * These values are saved in case we use the same @f$n@f$ the next time. They
 * give @f$p = n_0 \cdots n_{d-1} \ge n@f$. The batch size at level @f$j@f$
 * is again @f$n_{j+1} \cdots n_{d-1}@f$. When @f$n < p@f$, some batches (the
 * last ones) have fewer objects than the others.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class BatchSort<T extends MultiDimComparable<? super T>> implements MultiDimSortComparable<T> {
   int dimension;        // Number of coordinates.
   boolean useExponents; // True if we use proportion exponents.
   int[] batchNumbers;   // Number of batches n_j for each dimension.
   double[] batchExponents;  // The alpha_j.
   int batchProduct = 1; // Product p of numbers in batchNumbers.
   int nSaved = 0;       // Number n of objects last time sort was called.

   /**
    * Constructs a BatchSort that will always use the (fixed) batch
    * numbers given in `batchNumbers`. These batch numbers can be changed
    * only by creating new object, of by using proportion exponents via
    * `setBatchExponents`. The number of objects to sort should not exceed
    * the product of the numbers in `batchNumbers`.
    *  @param batchNumbers number of batches in each dimension
    */
   public BatchSort (int[] batchNumbers) {
      if (batchNumbers == null)
         throw new NullPointerException("batchNumbers is null");
      useExponents = false;
      this.batchNumbers = batchNumbers;
      dimension = batchNumbers.length;
      batchProduct = 1;
      for(int i = 0; i < batchNumbers.length; ++i) 
	batchProduct *= batchNumbers[i];
   }

   /**
    * Constructs a BatchSort that will use the proportion exponents
    * in `batchExponents`, which must contain non-negative numbers that
    * sum to 1.
    *  @param batchExponents proportion exponents to compute the batches
    *                        sizes
    */
   public BatchSort (double[] batchExponents) {
      if (batchExponents == null)
         throw new NullPointerException("batchExponents is null");
      this.batchExponents = batchExponents;
			useExponents = true;
      dimension = batchExponents.length;
      batchNumbers = new int[dimension];
      double sum = 0.0;
      for (int j=0; j < dimension; ++j)
         sum += batchExponents[j];
      if (Math.abs(sum-1) > 1.0e-10)
         throw new IllegalArgumentException("Sum of batchExponents not equal to 1");
   }

   /**
    * Resets the corresponding vector of batch numbers for the given
    * number @f$n@f$ of objects. This method can be invoked only when
    * proportion exponents are used and have been defined previously.
    *  @param n            number of objects
    */
   public void setBatchNumbers (int n) {
      if (batchExponents == null)
         throw new NullPointerException("batchExponents is null");
      if (useExponents == false)
         throw new IllegalArgumentException("method allowed only when using proportion exponents");
      nSaved = n;
			batchProduct = 1;       // Product of batch numbers.
      for (int dim = 0; dim < dimension; ++dim) {
         batchNumbers[dim] = (int) Math.ceil (Math.pow (n, batchExponents[dim])); // Round up.
         batchProduct *= batchNumbers[dim];
      }
   }

   /**
    * Sets the vector of proportion exponents
    * @f$\alpha_0,\alpha_1,…,\alpha_{d-1}@f$ to the values in
    * `batchExponents`. These values must be non-negative and sum to 1.
    * From now on, this BatchSort object will use proportion
    * exponents.
    *  @param batchExponents exponents in each dimension
    */
   public void setBatchExponents (double[] batchExponents) {
      if (batchExponents == null)
         throw new NullPointerException("batchExponents is null");
			if (batchExponents.length != dimension)
         throw new IllegalArgumentException("batchExponents has wrong dimension");
		  nSaved = 0;
      useExponents = true;
			this.batchExponents = batchExponents;
			double sumExponents = 0.0;  // Sum of exponents.
      for (int dim = 0; dim < dimension; ++dim) {
         sumExponents += batchExponents[dim];				
      }
			if (Math.abs(sumExponents-1.0) > 1.0e-10)
         throw new IllegalArgumentException("Sum of batchExponents not equal to 1");
   }

   /**
    * Returns the current vector of batch numbers @f$n_j@f$.
    */
   public int[] getBatchNumbers () {
			return batchNumbers;
   }

   /**
    * Returns the product @f$p@f$ of current batch numbers.
    */
   public int getBatchProduct () {
			return batchProduct;
   }

   /**
    * Returns the current vector of batch exponents @f$\alpha_j@f$.
    */
   public double[] getBatchExponents () {
			return batchExponents;
   }
public int dimension() {
      return dimension;
   }

/**
 * Sorts the subarray `a[iMin..iMax-1]` using this batch sort.
 */
public void sort (T[] a, int iMin, int iMax) {
      if (iMin+1 == iMax) return;
      if (useExponents && (nSaved != iMax-iMin))
				 setBatchNumbers (iMax-iMin);
      int i1, i2; 
			int bsize = iMax-iMin;  // Current batch size. At level 0: single batch.
      for (int j = 0; (j < dimension) && (bsize > 1); ++j) {
         MultiDimComparator<T> compar = new MultiDimComparator<T>(j);
         // Sort all the batches (entire array) according to coordinate j.
         if (batchNumbers[j]==1) continue;  // Can skip the sort for this dim j.
         // if (bsize <= 1) break;
         i1 = iMin;
         while (i1 < iMax) {
				    i2 = i1 + bsize;  
						if (i2 > iMax) i2 = iMax;
            Arrays.sort (a, i1, i2, compar);
            i1 += bsize;
         }
         bsize = (int) Math.ceil (bsize / batchNumbers[j]);  // New batch size.
      }
   }

   /**
    * Sorts the entire array.
    */
   public void sort (T[] a) {
      sort (a, 0, a.length);
   }

   /**
    * Sorts the subarray `a[iMin..iMax-1]` using this batch sort.
    */
   public void sort (double[][] a, int iMin, int iMax) {
      if (iMin+1 == iMax) return;
      if (useExponents && (nSaved != iMax-iMin))
				 setBatchNumbers (iMax-iMin);
      int i1, i2; 
			int bsize = iMax-iMin;  // Current batch size. At level 0: single batch.
      for (int j = 0; (j < dimension) && (bsize > 1); ++j) {
         DoubleArrayComparator compar = new DoubleArrayComparator (j);
         // Sort all the batches (entire array) according to coordinate j.
         if (batchNumbers[j]==1) continue;  // Can skip the sort for this dim j.
         i1 = iMin;
         while (i1 < iMax) {
				    i2 = i1 + bsize;  
						if (i2 > iMax) i2 = iMax;
            Arrays.sort (a, i1, i2, compar);
            i1 += bsize;
         }
         bsize = (int) Math.ceil (bsize / batchNumbers[j]);  // New batch size.
      }
   }

   /**
    * Sorts the entire array.
    */
   public void sort (double[][] a) {
      sort (a, 0, a.length);
   }

}
