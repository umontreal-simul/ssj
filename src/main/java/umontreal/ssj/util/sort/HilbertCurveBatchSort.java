/*
 * Class:        HilbertCurveBatchSort
 * Description:  performs a batch sort on the arrays
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
 import umontreal.ssj.util.*;

/**
 * This sort is similar to  @ref BatchSortPow2, except that after applying
 * the batch sort, the objects are given labels that map them to the
 * @f$d@f$-dimensional unit hypercube @f$[0,1)^d@f$ as explained below, and
 * then re-ordered by following a Hilbert curve as in the
 * @ref HilbertCurveSort.
 *
 * The batch sort uses positive integers @f$n_0, …, n_{d-1}@f$, where each
 * @f$n_j@f$ must be a power of 2, say @f$n_j = 2^{e_j}@f$, and one has @f$p
 * = n_0 n_1 \cdots n_{d-1} = 2^e \ge n@f$ where @f$n@f$ is the number of
 * objects to be sorted, @f$e = e_0 + \cdots+ e_{d-1}@f$, and the @f$e_j@f$
 * are non-negative integers. During the batch sort of @f$n@f$ objects, the
 * entire space is split into @f$p = 2^e \ge n@f$ rectangular boxes, where
 * @f$e@f$ is selected as in  @ref BatchSortPow2. The batch sort then places
 * at most one object per box, by separating the objects into @f$n_0@f$
 * batches using the first coordinate, then each batch into @f$n_1@f$ batches
 * using the second coordinate, and so on.
 *
 * Each object is then mapped to a point in @f$[0,1)^d@f$, and we perform a
 * @ref HilbertCurveSort with @f$m = \max(e_0,…,e_{d-1})@f$, as follows.
 * When separating a batch with respect to coordinate @f$j@f$, the @f$n_j@f$
 * sub-batches are numbered from 0 to @f$n_j-1@f$, which are @f$e_j@f$-bit
 * numbers. We multiply them by @f$2^{m-e_j}@f$ to obtain @f$m@f$-bit numbers
 * which are interpreted as the integer coordinates of the subcubes to which
 * the corresponding points belong for the Hilbert sort. All the objects have
 * distinct vectors of integer coordinates after the batch sort. We then
 * compute the Hilbert index of each subcube as explained in the
 * @ref HilbertCurveSort class, and we sort the points as in
 * @ref HilbertCurveSort.
 *
 * The latter is implemented by using an index that maps the position (in the
 * array) of an object after the batch sort to its final position after the
 * Hilbert curve sort. For any given @f$n@f$, this mapping is always the
 * same, so it is computed once for all and saved in an index which is used
 * to made the sort. This index is recomputed only when a sort is invoked
 * with a new value of @f$n@f$.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class HilbertCurveBatchSort <T extends MultiDimComparable<? super T>> implements MultiDimSortComparable<T> {
   int dimension;  // Dimension d of the points used for the sort. 
   int m;          // Number of bit retained for each coordinate. Must be <= 31.
   HilbertCurveMap hcMap;

   BatchSortPow2 bsort;   // Batch sort used.
   long[][] indexH;       // Table whose each line gives the number of an object after the 
                          // batch sort in first coord. and its Hilbert index in second coord.
                          // Here, the second coordinates are all distinct!
                          // For a given n, this index is computed once for all and saved.
   int nSavedIndex = 0;   // Value of n for which indexH is now computed.

   /**
    * Constructs a HilbertCurveBatchSort that will use the
    * @f$\alpha_j@f$ specified in `batchesExponents` for the batch sort.
    * They must be non-negative numbers and their sum must equal 1. The
    * dimension @f$d@f$ of this array will be the dimension for the sort.
    * The Hilbert sort will use @f$m@f$ bits in each coordinate, and thus
    * @f$2^{dm}@f$ subcubes. The constructor will initialize a
    * @ref HilbertCurveMap, which can be accessed with
    * #getHilbertCurveMap.
    *  @param batchExponents proportion exponents @f$\alpha_j@f$ to
    *                        compute the batches sizes
    *  @param m              the number of bits for each coordinate
    */
   public HilbertCurveBatchSort (double[] batchExponents, int m) {
      hcMap = new HilbertCurveMap(batchExponents.length, m);
      dimension = batchExponents.length;
      this.m = m;
      bsort = new BatchSortPow2 (batchExponents);
   }

   /**
    * This constructor is similar to  #HilbertCurveBatchSort(double[],int)
    * except that the mapping object is given by the user. This given
    * @ref HilbertCurveMap instance must have the same number of
    * dimensions (which is the value of <tt>batchExponents.length</tt>)
    * and same number of bits @f$m@f$.
    *  @param batchExponents proportion exponents @f$\alpha_j@f$ to
    *                        compute the batches sizes
    *  @param m              the number of bits for each coordinate
    *  @param map            the mapping object to use
    */
   public HilbertCurveBatchSort (double[] batchExponents, int m, HilbertCurveMap map) {
      dimension = batchExponents.length;
      this.m = m;
      if (map.dimension() != dimension)
         throw new IllegalArgumentException("HilbertCurveMap has a different dimension! Expecting: " + dimension);
      this.hcMap = map;
      bsort = new BatchSortPow2 (batchExponents);
   }

   /**
    * For the given @f$n@f$, @f$e_j@f$’s, and @f$m@f$, computes and sort
    * the table that gives the Hilbert index @f$r@f$ of each object from
    * its position after the batch sort. This index table has size @f$n@f$
    * by 2, and it is constructed so its row @f$i@f$ contains @f$i@f$ in
    * the first entry and the Hilbert index of object @f$i@f$ (for the
    * given @f$m@f$) in the second entry. This table is then sorted by the
    * second coordinate and saved.
    *  @param n            number of objects
    */
   public void computeIndexH (int n) {
      if (n == nSavedIndex) return;
      nSavedIndex = n;
      indexH = new long[n][2];   // To put Hilbert indexes of points.
			// For the following, we need e_0, ..., e_{d-1}.
      int twom = (1 << m);     // 2^m.
      int mdim = dimension-1;
      int[] ej = bsort.getBitNumbers();
      int stepj[] = new int[dimension];   // The integers 2^{m-e_j}.
      int[] icoord = new int[dimension];  // The integer coordinates of current point.
      for (int j = 0; j < dimension; ++j) {
         if (ej[j] > m)
            throw new RuntimeException("ej[j] is larger than m");
	 icoord[j] = 0;
	 stepj[j] = 1 << (m - ej[j]);
      }
      indexH[0][0] = 0; 
      indexH[0][1] = 0; 
      for (int i = 1; i < n; ++i) {
         icoord[mdim] += stepj[mdim];     // Increment last coordinate. 
         if (icoord[mdim] == twom) {        // Must propagate the carry.  
            for (int j = mdim; j >= 0; --j) {
	        if (icoord[j] < twom) break; 	 // Exit for loop if no carry.			
                icoord[j] = 0;
                ++icoord[j-1]; 
            }
	 }
         indexH[i][0] = i; 
         indexH[i][1] = hcMap.coordinatesToIndex (icoord);
      }

      // use static method of HilbertCurveSort to sort
      HilbertCurveSort.sortIndexOfLong2 (indexH, 0, n); // First coordinate of indexH gives the ordering 
			                                 // that must be applied.
   }

   /**
    * Sorts the subarray `a[iMin..iMax-1]` of  @ref MultiDimComparable<T>
    * objects with this sort.
    */
   @Override
   public void sort (T[] a, int iMin, int iMax) {
        if (iMin+1 == iMax) return;
        if (iMax-iMin != nSavedIndex) 
           computeIndexH (iMax-iMin); // This changes nSavedIndex.
        bsort.sort (a, iMin, iMax);   // This changes nSaved if needed.
        // Now use the sorted indexH to sort a.  
        // First make a copy, and pick the elements from this unmodified copy. 
        T[] aclone = a.clone(); // We should not clone all the objects in a, 
	                            // but only the array of pointers. 
		// T[] aclone = new Object[iMax];

        for (int i = iMin; i < iMax; ++i) {
           a[i] = aclone[(int) indexH[i][0]];
        }
   }

   /**
    * Sorts the entire array.
    */
   public void sort (T[] a) {
      sort (a, 0, a.length);
   }
   
   @Override
   public void sort (double[][] a, int iMin, int iMax) {
        if (iMin+1 == iMax) return;
        if (iMax-iMin != nSavedIndex) 
           computeIndexH (iMax-iMin); // This changes nSavedIndex.
	    bsort.sort (a, iMin, iMax);   // This changes nSaved if needed.
        // Now use the sorted indexH to sort a.  
	    // First make a copy, and pick the elements from this unmodified copy. 
        double[][] aclone = a.clone(); 
        for (int i = iMin; i < iMax; ++i) {
            a[i] = aclone[(int) indexH[i][0]];
        }
   }
    
   @Override
   public void sort (double[][] a) {
      sort (a, 0, a.length);
   }

   /**
    * Returns the dimension of the unit hypercube.
    */
   public int dimension() {
      return dimension;
   }

    /**
     * Returns the  @ref HilbertCurveMap used for the mapping.
     */
    public HilbertCurveMap getHilbertCurveMap() {
       return hcMap;
    }

}
