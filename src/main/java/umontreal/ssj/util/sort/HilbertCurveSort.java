/*
 * Class:        HilbertCurveSort
 * Description:  Sorts d-dimensional points in [0,1)^d based on Hilbert curve.
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2014  Pierre L'Ecuyer and Universite de Montreal
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

  /* IMPORTANT NOTE:
	* Much of this code has been taken (with adaptations) from  
  *     the hilbert.c  code  
  * Author:	Spencer W. Thomas
  * 		EECS Dept.
  * 		University of Michigan
  * Date:	Thu Feb  7 1991
  * Copyright (c) 1991, University of Michigan
  */
package umontreal.ssj.util.sort;
  import java.util.Comparator;
  import java.util.Arrays;

/**
 * This class implements a  @ref MultiDimSort01<T extends MultiDim01> that can sort an array of
 * points in the @f$d@f$-dimensional unit hypercube @f$[0,1)^d@f$, by
 * following a Hilbert curve, and using (at most) the first @f$m@f$ bits of
 * each point. See @cite iHAM07a&thinsp;. The objects sorted by this class
 * can only be points in @f$[0,1)^d@f$, represented as arrays of `double`.
 * This sort does not apply directly to more general  @ref MultiDimComparable<T>
 * objects. For that, see the class \ref HilbertCurveBatchSort instead. 
 * However, this sort can be applied to
 * points in another space if we first define a mapping between this space
 * and the unit hypercube. For example, to sort points in the real space, it
 * suffices to map each coordinate to @f$[0,1)@f$.
 *
 * This sort (conceptually) divides the unit hypercube @f$[0,1)^d@f$ in
 * @f$2^{dm}@f$ subcubes of equal sizes, by dividing each axis in @f$2^m@f$
 * equal parts, and uses the first @f$m@f$ bits of each of the @f$d@f$
 * coordinates to place each point in one of the subcubes. It then enumerates
 * the subcubes in the same order as a Hilbert curve in @f$[0,1)^d@f$ would
 * visit them, and orders the points accordingly. Each cube has an (integer)
 * Hilbert index @f$r@f$ from 0 to @f$2^{dm}-1@f$ and the cubes (and points)
 * are ordered according to this index. Two points that fall in the same
 * subcube can be placed in an unspecified (arbitrary) order.
 *
 * It may happen that some of the subcubes contain more than one point at the
 * end. To get a sense of the probability that this happens, in the case
 * where there are @f$n@f$ points and these points are independent and
 * uniformly distributed over @f$[0,1)^d@f$, it is known that the number
 * @f$C@f$ of collisions (the number of times that a point falls in a box
 * already occupied by another point when the points are generated one after
 * the other) has approximately a Poisson distribution with mean
 * @f$\lambda_c = n^2/(2k)@f$, where @f$k = 2^{md}@f$. See
 * @cite rLEC02c&thinsp;, for example. By taking @f$m@f$ such that @f$2^{md}
 * \geq n^2@f$, the expected number of collisions is less than 1/2 and then
 * one can often neglect them. Otherwise, one should increase @f$m@f$. Note
 * however that the assumption of uniformity and independence does not always
 * hold in practice.
 *
 * For the implementation of sorts based on the Hilbert curve (or Hilbert
 * index), we identify and sort the subcubes by their Hilbert index, but it
 * is also convenient to identify them (alternatively) with @f$m@f$-bit
 * integer coordinates: The subcube with coordinates @f$(i_1,…,i_d)@f$ is
 * defined as @f$\prod_{j=0}^{d-1} [i_j 2^{-m},  (i_j+1) 2^{-m})@f$. Note
 * that each interval is open on the right. That is, if we multiply the
 * coordinates of a point in the subcube by @f$2^m@f$ and truncate them to
 * integers, we obtain the integer coordinates of the subcube. For example,
 * if @f$d=2@f$ and @f$m=4@f$, we have @f$2^8 = 256@f$ subcubes, whose
 * integer coordinates go from 0 to 15, and the point @f$(0.1, 0.51)@f$
 * belongs to the subcube with integer coordinates @f$(1, 8)@f$.
 *
 * For given @f$d@f$ and @f$m@f$, this class offers methods to compute the
 * integer coordinates of the corresponding subcube from the real-valued
 * coordinates of a point in @f$[0,1)^d@f$, as well as the Hilbert index of a
 * subcube from its integer coordinates, and vice versa. The code that
 * computes the latter correspondences is taken (with slight adaptations)
 * from the `hilbert.c` program of Spencer W. Thomas, University of Michigan,
 * 1991.
 *
 * To sort a set of @f$n@f$ points in @f$[0,1)^d@f$, we first compute the
 * integer coordinates and then the Hilbert index of the subcube for each
 * point, then sort the points by order of Hilbert index. For the latter, we
 * construct an index of type `long[n][2]` whose first coordinate is the
 * point number and the second is its Hilbert index. The method
 * #sortIndexOfLong2 is provided to sort such an index by its second
 * coordinate. Points having the same Hilbert index are ordered arbitrarily.
 * After this sort, the first coordinate at position @f$i@f$ in the index is
 * the (original) number of the point that comes in position @f$i@f$ when the
 * points are sorted by Hilbert index. This index is used to reorder the
 * original points. It can be accessed (after each sort) via the method
 * #getIndexAfterSort(). This access can be convenient for example in case we
 * sort a `double[][]` array with the Hilbert sort and want to apply the
 * corresponding permutation afterward to another array of objects. Certain
 * subclasses of HilbertCurveSort use this.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class HilbertCurveSort implements MultiDimSort01<MultiDim01> {

    int dimension;  // Dimension d of the points used for the sort. 
    long[][] indexForSort;  // This is the index computed by the last Hilbert sort,
		          // sorted by the second coordinate. The order of the first coordinates 
							// gives the permutation made by that sort. 
							// This index is recomputed each time we sort.
    HilbertCurveMap hcMap; // The map used for sorting

   /**
    * Constructs a HilbertCurveSort object that will use the first
    * @f$m@f$ bits of each of the first `d` coordinates to sort the
    * points. The constructor will initialize a  @ref HilbertCurveMap with
    * arguments `d` and @f$m@f$. This map can be accessed with
    * #getHilbertCurveMap.
    *  @param d            maximum dimension
    *  @param m            number of bits used for each coordinate
    */
   public HilbertCurveSort (int d, int m) {      
      dimension = d;
      hcMap = new HilbertCurveMap(d, m);
   }

   /**
    * Constructs a  @ref HilbertCurveSort object with a given mapping of
    * the Hilbert curve in the unit hypercube.
    *  @param map          the mapping of the Hilbert curve
    */
   public HilbertCurveSort (HilbertCurveMap map) {
      this.hcMap = map;
      this.dimension = map.dimension();
    }

    /**
     * Sorts the subarray `a[iMin..iMax-1]` with this Hilbert curve sort.
     * The type `T` must actually be  MultiDimComparable01. This is
     * verified in the method.
     */
    @Override
    public void sort (MultiDim01[] a, int iMin, int iMax) {
			 // Copy the (0,1)^d transformations of the states in array b.
       double b[][] = new double[iMax][dimension];
       for (int i = iMin; i < iMax; ++i)
          b[i] = a[i].getPoint();
			 // Sort this array b by Hilbert sort.  The index 
			 // indexForSort will contain the permutation made by that sort.
       sort (b, iMin, iMax);
       // Now use indexForSort to sort a.
       // We do not want to clone all the objects in a, 
       // but only the array of pointers. 
       MultiDim01[] aclone = a.clone();    // new Object[iMax];
       for (int i = iMin; i < iMax; ++i)
           a[i] = aclone[(int) indexForSort[i][0]];
    }

    /**
     * Sorts the entire array: same as `sort (a, 0, a.length)`.
     */
    @Override
    public void sort (MultiDim01[] a) {
        sort (a, 0, a.length);
    }

    /**
     * Sort the array with Hilbert sort.
     */
    public void sort (double[][] a, int iMin, int iMax) {
        if (iMin+1 == iMax) return;
        indexForSort = new long[iMax][2];    // Index used for sort.
        int[] icoord = new int[dimension];   // To store integer coordinates.
        for (int i=0; i < a.length; ++i) {   
				    // Transform to integer coordinates.
            hcMap.pointToCoordinates (a[i], icoord);
            indexForSort[i][0] = i;
            indexForSort[i][1] = hcMap.coordinatesToIndex (icoord);  // Hilbert index of this point. 
        }
        // Sort the index based on the positions on the Hilbert curve
        sortIndexOfLong2 (indexForSort, iMin, iMax);
        // Now use the sorted index to sort a.
        double[][] aclone = a.clone();       // Save copy of a before the sort.
        for (int i= iMin; i< iMax; ++i) {
            a[i] = aclone[(int) indexForSort[i][0]];
        }
    }
public void sort (double[][] a) {
       sort (a, 0, a.length);
    }

/**
 * Returns the index computed by the last sort, which is sorted by the second
 * coordinate. It contains (in the first coordinates of its entries) the
 * permutation made by that sort.
 */
public long[][] getIndexAfterSort () {
       return indexForSort;
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

  // Compares two arrays of long according to their second coordinate.
  // This is used to sort an index of type long[][2] by the second coordinate.
  // The permutation can be recovered in the first coordinate.

/**
 * The comparator class used by  #sortIndexOfLong2.
 */
public static class LongIndexComparator2 implements Comparator<long[]> {
     @Override
     public int compare (long[] p1, long[] p2) {
        if (p1[1] > p2[1]) return 1;
        else if (p1[1] < p2[1])  return -1;
        else return 0;
     }
  }

  /**
   * Sorts the `index` table by its second coordinate.
   */
  public static void sortIndexOfLong2 (long[][] index, int iMin, int iMax) {
      // if (iMin==(iMax-1)) return;
      Arrays.sort (index, iMin, iMax, new LongIndexComparator2());        
  }

}
