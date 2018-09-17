/*
 * Class:        SplitSort
 * Description:  Performs a split sort on arrays.
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
 * Implements a  @ref MultiDimSortComparable that performs a *split sort* on
 * a  @ref MultiDimComparable<T> array based on its first @f$d@f$ dimensions.
 * The sort first separates the array in two, such that the objects in the
 * first part are smaller or equal to those in the second part for their
 * coordinate 0. If the number of objects is even, the two parts will contain
 * the same number of objects, otherwise there will be one more in the second
 * part. Then each part is recursively separated in two based on coordinate
 * 1, and so on, until each part contains a single object. If not done yet
 * after splitting on coordinate @f$d-1@f$, the procedure cycles back to
 * coordinate 0, and continues. The resulting order is the result of the
 * sort.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class SplitSort<T extends MultiDimComparable<? super T>> 
      implements MultiDimSortComparable<T> {
   private int dimension;

   /**
    * Constructs a SplitSort that will use the first `d` dimensions
    * to sort.
    *  @param d            number of dimensions to use for the sort
    */
   public SplitSort (int d) {
      this.dimension = d;
   }
   public void sort (T[] a, int iMin, int iMax) {
      if (dimension == 1) {
         Arrays.sort (a, iMin, iMax, new MultiDimComparator<T>(0));
      } else {
         splitSort (a, iMin, iMax, 0);
      }
   }

   private void splitSort (T[] a, int iMin, int iMax, int splitCoord){
      if (iMin==(iMax-1)) return;
      Arrays.sort (a, iMin, iMax, new MultiDimComparator<T>(splitCoord));
      int iMid = (iMin+iMax)/2;
      splitSort (a, iMin, iMid, (splitCoord+1)%dimension);
      splitSort (a, iMid, iMax, (splitCoord+1)%dimension);
   }

   public void sort (T[] a) {
      sort (a, 0, a.length);
   }

   public void sort (double[][] a, int iMin, int iMax) {
      if (dimension == 1) {
         Arrays.sort (a, iMin, iMax, new DoubleArrayComparator(0));
      }else{
         splitSort (a, iMin, iMax, 0);
      }
   }

   private void splitSort (double[][] a, int iMin, int iMax, int splitCoord){
      if (iMin==(iMax-1)) return;
      Arrays.sort (a,iMin,iMax,new DoubleArrayComparator(splitCoord));
      int iMid = (iMin+iMax)/2;
      splitSort (a,iMin,iMid,(splitCoord+1)%dimension);
      splitSort (a,iMid,iMax,(splitCoord+1)%dimension);
   }

   public void sort (double[][] a) {
      sort (a, 0, a.length);
   }

   public int dimension() {
      return dimension;
   }
}
