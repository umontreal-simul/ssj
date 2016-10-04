/*
 * Class:        OneDimSort
 * Description:  Sorts the arrays according to one dimension
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
 * This class implements a  @ref MultiDimSortComparable that simply sorts the
 * objects according to a given sorting coordinate @f$j \ge0@f$ specified in
 * the constructor. One must have @f$j \in\{0,…,d-1\}@f$ where @f$d@f$ is
 * the dimension of the  @ref MultiDimComparable<T> objects to be sorted. The
 * sort uses the method  java.util.Arrays.sort from class  Arrays. For
 * example, to sort objects based on the first coordinate, one should use
 * `OneDimSort(0)`. It is also possible to specify in the constructor the
 * dimension @f$d@f$ (<tt>dimState</tt>) of the objects to be sorted.
 * @remark **Pierre:** Is this useful for something? Give an example. Perhaps
 * when sorting a `CachedPointSet`?
 *
 * <div class="SSJ-bigskip"></div>
 */
public class OneDimSort<T extends MultiDimComparable<? super T>> 
      implements MultiDimSortComparable<T> {
   private int j;
   private int d;

   /**
    * Constructs a OneDimSort object that will sort the
    * @ref MultiDimComparable<T> objects according to coordinate @f$j@f$. The
    * sorted objects are assumed to have dimension `d`.
    */
   public OneDimSort (int j, int d) {
      if (j < 0 || j >= d)
         throw new IllegalArgumentException(
            "sortCoordinate from 0 to " + d + " - 1.");
      this.j = j;
      this.d = d;
   }

   /**
    * Constructs a  OneDimSort that will sort the objects with
    * respect to their coordinate @f$j@f$. The dimension returned by
    * #dimension() will then be @f$d=j+1@f$.
    */
   public OneDimSort (int j) {
      this (j, j + 1);
  //    this(d - 1, d);
   }
public void sort (T[] a, int iMin, int iMax) {
      Arrays.sort (a, iMin, iMax, new MultiDimComparator<T> (j));
   }

   public void sort (T[] a) {
      sort (a, 0, a.length);
   }

   public void sort (double[][] a, int iMin, int iMax) {
      Arrays.sort (a, iMin, iMax, new DoubleArrayComparator (j));
   }

   public void sort (double[][] a) {
      sort (a, 0, a.length);
   }

/**
 * Returns the coordinate @f$j@f$ used for sorting.
 *  @return coordinate j used for sorting
 */
public int getSortCoordinate() {
      return j;
   }

   /**
    * Returns the dimension @f$d@f$ of the states to be sorted.
    *  @return dimension d of states
    */
   public int dimension() {
      return d;
   }

   /**
    * Returns a `String` containing information about this object.
    *  @return a `String`
    */
   public String toString() {
      return this.getClass().getSimpleName() + ": d = " + d +
                ", sort coordinate (starts at 0) = " + j;
   }

}
