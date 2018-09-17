/*
 * Class:        MultiDimSort
 * Description:  Represents a sort on multidimensional arrays.
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
 * This interface is meant to be implemented by certain multivariate sorting
 * algorithms that sort objects based on different fields (or dimensions).
 * For example, in a list of animals, the animals can be sorted by both their
 * age, their weight, and their market value, in three dimensions. For
 * instance, one may want to regroup them based on their age, then split each
 * group based on the weight, and finally sort the subgroups based on market
 * value.
 *
 * The objects to be sorted have @f$d@f$ sorting fields (or coordinates),
 * which means that they can be seen as @f$d@f$-dimensional and can be sorted
 * according to any of the @f$d@f$ coordinates, which are numbered from 0 to
 * @f$d-1@f$. These objects can be simply @f$d@f$-dimensional vectors of
 * <tt>double</tt>’s, which are sorted in the natural way according to any of
 * the @f$d@f$ coordinates. They can also be more general types of objects
 * that implement the  @ref MultiDimComparable<T> interface, which is used to
 * sort them on any given coordinate. The ordering is then the one induced by
 * the method  \ref MultiDimComparable<T>#compareTo, and the
 * number of the largest coordinate used by the sort must not exceed
 * @f$d-1@f$, where the dimension @f$d@f$ refers to the value returned by
 * umontreal.ssj.util.MultiDimComparable.dimension. One can sort only a
 * subset of the objects, or all of them.
 *
 * Certain types of sorts can use several coordinates at the same time to
 * sort the objects; see the classes  @ref BatchSort and
 * @ref HilbertCurveSort, for example. These types of sorts have an
 * application in particular in the array-RQMC method to simulate Markov
 * chains @cite vLEC08a&thinsp;.
 *
 * <div class="SSJ-bigskip"></div>
 */
public interface MultiDimSort<T> {

   /**
    * Sorts the subarray of `a` made of the elements with indices from
    * `iMin` to `iMax-1`.
    *  @param a            array to sort
    *  @param iMin         index of first element to sort
    *  @param iMax         index of last element to sort is
    *                      @f$\mathtt{iMax}-1@f$
    */
   public void sort (T[] a, int iMin, int iMax);

   /**
    * Sorts the entire array `a`.
    *  @param a            array to sort
    */
   public void sort (T[] a);

   /**
    * Sorts the subarray of `a` made of the elements with indices from
    * `iMin` to `iMax-1`, using the natural order for real numbers for
    * each coordinate.
    *  @param a            array to sort
    *  @param iMin         index of first element to sort
    *  @param iMax         index of last element to sort is
    *                      @f$\mathtt{iMax}-1@f$
    */
   public void sort (double[][] a, int iMin, int iMax);

   /**
    * Same as above, but sorts the entire array.
    *  @param a            array to sort
    */
   public void sort (double[][] a);

   /**
    * Returns the number of dimensions used in the sort.
    *  @return number of dimensions used for the sort
    */
   public int dimension();

}
