/*
 * Class:        MultiDimSortComparable
 * Description:  Represents a sort on multidimensional arrays of MultiDimComparable.
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
 * This interface extends  @ref MultiDimSort<T> to implement multivariate
 * sorting algorithms that sort objects that are pairwise comparable. These
 * objects have @f$d@f$ sorting fields (or coordinates) and must implement
 * the interface  @ref MultiDimComparable, which is used to sort them on any
 * given coordinate, numbered from 0 to @f$d-1@f$.
 *
 * The ordering is induced by the method
 * umontreal.ssj.util.MultiDimComparable.compareTo, and the number of the
 * largest coordinate used by the sort must not exceed @f$d-1@f$, where the
 * dimension @f$d@f$ refers to the value returned by
 * umontreal.ssj.util.MultiDimComparable.dimension. One can sort only a
 * subset of the objects, or all of them.
 *
 * <div class="SSJ-bigskip"></div>
 */
public interface MultiDimSortComparable <T extends MultiDimComparable<? super T>> extends 
   MultiDimSort<T> {

}
