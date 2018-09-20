/*
 * Class:        MultiDimComparable
 * Description:  Represents an object which can be compared in many dimensions.
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

/**
 * This interface is an extension (or variant) of the  Comparable interface
 * in Java. In the  Comparable interface, the objects can be compared in only
 * one way, so there is a single total ordering between the objects. Here,
 * there is a choice of @f$d@f$ different ways (or criteria, or coordinates)
 * to compare the objects. These coordinates are indexed by @f$j =
 * 0,…,d-1@f$. The coordinate @f$j@f$ that is to used by the comparison must
 * be specified as the second argument in the method  {@link #compareTo()
 * compareTo(other, j)}. This interface must be implemented in particular by
 * objects that are to be sorted via a  @ref MultiDimSortComparable.
 *
 * The objects that are compared may really have a multidimensional behavior,
 * or they may not. For example, an object that contains an array can
 * implement this interface by defining  {@link #compareTo() compareTo(other,
 * j)} so that it compares the elements of index `j` in the arrays from the
 * two objects. As another example, if an object represents a person and
 * contains some information on that person, the method  #compareTo may
 * compare two persons based on the age, or height, er weight, etc., and
 * @f$j@f$ specifies which one. It can also be used to compare in different
 * ways, like comparing the birth date including the birth year to know who
 * is younger or comparing the birth date without the year to know who has
 * his next birthday first.
 *
 * For details on how to use these comparisons in algorithms that require an
 * ordinary Java  Comparator, such as the `sort` methods in class `Arrays`,
 * see  @ref MultiDimComparator.
 *
 * <div class="SSJ-bigskip"></div>
 */
public interface MultiDimComparable<T> {

   /**
    * This method returns the number @f$d@f$ of coordinates for which the
    * method  #compareTo can be invoked for this object. That is,
    * `compareTo (other, j)` can be called only for
    * @f$\mathtt{j}=0,…,\mathtt{dimension()-1}@f$.
    */
   public int getStateDimension();

   /**
    * Similar to  #compareTo in the class  Comparable, except that one
    * specifies the index @f$j@f$ of the criterion on which the objects
    * are to be compared. It compares objects of type `T` based on
    * coordinate @f$j@f$. The method must return a negative integer, zero,
    * or a positive integer depending on whether the current object is
    * less than, equal to, or greater than `o`, based on coordinate
    * @f$j@f$. If @f$j@f$ is outside the range @f$\{0,…,d-1\}@f$, this
    * method should throw an  IllegalArgumentException.
    */
   public int compareTo (T other, int j);

}
