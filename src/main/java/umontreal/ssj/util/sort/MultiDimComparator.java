/*
 * Class:        MultiDimComparator
 * Description:  
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
import java.lang.IllegalArgumentException;

/**
 * This class is useful if one wishes to perform an ordinary one-dimensional
 * sort on  @ref MultiDimComparable<T> objects based on a single coordinate
 * @f$j@f$, which is specified in the constructor. It defines a bridge
 * between the  @ref MultiDimComparable<T> interface and the classic Comparator
 * in Java. It implements  Comparator in a way that the method `compare(o1,
 * o2)` compares two  @ref MultiDimComparable<T> objects in the dimension
 * @f$j@f$ given in the constructor, by calling `o1.compareTo(o2, j)`.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class MultiDimComparator<T extends MultiDimComparable<? super T>>
                                  implements Comparator<T> {
   private int compareDim;

   /**
    * Constructs a comparator that uses coordinate `j` for the comparison
    * of  @ref MultiDimComparable<T> objects. One must have @f$j
    * \in\{0,…,d-1\}@f$.
    *  @param j            index used for comparison
    */
   public MultiDimComparator (int j) {
      compareDim = j;
   }

   /**
    * Calls `o1.compareTo(o2, j)` from class  @ref MultiDimComparable<T>.
    *  @param o1           first object to compare
    *  @param o2           second object to compare
    */
   public int compare (T o1, T o2) {
      if (compareDim >= o1.getStateDimension() || compareDim >= o2.getStateDimension())
         throw new IllegalArgumentException("Comparing in a "+
                    "dimension larger than object dimension");
      return o1.compareTo (o2, compareDim);
   }

 }
