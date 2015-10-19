/*
 * Class:        DoubleArrayComparator
 * Description:  Compares two double's arrays
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
package umontreal.ssj.util;
import java.util.Comparator;

/**
 * This provides an implementation of  Comparator in which arrays of `double`
 * in @f$d@f$ dimensions are compared by comparing their coordinate @f$j@f$
 * in the natural order of real numbers, where @f$j \in\{0,…,d-1\}@f$ is
 * given in the constructor. The method `compare(d1, d2)` returns @f$-1@f$,
 * @f$0@f$, or @f$1@f$ depending on whether `d1[j]` is less than, equal to,
 * or greater than `d2[j]`.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class DoubleArrayComparator implements Comparator<double[]> {
   private int j;

   /**
    * Constructs a comparator, where `j` is the coordinate used for the
    * comparisons.
    *  @param j            index used for comparison
    */
   public DoubleArrayComparator (int j) {
      this.j = j;
   }

   /**
    * Returns @f$-1@f$, @f$0@f$, or @f$1@f$ depending on whether `d1[j]`
    * is less than, equal to, or greater than `d2[j]`.
    *  @param d1           first array
    *  @param d2           second array
    */
   public int compare (double[] d1, double[] d2) {
      if (j >= d1.length || j >= d2.length)
         throw new IllegalArgumentException("Comparing in a"+
                    "dimension larger than array dimension");
      return (d1[j]< d2[j] ? -1 : (d1[j] > d2[j] ? 1 : 0));
   }

}