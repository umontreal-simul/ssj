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
 * An implementation of  Comparator which compares two `double` arrays by
 * comparing their <tt>i</tt>-<em>th</em> element, where `i` is given in the
 * constructor. Method `compare(d1, d2)` returns @f$-1@f$, @f$0@f$, or
 * @f$1@f$ depending on whether `d1[i]` is less than, equal to, or greater
 * than `d2[i]`.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class DoubleArrayComparator implements Comparator<double[]> {
   private int i;

   /**
    * Constructs a comparator, where `i` is the index used for the
    * comparisons.
    *  @param i            index used for comparison
    */
   public DoubleArrayComparator (int i) {
      this.i = i;
   }

   /**
    * Returns @f$-1@f$, @f$0@f$, or @f$1@f$ depending on whether `d1[i]`
    * is less than, equal to, or greater than `d2[i]`.
    *  @param d1           first array
    *  @param d2           second array
    */
   public int compare (double[] d1, double[] d2) {
      if (i >= d1.length || i >= d2.length)
         throw new IllegalArgumentException("Comparing in a"+
                    "dimension larger than arary dimension");
      return (d1[i]< d2[i] ? -1 : (d1[i] > d2[i] ? 1 : 0));
   }

}