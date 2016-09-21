/*
 * Class:        DoubleArrayComparator
 * Description:  Compares two double's arrays
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package umontreal.ssj.util.sort;
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
