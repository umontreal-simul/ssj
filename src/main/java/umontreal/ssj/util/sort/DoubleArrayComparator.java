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
