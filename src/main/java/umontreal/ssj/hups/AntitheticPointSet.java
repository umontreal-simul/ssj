/*
 * Class:        AntitheticPointSet
 * Description:  
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001--2018  Pierre L'Ecuyer and Universite de Montreal
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
package umontreal.ssj.hups;

import umontreal.ssj.util.PrintfFormat;

/**
 * This container class provides antithetic versions of the contained points.   
 * That is, @f$1 - u_{i,j}@f$ is returned in place of coordinate @f$u_{i,j}@f$. To generate
 * regular and antithetic variates with a point set `p`, e.g., for variance
 * reduction, one can define an  @ref AntitheticPointSet object `pa` that
 * contains `p`, and then generate the regular variates with `p` and the
 * antithetic variates with `pa`.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class AntitheticPointSet extends ContainerPointSet {

   /**
    * Constructs an antithetic point set from the given point set `p`.
    *  @param p            point set for which we want an antithetic version
    */
   public AntitheticPointSet (PointSet p) {
      init (p);
   }


   public double getCoordinate (int i, int j) {
      return 1.0 - P.getCoordinate (i, j);
   }

   public PointSetIterator iterator(){
      return new AntitheticPointSetIterator();
   }

   public String toString() {
      return "Antithetic point set of: {" + PrintfFormat.NEWLINE +
              P.toString() + PrintfFormat.NEWLINE + "}";
   }


   // ***************************************************************

   protected class AntitheticPointSetIterator
                   extends ContainerPointSetIterator {

      public double nextCoordinate() {
         return 1.0 - innerIterator.nextCoordinate();
      }

      public double nextDouble() {
         return 1.0 - innerIterator.nextCoordinate();
      }

      public void nextCoordinates (double p[], int d)  {
         innerIterator.nextCoordinates (p, d);
         for (int j = 0; j < d; j++)
            p[j] = 1.0 - p[j];
      }

      public int nextPoint (double p[], int d)  {
         innerIterator.nextPoint (p, d);
         for (int j = 0; j < d; j++)
            p[j] = 1.0 - p[j];
         return getCurPointIndex();
      }

   }
}