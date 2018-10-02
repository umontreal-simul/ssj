/*
 * Class:        BakerTransformedPointSet
 * Description:  
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001-2018  Pierre L'Ecuyer and Universite de Montreal
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
 * This container class embodies a point set to which a *baker's
 * transformation* (also called a *tent transform*) is applied 
 * (see, e.g., @cite rDIC10a, @cite vHIC02a, @cite vLEC09f). 
 * It transforms each coordinate @f$u@f$ of each point into @f$2u@f$ if @f$u \le1/2@f$ and
 * @f$2(1-u)@f$ if @f$u > 1/2@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class BakerTransformedPointSet extends ContainerPointSet {

   /**
    * Constructs a baker-transformed point set from the given point set
    * `p`.
    *  @param p   point set for which we want a baker-transformed version
    */
   public BakerTransformedPointSet (PointSet p) {
      init (p);
   }


   public double getCoordinate (int i, int j) {
      double u = P.getCoordinate (i, j);
      if (u < 0.5) return 2.0 * u;
      else return 2.0 * (1 - u);
   }

   public PointSetIterator iterator(){
      return new BakerTransformedPointSetIterator();
   }

   public String toString() {
      return "Baker transformed point set of: {" + PrintfFormat.NEWLINE
              + P.toString() + PrintfFormat.NEWLINE + "}";
   }

/*
   public String formatPoints() {
      try {
         return super.formatPoints();
      }
      catch (UnsupportedOperationException e) {
         return "The values are Baker transformed for each coordinate:" +
                 PrintfFormat.NEWLINE + " {" +
                 P.formatPoints() + PrintfFormat.NEWLINE + "}";
      }
   }
*/
   // ***************************************************************

   protected class BakerTransformedPointSetIterator
                   extends ContainerPointSetIterator {

      public double nextCoordinate() {
         double u = innerIterator.nextCoordinate();
         if (u < 0.5) return 2.0 * u;
         else return 2.0 * (1.0 - u);
      }

      // Same as nextCoordinate.
      public double nextDouble() {
         double u = innerIterator.nextCoordinate();
         if (u < 0.5) return 2.0 * u;
         else return 2.0 * (1.0 - u);
      }

      public void nextCoordinates (double p[], int d)  {
         innerIterator.nextCoordinates (p, d);
         for (int j = 0; j < d; j++)
            if (p[j] < 0.5) p[j] *= 2.0;
            else p[j] = 2.0 * (1.0 - p[j]);
      }

      public int nextPoint (double p[], int d)  {
         innerIterator.nextPoint (p, d);
         for (int j = 0; j < d; j++)
            if (p[j] < 0.5) p[j] *= 2.0;
            else p[j] = 2.0 * (1.0 - p[j]);
         return getCurPointIndex();
      }

   }
}