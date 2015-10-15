/*
 * Class:        BakerTransformedPointSet
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
package umontreal.ssj.hups;

import umontreal.ssj.util.PrintfFormat;

/**
 * This container class embodies a point set to which a *Baker
 * transformation* is applied (see, e.g., @cite vHIC02a&thinsp;). It
 * transforms each coordinate @f$u@f$ into @f$2u@f$ if @f$u \le1/2@f$ and
 * @f$2(1-u)@f$ if @f$u > 1/2@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class BakerTransformedPointSet extends ContainerPointSet {

   /**
    * Constructs a Baker-transformed point set from the given point set
    * `P`.
    *  @param P            point set for which we want a Baker-transfomed
    *                      version
    */
   public BakerTransformedPointSet (PointSet P) {
      init (P);
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