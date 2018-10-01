/*
 * Class:        PaddedPointSet
 * Description:  container class
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
package umontreal.ssj.hups;

import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.rng.RandomStream;

/**
 * This container class realizes *padded point sets*, constructed by taking
 * some coordinates from a point set @f$P_1@f$, other coordinates from a
 * point set @f$P_2@f$, and so on. This can be used to implement *Latin
 * Supercube Sampling* as proposed in @cite vOWE98a, for example. After calling the
 * constructor to create the structure, component point sets can be padded to
 * it by calling  `padPointSet` or  `padPointSetPermute`.
 *
 * Only sets with the same number of points can be padded. Point sets with
 * too many points or coordinates can be trimmed down by using the class
 * @ref SubsetOfPointSet before they are padded. Infinite-dimensional point
 * sets are allowed, but once one is padded, no additional point set can be padded.
 *
 * The points of each padded set can be permuted randomly, independently
 * across the padded sets. If such a random permutation is desired, the point
 * set should be padded via  `padPointSetPermute'. When calling  `randomize`,
 * random permutations are then generated for all point sets that have been padded
 * by  `padPointSetPermute`.  This is all what this `randomize` method does.
 * If other randomizations are desired, the must be applied separately.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class PaddedPointSet extends PointSet {
   protected int curPointSets = 0;     // Current number of padded point sets.
   protected int maxPointSets;         // Max. number of padded point sets.
   protected PointSet pointSet[];      // List of padded point sets
   protected int startDim[];           // Starting dim. for padded points sets.
   protected int permutation[][];      // One permutation for each point set.

   /**
    * Constructs a structure for padding at most `maxPointSets` point
    * sets. This structure is initially empty and will eventually contain
    * the different point sets that are padded, via the method `padPointSet`.
    *  @param maxPointSets maximum number of point sets authorized by the
    *  constructed object.
    */
   public PaddedPointSet (int maxPointSets) {
      this.maxPointSets = maxPointSets;
      pointSet = new PointSet[maxPointSets];
      startDim = new int[maxPointSets];
      permutation = new int[maxPointSets][];
   }

   /**
    * Pads the point set `P` to the present structure.  
    * This adds new coordinates to the points.
    *  @param P            point set being padded
    */
   public void padPointSet (PointSet P) {
      if (curPointSets == maxPointSets)
         throw new IllegalArgumentException
            ("Cannot pad more, increase maxPointSets parameter");

      if (dim == Integer.MAX_VALUE)
         throw new IllegalArgumentException
            ("Cannot pad more, dimension already infinite");

      if (curPointSets > 0 && numPoints != P.getNumPoints())
         throw new IllegalArgumentException
            ("Padded points must have same number of points");

      if (curPointSets == 0)
          numPoints = P.getNumPoints();
      
      if (P.getDimension() == Integer.MAX_VALUE)
          dim = Integer.MAX_VALUE;
      else
          dim += P.getDimension();

      pointSet[curPointSets] = P;
      startDim[curPointSets] = dim;
      ++curPointSets;
   }

   /**
    * Pads the point set `P`, which is assumed to be *finite*. A new random
    * permutation of these points will be generated when calling  `randomize` and 
    * will be used to access the coordinates taken from the points of `P` (i.e., these
    * points are randomly permuted).
    *  @param P            point set being padded
    */
   public void padPointSetPermute (PointSet P) {
      if (curPointSets == 0)
          numPoints = P.getNumPoints();
      if (numPoints == Integer.MAX_VALUE)
         throw new IllegalArgumentException
            ("Cannot generate infinite permutation");
      permutation[curPointSets] = new int[numPoints];
      for (int i = 0; i < numPoints; i++)
         permutation[curPointSets][i] = i;
      padPointSet (P);
   }


   /**
    *  @returns  Coordinate @f$u_{i,j}@f$ for the padded point set.
    */
   public double getCoordinate (int i, int j) {
      int set = 0;  // Will be the padded set to which this coordinate belongs,
      if (j >= dim)
         throw new IllegalArgumentException ("Not enough dimensions");
      while (j >= startDim[set])
         set++;
      /*
      if (permutation[set] == null)
         pointSet[set].resetPoint(i);
      else
         pointSet[set].resetPoint(permutation[set][i]);

      if (set == 0)
         pointSet[0].resetCoordinate (j);
      else
         pointSet[set].resetCoordinate (j - startDim[set - 1]);
      
      return pointSet[set].nextCoordinate();
      */
      if (permutation[set] != null)
         i = permutation[set][i];
      if (set != 0)
         j = j - startDim[set - 1];
      return pointSet[set].getCoordinate (i, j);
   }

   /**
    *  Erases all the random permutations of the padded point sets that are randomly permuted.
    *  All the permutations are reset to the identity.
    */
   public void unrandomize() {
      for (int set = 0; set < curPointSets; set++) {
         if (permutation[set] != null) {
            for (int i = 0; i < numPoints; i++)
               permutation[set][i] = i;
         }
      }
   }

   /**
    *  Randomly permutes the points for the padded point sets that have been added via
    *  `padPointSetPermute`.
    */
   public void randomize (RandomStream stream)  {
      // Executes the randomizations of the list
      //  super.randomize (stream);  // Removed this because PointSet.randomize(stream) does nothing.
      /* can also use lazy permutations */
      for (int set = 0; set < curPointSets; set++)
         if (permutation[set] != null) {
             for (int i = 0; i < numPoints - 1; i++) {
                 int u = stream.nextInt(0, numPoints - i - 1);
                 int h = permutation[set][i];
                 permutation[set][i] = permutation[set][i + u];
                 permutation[set][i + u] = h;
             }
         }
   }

   /**
    *  Returns a `PaddedIterator` for this padded point set.
    */
   public PointSetIterator iterator() {
      return new PaddedIterator();
   }

   public String toString() {
      StringBuffer sb = new StringBuffer ("Padded point set" +
                                           PrintfFormat.NEWLINE);
      sb.append ("Maximal number of point sets: " + maxPointSets +
                  PrintfFormat.NEWLINE);
      sb.append ("Current number of point sets: " + curPointSets +
                  PrintfFormat.NEWLINE);
      sb.append ("Number of points: " + numPoints + PrintfFormat.NEWLINE);
      for (int i = 0; i < curPointSets; i++) {
         if (i != 0)
            sb.append (PrintfFormat.NEWLINE);
         if (permutation[i] == null)
            sb.append ("Point set ");
         else
            sb.append ("Permuted point set ");
         sb.append ( i + " information: {" + PrintfFormat.NEWLINE 
                + pointSet[i].toString() + PrintfFormat.NEWLINE + "}");
      }
      return sb.toString();
   }


   // ************************************************************

   private class PaddedIterator extends DefaultPointSetIterator {

      private PointSetIterator[] pointSetIterators; // One for each padded set.
      private int currentSet = 0;
      private double[] temp;

      public PaddedIterator() {
         pointSetIterators = new PointSetIterator[curPointSets];
         int maxdim = 0;
         for (int i = 0; i < curPointSets; i++) {
            pointSetIterators[i] = pointSet[i].iterator();
            if (pointSet[i].getDimension() > maxdim)
               maxdim = pointSet[i].getDimension();
            if (permutation[i] != null)
               pointSetIterators[i].setCurPointIndex (permutation[i][0]);
         }
         if (maxdim == Integer.MAX_VALUE)
            temp = new double[16];
         else
            temp = new double[maxdim];
      }

      public void setCurCoordIndex (int j) {
         int set = 0;
         if (j >= dim)
            throw new IllegalArgumentException ("Not enough dimensions");
         while (j >= startDim[set])
            set++;
         currentSet = set;
         pointSetIterators[currentSet].setCurCoordIndex 
            (set == 0 ? j : j - startDim[set-1]);
         for (set = currentSet+1; set < pointSetIterators.length; set++)
            pointSetIterators[set].resetCurCoordIndex();
         curCoordIndex = j;
      }

      public void resetCurCoordIndex() {
         currentSet = 0;
         for (int i = 0; i < pointSetIterators.length; i++)
            pointSetIterators[i].resetCurCoordIndex();
         curCoordIndex = 0;
      }

      public double nextCoordinate() {
         if (curPointIndex >= numPoints || curCoordIndex >= dim)
            outOfBounds();
         if (curCoordIndex >= startDim[currentSet])
            currentSet++;
         double coord = pointSetIterators[currentSet].nextCoordinate();
         curCoordIndex++;
         return coord;
      }

      public void nextCoordinates (double[] p, int d) {
         if (curPointIndex >= numPoints || d > dim)
            outOfBounds();
         int i = 0;
         while (i < d) {
            int dimen = pointSet[currentSet].getDimension();
            if (dimen == Integer.MAX_VALUE)
               dimen = d - i;
            else
               dimen -= pointSetIterators[currentSet].getCurCoordIndex();
            pointSetIterators[currentSet].nextCoordinates (temp, dimen);
            System.arraycopy (temp, 0, p, i, dimen);
            i += dimen;
            curCoordIndex += dimen;
            if (i < d)
               currentSet++;
         }
      }

      public void setCurPointIndex (int i) {
         for (int it = 0; it < pointSetIterators.length; it++)
            pointSetIterators[it].setCurPointIndex 
               (permutation[it] == null ? i : permutation[it][i]);
         curPointIndex = i;  
         curCoordIndex = 0;
         currentSet = 0;
      }

      public void resetCurPointIndex() { 
         for (int i = 0; i < pointSetIterators.length; i++) {
            if (permutation[i] == null)
               pointSetIterators[i].resetCurPointIndex();
            else
               pointSetIterators[i].setCurPointIndex (permutation[i][0]);
         }
         curPointIndex = 0;  
         curCoordIndex = 0;
         currentSet = 0;
      }

      public int resetToNextPoint() {
         for (int i = 0; i < pointSetIterators.length; i++) {
            if (permutation[i] == null)
               pointSetIterators[i].resetToNextPoint();
            else
               pointSetIterators[i].setCurPointIndex 
                  (permutation[i][curPointIndex+1]);
         }
         currentSet = 0;
         curCoordIndex = 0;
         return ++curPointIndex;
      }

      public String formatState() {
         return super.formatState() + PrintfFormat.NEWLINE + 
                     "Current padded set: " + currentSet;
      }
   }
}
