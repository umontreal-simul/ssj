/*
 * Class:        SubsetOfPointSet
 * Description:  Subset of a point set
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

import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.util.PrintfFormat;

 /*Attention: No array index range tests neither for the dimension
   nor for the number of points is performed. This is left to JAVA. */

/**
 * This container class permits one to select a subset of a point set. This
 * is done by selecting a range or providing an array of either point or
 * coordinate indices. A typical application of a range selection is to make
 * the number of points or the dimension finite. It is also possible to
 * provide, for example, a random permutation in the selection of components.
 * It is possible also to take *projections* of coordinates for selected
 * dimensions.
 *
 * Selecting a new subset of points or coordinates overwrites the previous
 * selection. The specification of a subset with respect to the points is
 * independent from selecting a subset with respect to the coordinates. The
 * number of points and the dimension are adapted to the current selection
 * and all indices still start from 0, i.e., the subset works just like an
 * ordinary point set.
 *
 * When the points or coordinates ranges are changed, existing iterators
 * become invalid. They should be reconstructed or reset to avoid
 * inconsistencies.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 *
 * ## Additional Information
 *
 * @remark **Pierre:** Code à reviser.
 */
public class SubsetOfPointSet extends PointSet {
   protected PointSet P;                  // Source points
   protected int i_from, i_to, i_index[]; // Limits or lookup for row
   protected int j_from, j_to, j_index[]; // Limits or lookup for column

   /**
    * Constructs a new  @ref PointSet object, initially identical to `p`,
    * and from which a subset of the points and/or a subset of the
    * coordinates is to be extracted.
    *  @param p            point set for which a subset is constructed
    */
   public SubsetOfPointSet (PointSet p) {
      this.P = p;
      numPoints = p.getNumPoints();
      dim = p.getDimension();
      i_from = 0;
      i_to = p.getNumPoints();
      j_from = 0;
      j_to = p.getDimension();
   }

   /**
    * Selects the points numbered from "<tt>from</tt>" to "<tt>to -
    * 1</tt>" from the original point set.
    *  @param from         index of the point, in the contained point set,
    *                      corresponding to the point 0 of this point set
    *  @param to           index of the last point taken from the
    *                      contained point set
    */
   public void selectPointsRange (int from, int to) {
      if (0 > from || from >= to || to > P.getNumPoints())
         throw new IllegalArgumentException ("Invalid range for points");

      i_index = null;
      i_from = from;
      i_to = to;
      numPoints = to - from;
   }

   /**
    * Selects the `numPoints` points whose numbers are provided in the
    * array `pointIndices`.
    *  @param pointIndices array of point indices to be selected
    *  @param numPoints    number of points in the subset of point set
    */
   public void selectPoints (int[] pointIndices, int numPoints) {
      if (numPoints > P.getNumPoints() || numPoints > pointIndices.length)
         throw new IllegalArgumentException ("Number of indices too large");
      i_index = pointIndices;
      this.numPoints = numPoints;
   }

   /**
    * Selects the coordinates from "<tt>from</tt>" to "<tt>to - 1</tt>"
    * from the original point set.
    *  @param from         index of the coordinate, in the contained point
    *                      set, corresponding to the coordinate 0 of each
    *                      point of this point set
    *  @param to           index of the last coordinate taken for each
    *                      point from the contained point set
    */
   public void selectCoordinatesRange (int from, int to) {
      if (0 > from || from >= to || to > P.getDimension())
         throw new IllegalArgumentException("Invalid column range");
      j_index = null;
      j_from = from;
      j_to = to;
      dim = to - from;
   }

   /**
    * Selects the `numCoord` coordinates whose numbers are provided in the
    * array `coordIndices`.
    *  @param coordIndices array of coordinate indices to be selected
    *  @param numCoord     number of coordinatess for each point in the
    *                      subset of point set
    */
   public void selectCoordinates (int[] coordIndices, int numCoord) {
      if (numCoord > P.getDimension() || numCoord > coordIndices.length)
         throw new IllegalArgumentException ("Number of indices too large");

      j_index = coordIndices;
      this.dim = numCoord;
   }


   public double getCoordinate (int i, int j) {
      int access_i, access_j;

      // if no range check done: left to JAVA array index check
      
      if (i_index == null) {
         if (i < 0 || i >= numPoints)
            throw new IllegalArgumentException ("Row out of range");
         
         access_i = i + i_from;
      } else 
         access_i = i_index[i];

      if (j_index == null) {
         if (j < 0 || j > dim)
            throw new IllegalArgumentException("Column out of range");
         
         access_j = j + j_from;
      } else 
         access_j = j_index[j];

      return P.getCoordinate (access_i, access_j);
   }

   public PointSetIterator iterator() {
      return new SubsetIterator();
   }

   public String toString() {
      StringBuffer sb = new StringBuffer ("Subset of point set" +
              PrintfFormat.NEWLINE);
      sb.append ("Inner point set information {" + PrintfFormat.NEWLINE);
      sb.append (P.toString());
      sb.append (PrintfFormat.NEWLINE + "}" + PrintfFormat.NEWLINE);

      if (i_index == null)
         sb.append ("Points range from " + i_from + " to " + i_to + "." +
               PrintfFormat.NEWLINE);
      else {
         sb.append ("Point indices: [");
         boolean first = true;
         for (int i = 0; i < numPoints; i++) {
            if (first)
               first = false;
            else
               sb.append (", ");
            sb.append (i_index[i]);
         }
         sb.append ("]" + PrintfFormat.NEWLINE);
      }

      if (j_index == null)
         sb.append ("Coordinates range from " + j_from + " to " + j_to + ".");
      else {
         sb.append ("Coordinate indices: [");
         boolean first = true;
         for (int i = 0; i < dim; i++) {
            if (first)
               first = false;
            else
               sb.append (", ");
            sb.append (j_index[i]);
         }
         sb.append ("]");
      }

      return sb.toString();
   }


   // ***********************************************************

   private class SubsetIterator extends DefaultPointSetIterator {

      private PointSetIterator innerIterator;
/*
      private int i_from;
      private int i_to;
      private int j_from;
      private int j_to;
      private int[] i_index;
      private int[] j_index;
*/

      SubsetIterator() {
         // Since one can change range after construction, we
         // must save the current one.
         //this.i_from = SubsetOfPointSet.this.i_from;
         //this.i_to = SubsetOfPointSet.this.i_to;
         //this.j_from = SubsetOfPointSet.this.j_from;
         //this.j_to = SubsetOfPointSet.this.j_to;

         // Also recopy indices in case one has set the indices,
         // kept the array and modified it after construction.
         //if (SubsetOfPointSet.this.i_index == null)
         //   this.i_index = null;
         //else {
         //   this.i_index = new int[SubsetOfPointSet.this.i_index.length];
         //   System.arraycopy (SubsetOfPointSet.this.i_index, 0,
         //             this.i_index, 0, numPoints);
         //}

         //if (SubsetOfPointSet.this.j_index == null)
         //   this.j_index = null;
         //else {
         //   this.j_index = new int[SubsetOfPointSet.this.j_index.length];
         //   System.arraycopy (SubsetOfPointSet.this.j_index, 0, this.j_index, 0, dim);
         //}

         // Create the inner iterator and set its state according to the subset.
         innerIterator = P.iterator();
         if (i_index == null) {
            if (i_from != 0)
               innerIterator.setCurPointIndex (i_from);
         }
         else {
            if (i_index[0] != 0)
               innerIterator.setCurPointIndex (i_index[0]);
         }

         if (j_index == null) {
            if (j_from != 0)
               innerIterator.setCurCoordIndex (j_from);
         }
         else {
            if (j_index[0] != 0)
               innerIterator.setCurCoordIndex (j_index[0]);
         }
      }

      public void setCurCoordIndex (int j) {
         if (j_index == null)
            innerIterator.setCurCoordIndex (j + j_from);
         else
            innerIterator.setCurCoordIndex (j_index[j]);
         curCoordIndex = j;
      }

      public void resetCurCoordIndex() {
         if (j_index == null) {
            if (j_from == 0)
               innerIterator.resetCurCoordIndex();
            else
               innerIterator.setCurCoordIndex (j_from);
         }
         else {
            if (j_index[0] == 0)
               innerIterator.resetCurCoordIndex();
            else
               innerIterator.setCurCoordIndex (j_index[0]);
         }
         curCoordIndex = 0;
      }

      public double nextCoordinate() {
         if ((curPointIndex >= numPoints) || (curCoordIndex >= dim))
            outOfBounds();
         // The inner iterator could throw an exception.
         // If that happens, e must not alter the current coordinate.
         double coord = 0.0;

         if (j_index == null)
            coord = innerIterator.nextCoordinate();
         else {
            int currentIndex = j_index[curCoordIndex];
            int futureIndex = (curCoordIndex+1) == dim ? 
                              currentIndex+1 : j_index[curCoordIndex+1];
            coord = innerIterator.nextCoordinate();
            if (futureIndex != (currentIndex+1))
               innerIterator.setCurCoordIndex (futureIndex);
         }
         curCoordIndex++;
         return coord;
      }

      public void nextCoordinates (double[] p, int d) {
         if (curPointIndex >= numPoints || curCoordIndex + d > dim)
            outOfBounds();
         if (j_index != null) {
            super.nextCoordinates (p, d);
            return;
         }
         innerIterator.nextCoordinates (p, d);
         curCoordIndex += d;
      }

      public void setCurPointIndex (int i) {
         if (i_index == null)
            innerIterator.setCurPointIndex (i + i_from);
         else
            innerIterator.setCurPointIndex (i_index[i]);
         curPointIndex = i;
         resetCurCoordIndex();
      }

      public void resetCurPointIndex() { 
         if (i_index == null) {
            if (i_from == 0)
               innerIterator.resetCurPointIndex();
            else
               innerIterator.setCurPointIndex (i_from);
         }
         else {
            if (i_index[0] == 0)
               innerIterator.resetCurPointIndex();
            else
               innerIterator.setCurPointIndex (i_index[0]);
         }
         curPointIndex = 0;
         resetCurCoordIndex();
      }

      public int resetToNextPoint() {
         if (i_index == null)
            innerIterator.resetToNextPoint();
         else if (curPointIndex < (numPoints-1))
            innerIterator.setCurPointIndex (i_index[curPointIndex + 1]);
         curPointIndex++;  
         resetCurCoordIndex();
         return curPointIndex;
      }
   }
}