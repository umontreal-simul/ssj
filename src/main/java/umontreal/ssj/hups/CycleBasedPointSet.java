/*
 * Class:        CycleBasedPointSet
 * Description:  provides the basic structures for storing and manipulating
                 a highly uniform point set defined by a set of cycles
 defined by a set of cycles.
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

import umontreal.ssj.util.*;
import umontreal.ssj.rng.RandomStream;
import cern.colt.list.*;

/**
 * This abstract class provides the basic structures for storing and
 * manipulating a point set defined by a set of cycles. The @f$s@f$-dimensional 
 * points are all the vectors of @f$s@f$ successive
 * values found in any of the cycles, from any starting point. Since this is
 * defined for any positive integer @f$s@f$, the points effectively have an
 * infinite number of dimensions. The number of points, @f$n@f$, is the sum
 * of lengths of all the cycles. The cycles of the point set are simply
 * stored as a list of arrays, where each array contains the successive
 * values for a given cycle. By default, the values are stored in `double`.
 *
 * This structure is convenient for implementing recurrence-based point sets,
 * where the point set in @f$s@f$ dimensions is defined as the set of all
 * vectors of @f$s@f$ successive values of a periodic recurrence, from all
 * its possible initial states.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public abstract class CycleBasedPointSet extends PointSet {

   protected int numCycles = 0;     // Total number of cycles.
   // dim = Integer.MAX_VALUE;      // Dimension is infinite.
   private double[] shift;          // Random shift, initially null.
                                    // Entry j is for dimension j.
   protected ObjectArrayList cycles = new ObjectArrayList(); // List of cycles.


   public double getCoordinate (int i, int j) {
      // Find cycle that contains point i, then index in cycle.
      int l = 0;         // Length of next cycle.
      int n = 0;         // Total length of cycles added so far.
      int k;
      for (k = 0;  n <= i;  k++)
         n += l = ((AbstractList) cycles.get (k)).size();
      AbstractList curCycle = (AbstractList) cycles.get (k-1);
      int coordinate = (i - n + l + j) % curCycle.size();
//      double[] curCycleD = ((DoubleArrayList) curCycle).elements();
//      return curCycleD[coordinate];
      double x = ((DoubleArrayList) curCycle).get (coordinate);
      return x;
   }

/**
 * Same as the same method in `PointSet`.
 *  @param stream       Stream used to generate random numbers
 */
public void addRandomShift (int d1, int d2, RandomStream stream) {
      if (null == stream)
         throw new IllegalArgumentException (
              PrintfFormat.NEWLINE +
              "   Calling addRandomShift with null stream");
      if (0 == d2)
         d2 = Math.max (dim, 1);
      if (shift == null) {
         shift = new double[d2];
         capacityShift = d2;
      } else if (d2 > capacityShift) {
         int d3 = Math.max (4, capacityShift);
         while (d2 > d3)
            d3 *= 2;
         double[] temp = new double[d3];
         capacityShift = d3;
         for (int i = 0; i < d1; i++)
            temp[i] = shift[i];
         shift = temp;
      }
      dimShift = d2;
      for (int i = d1; i < d2; i++)
         shift[i] = stream.nextDouble ();
      shiftStream = stream;
   }

   /**
    * Erases the current random shift, if any.
    */
   public void clearRandomShift() {
      super.clearRandomShift();
      shift = null;
   }

   /**
    * Adds the cycle `c` to the list of all cycles. This method is used by
    * subclass constructors to fill up the list of cycles.
    */
   protected void addCycle (AbstractList c) {
      // Adds the cycle `c` to the list of all cycles.
      // Used by subclass constructors to fill up the list of cycles.
      cycles.add (c);
      numCycles++;
      numPoints += c.size();
   }


   public int getDimension() {
     return Integer.MAX_VALUE;
   }

   public PointSetIterator iterator(){
      return new  CycleBasedPointSetIterator();
   }

   public String toString() {
      String s = super.toString();
      return s + PrintfFormat.NEWLINE + "Number of cycles: " + numCycles;
   }

   public String formatPoints() {
      StringBuffer sb = new StringBuffer (toString());
      for (int c = 0; c < numCycles; c++) {
         AbstractList curCycle = (AbstractList)cycles.get (c);
         double[] cycle = ((DoubleArrayList)curCycle).elements();
         sb.append (PrintfFormat.NEWLINE + "Cycle " + c + ": (");
         boolean first = true;
         for (int e = 0; e < curCycle.size(); e++) {
            if (first)
               first = false;
            else
               sb.append (", ");
            sb.append (cycle[e]);
         }
         sb.append (")");
      }
      return sb.toString();
   }

   // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   public class CycleBasedPointSetIterator extends DefaultPointSetIterator {

      protected int startPointInCycle = 0;   // Index where the current point
                                             // starts in the current cycle.
      protected int curCoordInCycle = 0;     // Index of the current coordinate
                                             // in the current cycle.
      protected int curCycleIndex = 0;       // Index of the current cycle.
      protected AbstractList curCycle;       // The current cycle.
      protected double[] curCycleD;          // The array for current cycle


      public CycleBasedPointSetIterator () {
         init();
      }

      protected void init () {
         resetCurCycle(0);
      }

      public void resetCurCycle (int index) {
         curCycleIndex = index;
         curCycle = (AbstractList) cycles.get (index);
         curCycleD = ((DoubleArrayList) curCycle).elements();
      }

      public void setCurCoordIndex (int i) {
         curCoordIndex = i;
         curCoordInCycle = (i + startPointInCycle) % curCycle.size();
      }

      public void resetCurCoordIndex() {
         curCoordIndex = 0;
         curCoordInCycle = startPointInCycle;
      }

      public boolean hasNextCoordinate() {
         return true;
      }

      // We want to avoid generating 0 or 1
      public double nextDouble() {
         return nextCoordinate() + EpsilonHalf;
      }

      public double nextCoordinate() {
         // First, verify if there are still points....
         if (getCurPointIndex() >= getNumPoints ())
            outOfBounds();
         double x = curCycleD [curCoordInCycle];
         if (shift != null) {
             if (curCoordIndex >= dimShift)   // Extend the shift.
                addRandomShift (dimShift, curCoordIndex + 1, shiftStream);
             x += shift[curCoordIndex];
             if (x >= 1.0)
                x -= 1.0;
             if (x <= 0.0)
                x = EpsilonHalf;  // avoid x = 0
         }
         curCoordIndex++;
         curCoordInCycle++;
         if (curCoordInCycle >= curCycle.size())
            curCoordInCycle = 0;
         return x;
      }

      public void nextCoordinates (double p[], int dim) {
         // First, verify if there are still points....
         if (getCurPointIndex() >= getNumPoints ())
            outOfBounds();
         if (curCoordIndex + dim >= dimShift)
            addRandomShift (dimShift, curCoordIndex + dim + 1, shiftStream);
         // int j = curCoordInCycle;
         int maxj = curCycle.size();
         double x;
         for (int i = 0; i < dim; i++) {
            x = curCycleD [curCoordInCycle++];
            if (curCoordInCycle >= maxj) curCoordInCycle = 0;
            if (shift != null) {
               x += shift[curCoordIndex + i];
               if (x >= 1.0)
                  x -= 1.0;
               if (x <= 0.0)
                  x = EpsilonHalf;  // avoid x = 0
           }
            p[i] = x;
         }
         curCoordIndex += dim;
      }

      public void setCurPointIndex (int i) {
         int l = 0;
         int n = 0;
         int j ;
         for (j=0;  n <= i;  j++)
            n += l = ((AbstractList) cycles.get (j)).size();
         resetCurCycle (j-1);
         startPointInCycle = i - n + l;
         curPointIndex = i;
         curCoordIndex = 0;
         curCoordInCycle = startPointInCycle;
      }

      public void resetCurPointIndex() {
         resetCurCycle (0);
         startPointInCycle = 0;
         curPointIndex = 0;
         curCoordIndex = 0;
         curCoordInCycle = 0;
      }

      public int resetToNextPoint() {
         curPointIndex++;
         startPointInCycle++;
         if (startPointInCycle >= curCycle.size()) {
            startPointInCycle = 0;
            if (curCycleIndex < (numCycles - 1))
               resetCurCycle (curCycleIndex + 1);
         }
         curCoordIndex = 0;
         curCoordInCycle = startPointInCycle;
         return curPointIndex;
      }

      public int nextPoint (double p[], int dim) {
         // First, verify if there are still points....
         if (getCurPointIndex() >= getNumPoints ())
            outOfBounds();
         int j = startPointInCycle;
         int maxj = curCycle.size() - 1;
         for (int i = 0; i < dim; i++) {
            p[i] = curCycleD [j];
            if (j < maxj) j++;  else j = 0;
         }
         resetToNextPoint();
         return curPointIndex;
      }

      public String formatState() {
         return super.formatState() + PrintfFormat.NEWLINE +
           "Current cycle: " + curCycleIndex;
      }
   }

}