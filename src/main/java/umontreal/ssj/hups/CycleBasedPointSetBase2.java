/*
 * Class:        CycleBasedPointSetBase2
 * Description:  
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
import cern.colt.list.*;

/**
 * Similar to  @ref CycleBasedPointSet, except that the successive values in
 * the cycles are stored as integers in the range @f$\{0,\dots,2^k-1\}@f$, where
 * @f$1\le k \le31@f$. The output values @f$u_{i,j}@f$ are obtained by
 * dividing these integer values by @f$2^k@f$. Point sets where the
 * successive coordinates of each point are obtained via linear recurrences
 * modulo 2 (e.g., linear feedback shift registers or Korobov-type polynomial
 * lattice rules) are naturally expressed in this form. Storing the integers
 * @f$2^k u_{i,j}@f$ instead of the @f$u_{i,j}@f$ themselves makes it easier
 * to apply randomizations such as digital random shifts in base 2, which are
 * applied to the bits *before* transforming the value to a real number
 * @f$u_{i,j}@f$. When a random digital shift is performed, it applies a
 * bitwise exclusive-or of all the points with a single random point.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public abstract class CycleBasedPointSetBase2 extends CycleBasedPointSet {

// dim = Integer.MAX_VALUE;     // Dimension is infinite.
   private int[] digitalShift;  // Digital shift, initially zero (null).
                                // Entry j is for dimension j.
   protected int numBits;       // Number of bits in stored values.
   protected double normFactor; // To convert output to (0,1); 1/2^numBits.


   public double getCoordinate (int i, int j) {
      // Find cycle that contains point i, then index in cycle.
      int l = 0;         // Length of next cycle.
      int n = 0;         // Total length of cycles added so far.
      int k;
      for (k = 0;  n <= i;  k++)
         n += l = ((AbstractList) cycles.get (k)).size();
      AbstractList curCycle = (AbstractList) cycles.get (k-1);
      int[] curCycleI = ((IntArrayList) curCycle).elements();
      int coordinate = (i - n + l + j) % curCycle.size();
      int shift = 0;
      if (digitalShift != null) {
         shift = digitalShift[j];
         return (shift ^ curCycleI[coordinate]) * normFactor + EpsilonHalf;
      } else
         return (shift ^ curCycleI[coordinate]) * normFactor;
   }

   public PointSetIterator iterator() {
      return new CycleBasedPointSetBase2Iterator ();
   }

/**
 * Adds a random digital shift in base 2 to all the points of the point set,
 * using stream `stream` to generate the random numbers, for coordinates `d1`
 * to `d2 - 1`. This applies a bitwise exclusive-or of all the points with a
 * single random point.
 */
public void addRandomShift (int d1, int d2, RandomStream stream) {
      if (null == stream)
         throw new IllegalArgumentException (
              PrintfFormat.NEWLINE +
              "   Calling addRandomShift with null stream");
      if (0 == d2)
         d2 = Math.max (1, dim);
      if (digitalShift == null) {
         digitalShift = new int[d2];
         capacityShift = d2;
      } else if (d2 > capacityShift) {
         int d3 = Math.max (4, capacityShift);
         while (d2 > d3)
            d3 *= 2;
         int[] temp = new int[d3];
         capacityShift = d3;
         for (int i = 0; i < dimShift; i++)
            temp[i] = digitalShift[i];
         digitalShift = temp;
      }
      dimShift = d2;
      int maxj;
      if (numBits < 31) {
         maxj = (1 << numBits) - 1;
      } else {
         maxj = 2147483647;
      }
      for (int i = d1; i < d2; i++)
         digitalShift[i] = stream.nextInt (0, maxj);
      shiftStream = stream;

   }

   /**
    * Erases the current digital random shift, if any.
    */
   public void clearRandomShift() {
      super.clearRandomShift();
      digitalShift = null;
   }


   public String formatPoints() {
      StringBuffer sb = new StringBuffer (toString());
      for (int c = 0; c < numCycles; c++) {
         AbstractList curCycle = (AbstractList)cycles.get (c);
         int[] cycle = ((IntArrayList)curCycle).elements();
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

   public class CycleBasedPointSetBase2Iterator
                   extends CycleBasedPointSetIterator {

      protected int[] curCycleI;          // The array for current cycle

      public CycleBasedPointSetBase2Iterator () {
         super ();
         resetCurCycle (0);
      }

      protected void init() { }

      public void resetCurCycle (int index) {
         curCycleIndex = index;
         curCycle = (AbstractList) cycles.get (index);
         curCycleI = ((IntArrayList) curCycle).elements();
      }

      public double nextCoordinate() {
          // First, verify if there are still points....
          if (curPointIndex >= numPoints)
             outOfBounds();
          int x = curCycleI [curCoordInCycle];
          if (digitalShift != null) {
             if (curCoordIndex >= dimShift)   // Extend the shift.
                addRandomShift (dimShift, curCoordIndex + 1, shiftStream);
             x ^= digitalShift[curCoordIndex];
          }
          curCoordIndex++;
          curCoordInCycle++;
          if (curCoordInCycle >= curCycle.size())
             curCoordInCycle = 0;
          if (digitalShift == null)
             return x * normFactor;
          else
             return x * normFactor + EpsilonHalf;
     }

      public void nextCoordinates (double p[], int dim) {
         // First, verify if there are still points....
         if (curPointIndex >= numPoints)
            outOfBounds();
         if (curCoordIndex + dim >= dimShift)
            addRandomShift (dimShift, curCoordIndex + dim + 1, shiftStream);
         // int j = curCoordInCycle;
         int maxj = curCycle.size();
         int x;
         for (int i = 0; i < dim; i++) {
            x = curCycleI [curCoordInCycle++];
            if (curCoordInCycle >= maxj) curCoordInCycle = 0;
            if (digitalShift == null)
               p[i] = x * normFactor;
            else
               p[i] = (digitalShift[curCoordIndex + i] ^ x) * normFactor + EpsilonHalf;
         }
         curCoordIndex += dim;
      }

      public int nextPoint (double p[], int dim) {
         if (getCurPointIndex() >= getNumPoints ())
            outOfBounds();
         curCoordIndex = 0;
         curCoordInCycle = startPointInCycle;
         nextCoordinates (p, dim);
         resetToNextPoint();
         return curPointIndex;
      }
   }
}