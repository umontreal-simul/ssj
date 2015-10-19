/*
 * Class:        ContainerPointSet
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
import umontreal.ssj.rng.RandomStream;

/**
 * This acts as a generic base class for all *container classes* that contain
 * a point set and apply some kind of transformation to the coordinates to
 * define a new point set. One example of such transformation is the
 * *antithetic* map, applied by the container class  @ref AntitheticPointSet,
 * where each output coordinate @f$u_{i,j}@f$ is transformed into
 * @f$1-u_{i,j}@f$. Another example is  @ref RandShiftedPointSet.
 *
 * The class implements a specialized type of iterator for container point
 * sets. This type of iterator contains itself an iterator for the contained
 * point set and uses it to access the points and coordinates internally,
 * instead of maintaining itself indices for the current point and current
 * coordinate.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public abstract class ContainerPointSet extends PointSet {
   protected PointSet P;                 // contained point set

   /**
    * Initializes the container point set which will contain point set
    * `P0`. This method must be called by the constructor of any class
    * inheriting from  @ref ContainerPointSet.
    *  @param P0           contained point set
    */
   protected void init (PointSet p0) {
      P = p0;
//      this.dim = P.getDimension();
//      this.numPoints = P.getNumPoints();
   }

   /**
    * Returns the (untransformed) point set inside this container.
    *  @return the point set inside this container
    */
   public PointSet getOriginalPointSet() {
      return P;
   }

   /**
    * Returns the dimension of the contained point set.
    *  @return the dimension of the contained point set
    */
   public int getDimension() {
      return P.getDimension();
   }

   /**
    * Returns the number of points of the contained point set.
    *  @return the number of points of the contained point set
    */
   public int getNumPoints() {
      return P.getNumPoints();
   }


   public double getCoordinate(int i, int j) {
      return P.getCoordinate (i, j);
   }

   public PointSetIterator iterator(){
      return new ContainerPointSetIterator();
   }

/**
 * Randomizes the contained point set using `rand`.
 *  @param rand          @ref PointSetRandomization to use
 */
public void randomize (PointSetRandomization rand) {
       P.randomize(rand);
   }

   /**
    * Calls `addRandomShift(d1, d2, stream)` of the contained point set.
    *  @param d1           lower dimension of the random shift
    *  @param d2           upper dimension of the random shift
    *  @param stream       the random stream
    */
   public void addRandomShift (int d1, int d2, RandomStream stream) {
      P.addRandomShift (d1, d2, stream);
   }

   /**
    * Calls `addRandomShift(stream)` of the contained point set.
    *  @param stream       the random stream
    */
   public void addRandomShift (RandomStream stream) {
      P.addRandomShift (stream);
   }

   /**
    * Calls `clearRandomShift()` of the contained point set.
    */
   public void clearRandomShift() {
      P.clearRandomShift ();
   }


   public String toString() {
      return "Container point set of: {" + PrintfFormat.NEWLINE
              + P.toString() + PrintfFormat.NEWLINE + "}";
   }


   // ********************************************************
   protected class ContainerPointSetIterator extends DefaultPointSetIterator {

      protected PointSetIterator innerIterator = P.iterator();

      public void setCurCoordIndex (int j) {
         innerIterator.setCurCoordIndex (j);
      }

      public void resetCurCoordIndex() {
         innerIterator.resetCurCoordIndex();
      }

      public boolean hasNextCoordinate() {
         return innerIterator.hasNextCoordinate();
      }

      public double nextCoordinate() {
         return innerIterator.nextCoordinate();
      }

      public void setCurPointIndex (int i) {
         innerIterator.setCurPointIndex(i);
      }

      public void resetCurPointIndex() {
         innerIterator.resetCurPointIndex();
      }

      public int resetToNextPoint() {
         return innerIterator.resetToNextPoint();
      }

      public boolean hasNextPoint() {
        return innerIterator.hasNextPoint();
      }

   }
}