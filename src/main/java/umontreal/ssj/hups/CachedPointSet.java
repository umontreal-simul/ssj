/*
 * Class:        CachedPointSet
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
import java.util.Arrays;
import umontreal.ssj.util.DoubleArrayComparator;
import umontreal.ssj.util.MultiDimSort;
import umontreal.ssj.rng.RandomStream;

/**
 * This container class caches a point set by precomputing and storing its
 * points locally in an array. This can be used to speed up computations when
 * using a small low-dimensional point set more than once. Some types of
 * point sets are also implemented directly as subclasses of this one,
 * without having a contained point set `P`.
 *
 *  After the points are stored in the array, this class uses the default
 * methods and the default iterator type provided by the base class
 * @ref PointSet. This iterator relies exclusively on the  #getCoordinate
 * method. This is one (simple) special case where direct use of
 * #getCoordinate is efficient.
 * @remark **Pierre:** We could also implement an iterator that directly
 * returns `x[i][j]` instead of calling `getCoordinate`, for slightly better
 * efficiency. On the other hand, even better efficiency can be achieved by
 * getting an entire point at a time in an array.
 *
 *  However, it may require too much memory for a large point set.
 *
 * There is a  #sort(m) ethod available to sort the cached points via a and
 * @ref umontreal.ssj.util.MultiDimSort, and a  #stripCoordinates(m) ethod
 * that remove some coordinates of the cached points. Those methods affect
 * only the array where the points are cached; they have no impact on the
 * (underlying) original point set. The  #randomize(m) ethod, on the other
 * hand, randomizes the underlying point set `P` and re-caches them.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class CachedPointSet extends PointSet {
   protected PointSet P;        // Original PointSet which is cached in the present object.
   protected double x[][];      // Cached points.
   int fromPoint = 0;           // Number of skipped points (usually 0).
   int fromDim = 0;             // Number of skipped coordinates (usually 0).
   // int numPoints;            // Number of retained points, inherited from PointSet.
	 // int dim;                  // Dimension of the *retained* points, inherited from PointSet.

   protected CachedPointSet() {}   // Constructs an empty cache for a point set.

   /**
    * Creates a new `PointSet` object that contains an array storing the
    * coordinates `fromDim` (inclusive) to `toDim` (exclusive) of the
    * points numbered `fromPoint` (inclusive) to `toPoint` (exclusive) of
    * `P`. Thus, the cached points have `dim = fromDim - toDim`
    * coordinates and there are @f$n =
    * \mathtt{toPoint}-\mathtt{fromPoint}@f$ of them. The original point
    * set `P` itself is not modified.
    *  @param P            point set to be cached
    *  @param fromPoint    number of skipped points (usually 0)
    *  @param toPoint      1 + number of the last point selected from `P`
    *  @param fromDim      number of skipped coordinates (often 0)
    *  @param toDim        1 + index of last coordinate in the original
    *                      points
    */
   public CachedPointSet (PointSet P, int fromPoint, int toPoint, int fromDim, int toDim) {
      if (P.getNumPoints() < toPoint - fromPoint)
         throw new IllegalArgumentException(
            "Cannot cache more points than in point set P.");
      if (P.getDimension() < toDim - fromDim)
         throw new IllegalArgumentException(
            "Attempt to cache points using more coordinates than the dimension.");
      if (toPoint == Integer.MAX_VALUE)
         throw new IllegalArgumentException(
            "Cannot cache infinite number of points");
      if (toDim == Integer.MAX_VALUE)
         throw new IllegalArgumentException(
            "Cannot cache infinite dimensional points");
      this.fromPoint = fromPoint;
      numPoints = toPoint - fromPoint;
      this.fromDim = fromDim;
      dim = toDim - fromDim;
      this.P = P;
      x = new double[numPoints][dim];
      fillCache (fromDim, dim);
   }

   // Cache the points in x, skipping fromDim coordinates and taking the next dim ones.
   protected void fillCache (int fromDim, int dim) {
      PointSetIterator itr = P.iterator();
      if (fromPoint > 0) itr.setCurPointIndex (fromPoint);
      for (int i = 0; i < numPoints; i++)
         itr.nextPoint (x[i], fromDim, dim);
   }

   /**
    * Same as  {@link #CachedPointSet(PointSet,int,int,int,int)
    * CachedPointSet(P, 0, n, 0, dim)}.
    */
   public CachedPointSet (PointSet P, int n, int dim) {
      this (P, 0, n, 0, dim);
   }

   /**
    * Creates a new `PointSet` object that contains an array storing the
    * points of `P`. The number of points and their dimension are the same
    * as in the original point set. Both must be finite.
    *  @param P            point set to be cached
    */
   public CachedPointSet (PointSet P) {
      this (P, 0, P.getNumPoints(), 0, P.getDimension());
   }

   /**
    * Constructs and returns a point set iterator which gets the values
    * directly from the array.
    */
   public PointSetIterator iterator() {
      return new CachedPointSetIterator();
   }

/**
 * Add the shift to the contained point set and re-caches the points. See the
 * doc of the overridden method
 * {@link umontreal.ssj.hups.PointSet.addRandomShift(int,int,RandomStream)
 * addRandomShift(d1, d2, stream)} in  @ref PointSet. In case there is no
 * underlying point set, this method must be redefined to randomize the
 * cached points.
 */
public void addRandomShift(int d1, int d2, RandomStream stream) {
        P.addRandomShift(d1, d2, stream);
        fillCache (fromDim, dim);
   }

/**
 * Randomizes the underlying point set using `rand` and re-caches the points.
 * In case there is no underlying point set, this method must be redefined to
 * randomize the cached points.
 */
public void randomize (PointSetRandomization rand) {
      P.randomize(rand);
      fillCache (fromDim, dim);
   }

/**
 * Sorts the *cached* points by increasing order of coordinate `j`. This is
 * useful in the ArrayRQMC simulation method, for example. Note that the sort
 * applies only to the cached points, and not to the original points in `P`.
 */
public void sortByCoordinate (int j) {
   Arrays.sort (x, new DoubleArrayComparator (j));
}

/**
 * Sorts the cached points (only) with the given
 * @ref umontreal.ssj.util.MultiDimSort sorting algorithm `sort`. This does
 * not affect the underlying point set `P`.
 */
public <T> void sort (MultiDimSort<T> sort) {
      sort.sort(x);
      // sort.sort (P);   init(); ?   No.
   }

/**
 * Removes the first `d` coordinates of each cached point. This does not
 * affect the underlying point set `P`. This could be useful in particular
 * for the ArrayRQMC simulation method.
 */
public void stripCoordinates (int d) {
      for (int i = 0; i < numPoints; i++)
         for (int j = 0; j < d; j++)
            x[i][j] = x[i][j+d];
      dim = dim - d;
   }


   public String toString() {
     StringBuffer sb = new StringBuffer ("Cached point set" +
          PrintfFormat.NEWLINE);
     sb.append (super.toString());
     sb.append (PrintfFormat.NEWLINE + "Cached point set information {"
                + PrintfFormat.NEWLINE);
     sb.append (P.toString());
     sb.append (PrintfFormat.NEWLINE + "}");
     return sb.toString();
   }

   public double getCoordinate (int i, int j) {
      return x[i][j];
   }

   public double[][] getArray () {
      return x;
   }

   /**
    * Returns the reference point set that was passed to the constructor.
    */
   public PointSet getRefPointSet() {
      return P;
   }

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// This class implements a CachedPointSet iterator.
// which takes the value in x rather than calling the function
// getCoordinate inherited from CachedPointSet.

   protected class CachedPointSetIterator extends DefaultPointSetIterator {

      public double nextCoordinate() {
         if (getCurPointIndex() >= numPoints || getCurCoordIndex() >= dim)
            outOfBounds();
         return x[curPointIndex][curCoordIndex++];
      }

      public void nextCoordinates (double p[], int d)  {
         if (getCurCoordIndex() + d > getDimension()) outOfBounds();
         for (int j = 0; j < d; j++)
            p[j] = x[curPointIndex][curCoordIndex++];
      }

   }

}
