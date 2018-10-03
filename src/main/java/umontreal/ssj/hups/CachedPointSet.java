/*
 * Class:        CachedPointSet
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
import java.util.Arrays;
import umontreal.ssj.util.sort.DoubleArrayComparator;
import umontreal.ssj.util.sort.MultiDimSort;
import umontreal.ssj.rng.RandomStream;

/**
 * This container class caches a point set by precomputing and storing its
 * points locally in an array. This can be used to speed up computations when
 * using a small low-dimensional point set more than once. Some types of
 * point sets are also implemented directly as subclasses of this one,
 * without having a contained point set `p`.  Examples include 
 * `LatinHypercube` and `StratifiedUnitCube`.
 *
 * After the points are stored in the array, this class uses the default
 * methods and the default iterator type provided by the base class
 * @ref PointSet. This iterator relies exclusively on the  #getCoordinate
 * method. This is one (simple) special case where direct use of
 * #getCoordinate is efficient.
 * **Pierre:** We could also implement an iterator that directly
 * returns `x[i][j]` instead of calling `getCoordinate`, for slightly better
 * efficiency. On the other hand, even better efficiency can be achieved by
 * getting an entire point at a time in an array.
 *
 * There is a  #sort  method available to sort the cached points via a
 * @ref umontreal.ssj.util.MultiDimSort, and a  #stripCoordinates  method
 * that remove some coordinates of the cached points. Those methods affect
 * only the array where the points are cached; they have no impact on the
 * (underlying) original point set. The  #randomize  method, on the other
 * hand, randomizes the underlying point set `p` and re-caches them.
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
   boolean randomizeParent = true;

   protected CachedPointSet() {}   // Constructs an empty cache for a point set.

   /**
    * Creates a new `PointSet` object that contains an array storing the
    * coordinates `fromDim` (inclusive) to `toDim` (exclusive) of the
    * points numbered `fromPoint` (inclusive) to `toPoint` (exclusive) of
    * `p`. Thus, the cached points have `dim = fromDim - toDim`
    * coordinates and there are @f$n =
    * \mathtt{toPoint}-\mathtt{fromPoint}@f$ of them. The original point
    * set `p` itself is not modified.
    *  @param p            point set to be cached
    *  @param fromPoint    number of skipped points (usually 0)
    *  @param toPoint      1 + number of the last point selected from `p`
    *  @param fromDim      number of skipped coordinates (often 0)
    *  @param toDim        1 + index of last coordinate in the original
    *                      points
    */
   public CachedPointSet (PointSet p, int fromPoint, int toPoint, int fromDim, int toDim) {
      if (p.getNumPoints() < toPoint - fromPoint)
         throw new IllegalArgumentException(
            "Cannot cache more points than in point set p.");
      if (p.getDimension() < toDim - fromDim)
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
      this.P = p;
      x = new double[numPoints][dim];
      fillCache (fromDim, dim);
   }

   /**
    * Caches the points, skipping `fromDim` coordinates and taking the next `dim` ones.
    */
   protected void fillCache (int fromDim, int dim) {
      PointSetIterator itr = P.iterator();
      if (fromPoint > 0) itr.setCurPointIndex (fromPoint);
      for (int i = 0; i < numPoints; i++)
         itr.nextPoint (x[i], fromDim, dim);
   }

   /**
    * Same as  {@link #CachedPointSet(PointSet,int,int,int,int)
    * CachedPointSet(p, 0, n, 0, dim)}.
    */
   public CachedPointSet (PointSet p, int n, int dim) {
      this (p, 0, n, 0, dim);
   }

   /**
    * Creates a new `PointSet` object that contains an array storing the
    * points of `p`. The number of points and their dimension are the same
    * as in the original point set. Both must be finite.
    *  @param p            point set to be cached
    */
   public CachedPointSet (PointSet p) {
      this (p, 0, p.getNumPoints(), 0, p.getDimension());
   }

   /**
    * Constructs and returns a point set iterator which gets the values
    * directly from the array.
    */
   public PointSetIterator iterator() {
      return new CachedPointSetIterator();
   }

   /**
    * If `randomizeParent` is `true`, calls to randomize() will be defered to the
    * parent point set (this is the default); otherwise, the randomize method of
    * the PointSetRandomization instance is invoked with this CachedPointSet
    * instance as its argument.
    */
   public void setRandomizeParent(boolean randomizeParent) {
      this.randomizeParent = randomizeParent;
   }

   /**
    * Generates a random shift, adds it to the contained point set, and re-caches 
    * the points. See the doc of the overridden method
    * {@link umontreal.ssj.hups.PointSet.addRandomShift(int,int,RandomStream)
    * addRandomShift(d1, d2, stream)} in  @ref PointSet. In case there is no
    * underlying point set, this method should be redefined to randomize the
    * cached points.
    */
   public void addRandomShift(int d1, int d2, RandomStream stream) {
        P.addRandomShift(d1, d2, stream);
        fillCache (fromDim, dim);
   }

  /**
   * Randomizes the underlying point set using `rand` and re-caches the points.
   * In case there is no underlying point set, this method must be redefined to
   * randomize the cached points in the appropriate way.
   * If `setRandomizeParent(false)` was called, this method invokes
   * `rand.randomize(this)` instead of on the parent point set.
   */
   public void randomize (PointSetRandomization rand) {
      if (randomizeParent) {
         rand.randomize(P);
         fillCache (fromDim, dim);
      }
      else {
         rand.randomize(this);
      }
   }

    /**
     * Sorts the *cached* points by increasing order of coordinate `j`. This is
     * useful in the ArrayRQMC simulation method, for example. Note that the sort
     * applies only to the cached points, and not to the original contained point set `p`.
     */
    public void sortByCoordinate (int j) {
       Arrays.sort (x, new DoubleArrayComparator (j));
    }

   /**
    * Sorts the cached points (only) with the given
    * @ref umontreal.ssj.util.MultiDimSort sorting algorithm `sort`. This does
    * not affect the underlying point set `p`.
    */
    public <T> void sort (MultiDimSort<T> sort) {
      sort.sort(x);
      // sort.sort (P);   init(); ?   No.
   }

   /**
    * Removes the first `d` coordinates of each cached point. This does not
    * affect the underlying point set `p`. This could be useful in particular
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
   public PointSet getParentPointSet() {
      return P;
   }

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 

   /**
    * This class implements a `CachedPointSet` iterator which takes the value 
    * directly in the array `x` in which the points are cached, rather than 
    * calling `getCoordinate` inherited from `CachedPointSet`.
    */
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
