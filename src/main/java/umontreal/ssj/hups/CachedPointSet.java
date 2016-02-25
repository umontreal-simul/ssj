/*
 * Class:        CachedPointSet
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

/**
 * This container class caches a point set by precomputing and storing its
 * points locally in an array. This can be used to speed up computations when
 * using a small low-dimensional point set more than once.
 *
 *  After the points are stored in the array, this class uses the default
 * methods and the default iterator type provided by the base class
 * @ref PointSet. This is one of the rare cases where direct use of the
 * #getCoordinate method is efficient.
 * @remark **Pierre:** We could also implement an iterator that directly
 * returns `x[i][j]` instead of calling `getCoordinate`, for slightly better
 * efficiency. On the other hand, even better efficiency can be achieved by
 * getting an entire point at a time in an array.
 *
 *  However, it might require too much memory for a large point set.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class CachedPointSet extends PointSet {
   protected PointSet P;        // Original PointSet which is cached here.
   protected double x[][];      // Cached points.
   protected CachedPointSet() {}

   /**
    * Creates a new `PointSet` object that contains an array storing the
    * first `dim` coordinates of the first `n` points of `P`. The original
    * point set `P` itself is not modified.
    *  @param P            point set to be cached
    *  @param n            number of points
    *  @param dim          number of dimensions of the points
    */
   public CachedPointSet (PointSet P, int n, int dim) {
      if (P.getNumPoints() < n)
         throw new IllegalArgumentException(
            "Cannot cache more points than in point set P.");
      if (P.getDimension() < dim)
         throw new IllegalArgumentException(
            "Cannot cache points with more coordinates than the dimension.");
      numPoints = n;
      this.dim = dim;
      this.P = P;
      init ();
   }

   protected void init () {
      PointSetIterator itr = P.iterator();
      x = new double[numPoints][dim];
      for (int i = 0; i < numPoints; i++)
         itr.nextPoint (x[i], dim);
   }

   /**
    * Creates a new `PointSet` object that contains an array storing the
    * points of `P`. The number of points and their dimension are the same
    * as in the original point set. Both must be finite.
    *  @param P            point set to be cached
    */
   public CachedPointSet (PointSet P) {
      numPoints = P.getNumPoints();
      dim = P.getDimension();
      if (numPoints == Integer.MAX_VALUE)
         throw new IllegalArgumentException(
            "Cannot cache infinite number of points");
      if (dim == Integer.MAX_VALUE)
         throw new IllegalArgumentException(
            "Cannot cache infinite dimensional points");
      this.P = P;
      init ();
   }

   /**
    * Add the shift to the contained point set and recaches the points.
    * See the doc of the overridden method
    * {@link umontreal.ssj.hups.PointSet.addRandomShift(int,int,RandomStream)
    * addRandomShift(d1, d2, stream)} in  @ref PointSet.
    */
   public void addRandomShift(int d1, int d2, RandomStream stream) {
        P.addRandomShift(d1, d2, stream);
        init();
   }

   /**
    * Randomizes the underlying point set using `rand` and recaches the
    * points.
    */
   public void randomize (PointSetRandomization rand) {
      P.randomize(rand);
      init();
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

}