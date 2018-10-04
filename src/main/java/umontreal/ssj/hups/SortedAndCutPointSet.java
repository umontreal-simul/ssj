/*
 * Class:        SortedAndCutPointSet
 * Description:  
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Adam L'Archevêque Gaudet  and Pierre L'Ecuyer
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
 import umontreal.ssj.util.sort.MultiDimSort;
 import umontreal.ssj.rng.RandomStream;
 import umontreal.ssj.util.PrintfFormat;
 // import java.util.Comparator;
 // import java.util.Arrays;

/**
 * This class is useful for the Array-RQMC method, in the situation where
 * the Markov chain has a multidimensional state, the RQMC points are sorted
 * once for all, based on their first @f$\ell@f$ coordinates, then these
 * coordinates are removed and only the remaining coordinates are kept and
 * randomized at each step.
 *
 * It builds a sorted point set based on a given original point set `p`. The
 * original points are first sorted based on their first @f$\ell@f$
 * coordinates, via a  @ref umontreal.ssj.util.MultiDimSort in @f$\ell@f$
 * dimensions, and this is used to build an index that stores the
 * corresponding permutation, once for all. Then those @f$\ell@f$
 * coordinates are removed from the points and never seen again. After that,
 * when the points are enumerated, only the remaining coordinates are used
 * and the points are enumerated by following the index. When they are
 * randomized, the randomization is applied to the original point set, but
 * only the remaining coordinates are randomized. If the original point set
 * has dimension @f$d_0@f$, the retained point set will have dimension
 * @f$d_0 - \ell@f$.
 *
 * In this implementation, the retained coordinates of the points are cached,
 * and they are re-cached after reach randomization.
 * @remark **Pierre:** This could be inefficient (and unnecessary) for very
 * large point sets. One could implement an iterator that uses directly the
 * randomized points in `p` without caching them.
 *
 *  To perform the sort, we first cache only the sorting coordinates of the
 * points, with an extra coordinate that holds the original number of each
 * point, and we apply the sort to those cached points. This extra coordinate
 * is used to recover the permutation made by the sort and to produce the
 * index.
 * 
 * The `CachedPointSet` class also offers methods for
 * multi-dimensional sort and for skipping coordinates.
 *
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class SortedAndCutPointSet extends CachedPointSet {
   protected int[] index;        // Index that represents the permutation of the sort.
   protected MultiDimSort sort;  // Seems to be needed only by the constructor.
   protected int numSortCoord;   // Number of coordinates used for the sort.

   /**
    * Takes the first `dim` coordinates of the first `n` points of `p` and
    * creates a `SortedAndCutPointSet` from these points by sorting them
    * according to `sort`. If `sort` is a
    * @ref umontreal.ssj.util.MultiDimSort of dimension @f$d_s@f$, then
    * the ordering of the points is determined from this `sort` acting on
    * the first @f$d_s@f$ coordinates of the original subset. The
    * dimension of the resulting `SortedAndCutPointSet` is `toDim -
    * fromDim`. Typically, in Array-RQMC applications, `fromDim` will be
    * the dimension used by the sort.
    *  @param p            point set to be cached
    *  @param n            number of points from `p` to use
    *  @param dim          dimension to use for the points of `p`, before
    *                      sorting.
    *  @param sort          @ref umontreal.ssj.util.MultiDimSort to use
    *                       for sorting
    */
   public SortedAndCutPointSet (PointSet p, MultiDimSort sort) {
      super (p, 0, p.getNumPoints(), 0, sort.dimension() + 1);
      numPoints = p.getNumPoints();
      numSortCoord = sort.dimension();
      this.sort = sort;
      this.P = p;
      makeIndex();
      dim = p.getDimension() - numSortCoord;
      x = new double[numPoints][dim]; 
      fillCacheByIndex (numSortCoord, dim);
   }


   /** Sorts the cached points in x[][] according to the first cached coordinates
      (the dimension of the sort) and constructs an index that retains the permutation.
			The coordinate after those used by the sort is used to identify the original 
			point number and recover the permutation.  
			After this, one must call fillCacheByIndex to restore the cached points in proper order in x[][]. 
   */
   protected int[] makeIndex() {
      for (int i = 0; i < numPoints; ++i)
         x[i][numSortCoord] = i;  // Adds extra coordinate that saves the point number.
      sort.sort(x, 0, numPoints); 
      index = new int[numPoints];
      for (int i = 0; i < numPoints; ++i)
         index[(int)x[i][numSortCoord]] = i;
      return index;
   }


   /** Called by the constructors and also by \texttt{randomize}. 
       Cache the points in x, skipping the first fromDim coordinates 
			 and taking the next dim ones.
	     The points are enumerated by order of the index.
   */	
   protected void fillCacheByIndex (int fromDim, int dim) {
      PointSetIterator itr = P.iterator();
      for (int i = 0; i < numPoints; ++i){
         itr.nextPoint (x[index[i]], fromDim, dim);
      }
   }

   /**
    * Returns the number of coordinates of each point, which is the dimension of
    * the original point set minus the dimension of the sort.
    */
   public int getDimension() {
      return dim;
   }

   /**
    * Returns the sort used.
    */
   public MultiDimSort getSort() {
      return sort;
   }

   /**
    * Constructs and returns a point set iterator which gets the values
    * directly from the array.
    */
   public PointSetIterator iterator() {
      return new SortedAndCutPointSetIterator();
   }

   /**
    * Add the shift to the contained point set.
    */
   public void addRandomShift(int d1, int d2, RandomStream stream) {
        P.addRandomShift(d1 + numSortCoord, d2 + numSortCoord, stream);
				fillCacheByIndex (d1 + numSortCoord, d2 + numSortCoord);
   }

   /**
    * Randomizes the contained point (all coordinates) set using `rand`.
    *  @param rand          @ref PointSetRandomization to use
    */
   public void randomize (PointSetRandomization rand) {
      P.randomize (rand);
			// Actually we could randomize only from coordinate numSortCoord.    ****
      fillCacheByIndex (numSortCoord, dim);
   }

   /**
    * Formats a string that contains the information about this point set.
    *  @return string representation of the point set information
    */
   public String toString() {
     StringBuffer sb = new StringBuffer ("SortAndCutPointSet " +
          PrintfFormat.NEWLINE);
     sb.append (PrintfFormat.NEWLINE + "Cached point set information {"
                + PrintfFormat.NEWLINE);
     sb.append (P.toString());
     sb.append (PrintfFormat.NEWLINE + "}");
     return sb.toString();
   }


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// This class implements a SortedAndCutPointSet iterator.
// which takes the value in x rather than calling the function
// getCoordinate inherited from CachedPointSet.

   protected class SortedAndCutPointSetIterator extends DefaultPointSetIterator {

      public double nextCoordinate() {
         if (getCurPointIndex() >= numPoints || getCurCoordIndex() >= dim)
            outOfBounds();
         return x[curPointIndex][curCoordIndex++];
      }
   }

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// This class implements a SortedAndCutPointSet iterator that takes the points
// directly in P instead of using the cache. 
/*
   protected class SortedAndCutPointSetIteratorNoCache extends DefaultPointSetIterator {

      public double nextCoordinate() {
         if (getCurPointIndex() >= numPoints || getCurCoordIndex() >= dim)
            outOfBounds();
				 //  Not implemented!!   To be done...
      }
   }
*/
}
