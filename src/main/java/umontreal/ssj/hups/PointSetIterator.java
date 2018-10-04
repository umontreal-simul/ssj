/*
 * Interface:    PointSetIterator
 * Description:  Iterator over point sets
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

/**
 * This is the interface for *iterators* that permit one to go through the points of a #PointSet and
 * the successive coordinates of these points.
 * 
 * Each @ref PointSetIterator is associated with a given point set and maintains a *current point*
 * index @f$i@f$ and a *current coordinate* index @f$j@f$, which are both initialized to zero.
 * Successive coordinates can be accessed one or many at a time by the methods #nextCoordinate and
 * #nextCoordinates, respectively. The current coordinate index @f$j@f$ can be set explicitely by
 * #setCurCoordIndex and #resetCurCoordIndex. This could be used to skip some coordinates for each
 * point, for example. Similar methods are available for resetting and accessing the current point.
 * The method #nextPoint permits one to enumerate the successive points in natural order.
 *
 * This interface also extends the @ref umontreal.ssj.rng.RandomStream interface. This permits one
 * to replace random numbers by the coordinates of (randomized) quasi-Monte Carlo points without
 * changing the code that calls the generators in a simulation program. That is, the same simulation
 * program can be used for both Monte Carlo and quasi-Monte Carlo simulations. The method
 * #nextDouble does exactly the same as #nextCoordinate, it returns the current coordinate of the
 * current point and advances the current coordinate by one. The substreams correspond to the
 * points, so #resetStartSubstream resets the current point coordinate to zero, #resetNextSubstream
 * resets the iterator to the next point, and #resetStartStream resets the iterator to the first
 * point of the point set.
 *
 * There can be several iterators over the same point set. These iterators are independent from each
 * other. Classes that implement this interface must maintain enough information so that each
 * iterator is unaffected by other iterator’s operations. However, the iterator does not need to be
 * independent of the underlying point set. If the point set is modified (e.g., randomized), the
 * iterator may continue to work as usual. All iterators on a given point set will see the same
 * randomized points, because the randomizations are in the point sets, not in the iterators.
 *
 * Point set iterators are implemented as inner classes because this gives a direct access to the
 * private members (or variables) of the class. This is important for efficiency. They are quite
 * similar to the iterators in Java *collections*.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public interface PointSetIterator extends RandomStream {

	/**
	 * Sets the current coordinate index to @f$j@f$, so that the next calls to #nextCoordinate or
	 * #nextCoordinates will return the values @f$u_{i,j}, u_{i,j+1}, \dots @f$, where @f$i@f$ is the
	 * index of the current point. This can be useful to skip certain coordinates for each point,
	 * for example.
	 * 
	 * @param j
	 *            index of the new current coordinate
	 */
	public void setCurCoordIndex(int j);

	/**
	 * Equivalent to {@link #setCurCoordIndex() setCurCoordIndex(0)}.
	 */
	public void resetCurCoordIndex();

	/**
	 * Returns the index @f$j@f$ of the current coordinate. This may be useful, e.g., for testing if
	 * enough coordinates are still available.
	 * 
	 * @return index of the current coordinate
	 */
	public int getCurCoordIndex();

	/**
	 * Returns `true` if the current point has another coordinate. This can be useful for testing if
	 * coordinates are still available.
	 * 
	 * @return `true` if the current point has another coordinate
	 */
	public boolean hasNextCoordinate();

	/**
	 * Returns the current coordinate @f$u_{i,j}@f$ and advances to the next one. If no current
	 * coordinate is available (either because the current point index has reached the number of
	 * points or because the current coordinate index has reached the number of dimensions), it
	 * throws a #NoSuchElementException.
	 * 
	 * @return value of the current coordinate
	 *
	 * @exception NoSuchElementException
	 *                if no such coordinate is available
	 */
	public double nextCoordinate();

	/**
	 * Returns  in `p` the next `d` coordinates of the current point and advances the current
	 * coordinate index by `d`. If the remaining number of coordinates is too small, a
	 * `NoSuchElementException` is thrown, as in #nextCoordinate.
	 * 
	 * @param p
	 *            array to be filled with the coordinates, starting at index 0
	 * @param d
	 *            number of coordinates to get
	 * @exception NoSuchElementException
	 *                if there are not enough remaining coordinates in the current point
	 */
	public void nextCoordinates(double[] p, int d);

	/**
	 * Resets the current point index to @f$i@f$ and the current coordinate index to zero. If `i` is
	 * larger or equal to the number of points, an exception will *not* be raised here, but only
	 * later if we ask for a new coordinate or point.
	 * 
	 * @param i
	 *            new index of the current point
	 */
	public void setCurPointIndex(int i);

	/**
	 * Equivalent to {@link #setCurPointIndex() setCurPointIndex(0)}.
	 */
	public void resetCurPointIndex();

	/**
	 * Increases the current point index by 1 and returns its new value. If there is no more point,
	 * an exception will be raised only if we ask for a new coordinate or point later on.
	 * 
	 * @return index of the new current point
	 */
	public int resetToNextPoint();

	/**
	 * Returns the index @f$i@f$ of the current point. This can be useful, e.g., for caching point
	 * sets.
	 * 
	 * @return index of the current point
	 */
	public int getCurPointIndex();

	/**
	 * Returns `true` if there is a next point. This can be useful for testing if points are still
	 * available.
	 * 
	 * @return `true` if a next point is available from the iterated point set
	 */
	public boolean hasNextPoint();

	/**
	 * Returns in `p` the next `d` coordinates of the *current* point, starting at coordinate
	 * `fromDim` (i.e., after skipping `fromDim` coordinates), then advances to the next point and
	 * returns the index of the *new* current point. Regardless of the current coordinate index, the
	 * point returned starts from coordinate `fromDim`. After obtaining the last point via this
	 * method, the current point index (returned by the method) is equal to the number of points, so
	 * it is no longer a valid point index. An exception will then be raised if we attempt to
	 * generate additional points or coordinates.
	 *
	 * Specialized implementations of this method often allow for increased efficiency, e.g., for
	 * cycle-based point sets where the cycles (but not the points) are stored explicitly or for
	 * digital nets by allowing non-incremental point enumerations via Gray-code counters.
	 * 
	 * @param p
	 *            array to be filled with the coordinates, starting from array index 0
	 * @param fromDim
	 *            number of coordinates to be skipped
	 * @param d
	 *            number of coordinates to return
	 * @return index of the new current point
	 *
	 * @exception NoSuchElementException
	 *                if there are not enough coordinates available in the current point for filling
	 *                `p`
	 */
	public int nextPoint(double[] p, int fromDim, int d);

	/**
	 * Same as {@link #nextPoint(double[],int,int) nextPoint(p, 0, d)}.
	 */
	public int nextPoint(double[] p, int d);

}