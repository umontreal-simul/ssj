/*
 * Class:        ContainerPointSet
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
 * This acts as a generic base class for all *container classes* that contain a point set and apply
 * a specific kind of transformation to the coordinates @f$u_{i,j}@f$ when the points are generated
 * by the iterator.
 * 
 * One example of such transformation is the *antithetic* map, applied by the container class @ref
 * AntitheticPointSet, where each output coordinate @f$u_{i,j}@f$ is transformed
 * into @f$1-u_{i,j}@f$. Another example is @ref RandShiftedMod1PointSet, which applies a random
 * shift modulo 1 (the same shift) to all the points. Thus, the transformation can be deterministic
 * in some cases, and in other cases it may implement a randomization.
 * 
 * The `PointSetIterator` implemented in this `ContainerPointSet` class contains itself an iterator
 * for the contained point set and uses it to access the points and coordinates internally, instead
 * of maintaining itself indices for the current point and current coordinate. It gets the
 * coordinates @f$u_{i,j}@f$ from the internal iterator, applies the appropriate transformation, and
 * returns the transformed coordinates.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public abstract class ContainerPointSet extends PointSet {
	protected PointSet P;                 // contained point set

	/**
	 * Initializes this container point set so that it will contain the point set `p0`. This method
	 * must be called by the constructor of any subclass of the present abstract class.
	 * 
	 * @param p0
	 *            contained point set
	 */
	protected void init(PointSet p0) {
		P = p0;
		// this.dim = P.getDimension();
		// this.numPoints = P.getNumPoints();
	}

	/**
	 * Returns the (untransformed) point set contained in this container.
	 * 
	 * @return the point set inside this container
	 */
	public PointSet getOriginalPointSet() {
		return P;
	}

	/**
	 * Returns the dimension of the contained point set.
	 * 
	 * @return the dimension of the contained point set
	 */
	public int getDimension() {
		return P.getDimension();
	}

	/**
	 * Returns the number of points in the contained point set.
	 * 
	 * @return the number of points in the contained point set
	 */
	public int getNumPoints() {
		return P.getNumPoints();
	}

	/**
	 * By default, returns the *untransformed* coordinate @f$u_{i,j}@f$ for the contained point set.
	 * 
	 * @return @f$u_{i,j}@f$ for the contained point set
	 */
	public double getCoordinate(int i, int j) {
		return P.getCoordinate(i, j);
	}

	/**
	 * Return an iterator for this `ContainerPointSet`.
	 * 
	 * @return a new `ContainerPointSetIterator`
	 */
	public PointSetIterator iterator() {
		return new ContainerPointSetIterator();
	}

	/**
	 * Randomizes the *contained* point set using `rand`.
	 * 
	 * @param rand
	 * @ref PointSetRandomization to use
	 */
	public void randomize(PointSetRandomization rand) {
		rand.randomize(P);
	}

	/**
	 * Calls `addRandomShift(d1, d2, stream)` of the contained point set.
	 * 
	 * @param d1
	 *            lower dimension of the random shift
	 * @param d2
	 *            upper dimension of the random shift
	 * @param stream
	 *            the random stream
	 */
	public void addRandomShift(int d1, int d2, RandomStream stream) {
		P.addRandomShift(d1, d2, stream);
	}

	/**
	 * Calls `addRandomShift(stream)` of the contained point set.
	 * 
	 * @param stream
	 *            the random stream
	 */
	public void addRandomShift(RandomStream stream) {
		P.addRandomShift(stream);
	}

	/**
	 * Calls `clearRandomShift()` of the contained point set.
	 */
	public void clearRandomShift() {
		P.clearRandomShift();
	}

	public String toString() {
		return "Container point set of: {" + PrintfFormat.NEWLINE + P.toString()
		        + PrintfFormat.NEWLINE + "}";
	}

	// ********************************************************
	protected class ContainerPointSetIterator extends DefaultPointSetIterator {

		protected PointSetIterator innerIterator = P.iterator();

		public void setCurCoordIndex(int j) {
			innerIterator.setCurCoordIndex(j);
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

		public void setCurPointIndex(int i) {
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