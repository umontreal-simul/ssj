/*
 * Class:        IndependentPointsCached
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

import umontreal.ssj.rng.RandomStream;

/**
 * Similar to IndependentPoints, but the points are all generated and stored (cached) when the point
 * set is randomized. The points are independent and uniformly distributed over
 * the @f$s@f$-dimensional unit hypercube. They are enumerated in the order in which they are
 * generated.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class IndependentPointsCached extends CachedPointSet {

	/**
	 * Constructs the structure for `n` points in `dim` dimensions.
	 * 
	 * @param n
	 *            number of points
	 * @param dim
	 *            dimension of the points set
	 */
	public IndependentPointsCached(int n, int dim) {
		this.dim = dim;
		numPoints = n;
		x = new double[numPoints][dim];
	}

	/**
	 * This randomization generates and stores independent random points.
	 * 
	 * @param stream
	 *            Random stream used to generate the random points
	 */
	public void randomize(RandomStream stream) {
		for (int j = 0; j < dim; j++) {
			for (int i = 0; i < numPoints; i++)
				x[i][j] = stream.nextDouble();
		}
	}

	/**
	 * Random shifts and partial randomizations are irrelevant here, so this method is redefined to
	 * be equivalent to `randomize (stream)`. The parameters `fromDim` and `toDim` are *not used*.
	 */
	public void addRandomShift(int fromDim, int toDim, RandomStream stream) {
		randomize(stream);
	}

	/**
	 * Generates a new set of @f$n@f$ independent points, regardless of what `rand` is. 
	 * Equivalent to `randomize(rand.getStream())`.
	 */
	public void randomize(PointSetRandomization rand) {
		randomize(rand.getStream());
	}

	public String toString() {
		return "IndependentPointsCached: independent points in " + dim + "dimensions.";
	}
}