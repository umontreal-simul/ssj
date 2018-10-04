/*
 * Class:        LatinHypercube
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
import umontreal.ssj.rng.*;

/**
 * Implements Latin Hypercube Sampling (LHS) with @f$n@f$ points in the
 * @f$s@f$-dimensional unit hypercube. Each axis of the hypercube is divided
 * into @f$n@f$ intervals of length @f$1/n@f$. The randomized point set has
 * the property that for each coordinate @f$j@f$, there is exactly one point
 * with coordinate @f$j@f$ in each interval on length @f$1/n@f$.
 *
 * To generate the points, for each coordinate @f$j@f$ we generate a random
 * permutation @f$\pi_j@f$ of the @f$n@f$ intervals. Then, for point
 * @f$i@f$, we generate coordinate @f$j@f$ randomly uniformly in the interval
 * number @f$\pi_j(i)@f$. The points are enumerated by order of their first
 * coordinate.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class LatinHypercube extends CachedPointSet {

    protected double delta;   // Size of divisions in each dimension, = 1/n

    /**
     * Constructs the structure for a Latin hypercube with `n` points in
     * <tt>dim</tt> dimensions.
     *  @param n            number of points
     *  @param dim          dimension of the points set
     */
    public LatinHypercube (int n, int dim) {
        this.dim = dim;
        delta = 1.0 / n;
        numPoints=n;
   	    x = new double[numPoints][dim];  
    }

	/**
	 * This randomization generates a random LHS point set. The points are sorted by their first
	 * coordinate. The LHS points are defined only after this method has been called.
	 * 
	 * @param stream
	 *            Random stream to generate the permutations and the random points in the intervals
	 */
	public void randomize(RandomStream stream) {
		int[] permutation = new int[numPoints];
		// RandomPermutation.init (permutation, numPoints);
		// Generate one random number uniformly in each interval for first coord.
		for (int i = 0; i < numPoints; i++) {
			permutation[i] = i;
			x[i][0] = (i + stream.nextDouble()) * delta;
		}
		for (int j = 1; j < dim; j++) {
			RandomPermutation.shuffle(permutation, stream);
			// Generate one random number uniformly in each interval for coord. j.
			for (int i = 0; i < numPoints; i++)
				x[permutation[i]][j] = (i + stream.nextDouble()) * delta;
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
	 * Randomizes the points using LHS, regardless of what `rand` is. Equivalent to `randomize
	 * (rand.getStream)`.
	 */
	public void randomize(PointSetRandomization rand) {
		randomize(rand.getStream());
	}

   public String toString() {
      return "LatinHypercube: LHS over the unit cube in " + dim + "dimensions.";
   }
}