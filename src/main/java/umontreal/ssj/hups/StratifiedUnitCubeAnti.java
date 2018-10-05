/*
 * Class:        StratifiedUnitCubeAnti
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
 * This class implements a stratification of the unit cube in rectangular
 * boxes of same size and orientation, similar to  @ref StratifiedUnitCube.
 * The difference is that in each rectangular box, there are two antithetic
 * points with respect to the opposite (smallest and largest) corners of the
 * box. The points are numbered as in  @ref StratifiedUnitCube, except that
 * point @f$i@f$ in  @ref StratifiedUnitCube becomes point @f$2i@f$ here, and
 * its antithetic point in the same box is point number @f$2i+1@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class StratifiedUnitCubeAnti extends CachedPointSet {

    /**
     * Number of divisions in each dimension, `numDiv[j]` for coordinate `j`.
     */
	protected int[] numDiv;   // Number of divisions in each dimension, k_j

	 /**
     * Size of divisions in each dimension, `delta[j] = 1/numDiv[j]`.
     */
	protected double[] delta;

    /**
     * Builds a stratified points set in <tt>dim</tt> dimensions, with
     * <tt>k[j]</tt> intervals in dimension <tt>j</tt>.
     *  @param k            number of intervals in each dimension
     *  @param dim          integer, dimension of the points set
     */
    public StratifiedUnitCubeAnti (int[] k, int dim) {
        this.dim = dim;
        delta = new double[dim];
        numDiv = new int[dim];
        numPoints=2;
        for (int j=0; j<dim; j++) {
            numDiv[j] = k[j];
            numPoints *= k[j];
            delta[j] = 1.0/(double)numDiv[j];
        }
   	  x = new double[numPoints][dim];  
    }

    /**
     * Same as <tt>StratifiedUnitCube (int[] k, int dim)</tt> with all
     * coordinates of the vector <tt>k</tt> equal to the integer
     * <tt>k</tt>. This gives @f$k^s@f$ strata (and @f$n = 2k^s@f$
     * points), with @f$s = {\mathtt{dim}}@f$.
     */
    public StratifiedUnitCubeAnti (int k, int dim) {
        this.dim = dim;
        delta = new double[dim];
        numDiv = new int[dim];
        numPoints=2;
        for (int j=0; j<dim; j++) {
            numDiv[j] = k;
            numPoints *= k;
            delta[j] = 1.0/(double)k;
        }
   	  x = new double[numPoints][dim];  
    }

	/**
	 * This randomization generates one point randomly in its corresponding box, and the compute a
	 * locally antithetic point in the same box, for each of the @f$n@f$ boxes. The stratified
	 * points are defined only after this method has been called.
	 * 
	 * @param stream
	 *            Random stream to generate the
	 * @f$n\times{\mathtt{dim}}@f$ uniforms required to randomize the points
	 */
	public void randomize(RandomStream stream) {
		int i;
		double u;
		int[] current = new int[dim];  // current[j] = current division for
		                               // dim j when we enumerate the points
		for (int j = 0; j < dim; j++)
			current[j] = 0;
		int numBoxes = numPoints / 2;
		for (int b = 0; b < numBoxes; b++) {
			// Generate random point i=2b in box b, and compute locally antithetic point i+1.
			i = 2 * b;
			for (int j = 0; j < dim; j++) {
				u = stream.nextDouble();
				x[i][j] = (current[j] + u) * delta[j];
				x[i + 1][j] = (current[j] + 1.0 - u) * delta[j];
			}
			// Find the next box.
			for (int l = 0; l < dim; l++) {
				current[l]++;
				if (current[l] >= numDiv[l])
					current[l] = 0;   // next, we will add the carry to current[l+1].
				else
					l = dim;  // Exit loop.
			}
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
	 * Randomizes the points using the locally antithetic stratification, regardless of what `rand`
	 * is. Equivalent to `randomize (rand.getStream)`.
	 */
	public void randomize(PointSetRandomization rand) {
		randomize(rand.getStream());
	}

	public String toString() {
		return "StratifiedUnitCubeAnti: locally antithetic stratified point set over the unit cube in "
		        + dim + "dimensions.";
	}
}