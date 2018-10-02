/*
 * Class:        StratifiedUnitCube
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
 * boxes of same size and orientation. In @f$s@f$ dimensions, it divides the
 * @f$j@f$th axis in @f$k_j \ge1@f$ equal parts of length @f$d_j = 1/k_j@f$,
 * for @f$j=0,…,s-1@f$. This partitions the unit cube @f$[0,1)^s@f$ into @f$n
 * = k_0\cdots k_{s-1}@f$ rectangular boxes of the same size. The point set
 * is defined by generating one random point uniformly in each of those
 * @f$n@f$ boxes.
 *
 * We number the boxes (or strata) as follows: The box determined by the
 * interval @f$[c_j d_j, (c_j+1)d_j)@f$ for @f$j=0,…,s-1@f$ has number @f$i =
 * c_0 + c_1 k_0 + \cdots+ c_{s-1} k_0 \cdots k_{s-2}@f$.  By default, the
 * points are enumerated by order of their box number. One can also create
 * iterators that enumerate the points by order of any given coordinate.
 *
 * Although this class extends `CachedPointSet`, here there is no underlying
 * point set `p` that is cached.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class StratifiedUnitCube extends CachedPointSet {

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
    public StratifiedUnitCube (int[] k, int dim) {
        this.dim = dim;
        delta = new double[dim];
        numDiv = new int[dim];
        numPoints=1;
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
     * <tt>k</tt>. This gives @f$n = k^s@f$ strata (and points), with
     * @f$s = {\mathtt{dim}}@f$.
     */
    public StratifiedUnitCube (int k, int dim) {
        this.dim = dim;
        delta = new double[dim];
        numDiv = new int[dim];
        numPoints=1;
        for (int j=0; j<dim; j++) {
            numDiv[j] = k;
            numPoints *= k;
            delta[j] = 1.0/(double)k;
        }
   	  x = new double[numPoints][dim];  
    }

	/**
	 * This randomization generates one point randomly in its corresponding box, for each of
	 * the @f$n@f$ boxes. The stratified points are defined only after this method has been called.
	 * 
	 * @param stream
	 *            Random stream to generate the
	 * @f$n\times{\mathtt{dim}}@f$ uniforms required to randomize the points
	 */
	public void randomize(RandomStream stream) {
		int[] current = new int[dim];  // current[j] = current division for
		                               // dim j when we enumerate the points
		for (int j = 0; j < dim; j++)
			current[j] = 0;
		for (int i = 0; i < numPoints; i++) {
			// Generate random point in current box; this is point i.
			for (int j = 0; j < dim; j++)
				x[i][j] = (current[j] + stream.nextDouble()) * delta[j];
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
	 * Randomizes the points using stratification, regardless of what `rand` is. Equivalent to
	 * `randomize (rand.getStream)`.
	 */
	public void randomize(PointSetRandomization rand) {
		randomize(rand.getStream());
	}

	public String toString() {
		return "StratifiedUnitCube: stratified point set over the unit cube in " + dim
		        + "dimensions.";
	}
}