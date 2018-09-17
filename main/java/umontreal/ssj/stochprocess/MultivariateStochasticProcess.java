/*
 * Class:        MultivariateStochasticProcess
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
package umontreal.ssj.stochprocess;

/**
 * This class is a multivariate version of @ref StochasticProcess where the process evolves in
 * the @f$c@f$-dimensional real space. It is an abstract (base) class for a multivariate stochastic
 * process @f$\{\mathbf{X}(t) = (X_1(t),…,X_c(t)), t \geq0 \}@f$, sampled (or observed) at a finite
 * number of time points, @f$0 = t_0 < t_1 < \cdots< t_d@f$. The observation times can be specified
 * by #setObservationTimes. The method #generatePath
 * generates @f$\mathbf{X}(t_1),…,\mathbf{X}(t_d)@f$ and memorizes them in a (one-dimensional)
 * vector, which can be recovered by #getPath. The element @f$cj+i-1@f$ of this vector
 * contains @f$X_i(t_j)@f$, for @f$j=0,…,d@f$ and
 * 
 * @f$i=1,…,c@f$. Alternatively, in some cases, the observations
 * @f$\mathbf{X}(t_j)@f$ can be generated sequentially, one at a time, by invoking
 *                       #resetStartProcess first, and then #nextObservationVector repeatedly.
 *
 *                       <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public abstract class MultivariateStochasticProcess extends StochasticProcess {

	protected double[] x0;         // Default value of the process at time t_0.
	protected int c = 1;                // Dimension of the process.

	/**
	 * Generates, returns, and saves the sample path
	 * 
	 * @f$\{\mathbf{X}(t_0), \mathbf{X}(t_1), …, \mathbf{X}(t_d)\}@f$, which can then be accessed
	 *                       via `getPath`, `getSubpath`, or `getObservation`. The generation method
	 *                       depends on the process type. If `path[]` denotes the returned array,
	 *                       then `path[cj + i-1]` contains @f$X_i(t_j)@f$ for @f$j=0,…,d@f$
	 *                       and @f$i=1,…,c@f$.
	 */
	public abstract double[] generatePath();

	/**
	 * Returns in `subpath` the values of the process at a subset of the observation times,
	 * specified as the times @f$t_j@f$ whose indices
	 * 
	 * @f$j@f$ are in the array `pathIndices`. The size of `pathIndices` should be at least as much
	 *         as that of `subpath`.
	 */
	public void getSubpath(double[] subpath, int[] pathIndices) {
		for (int j = 0; j < pathIndices.length; j++) {
			for (int i = 0; i < c; i++) {
				subpath[c * j + i] = path[c * pathIndices[j] + i];
			}
		}
	}

	/**
	 * Sets the observation times of the process to a copy of `t`, with
	 * 
	 * @f$t_0 = \mathtt{t[0]}@f$ and @f$t_d = \mathtt{t[d]}@f$. The size of `t` must be @f$d+1@f$.
	 */
	public void setObservationTimes(double[] t, int d) {
		if (d <= 0)
			throw new IllegalArgumentException("Number of observation times d <= 0");

		this.d = d;
		observationTimesSet = true;

		/*** Copy of the observation times ***/
		this.t = new double[d + 1];
		System.arraycopy(t, 0, this.t, 0, d + 1);

		/*** Test chronological order ***/
		for (int i = 0; i < d; i++) {
			if (t[i + 1] < t[i])
				throw new IllegalArgumentException(
				        "Observation times t[] are not ordered chronologically");
		}

		/***
		 * Process specific initialization; usually precomputes quantities that depend on the
		 * observation times.
		 ***/
		init();
	}

	/**
	 * Returns @f$\mathbf{X}(t_j)@f$ in the @f$c@f$-dimensional vector `obs`.
	 */
	public void getObservation(int j, double[] obs) {
		for (int i = 0; i < c; i++)
			obs[i] = path[c * j + i];
	}

	/**
	 * Returns @f$X_i(t_j)@f$ from the current sample path.
	 */
	public double getObservation(int j, int i) {
		return path[c * j + i];
	}

	/**
	 * Generates and returns in `obs` the next observation
	 * 
	 * @f$\mathbf{X}(t_j)@f$.
	 */
	public abstract void nextObservationVector(double[] obs);

	/**
	 * Returns the value of the last generated observation
	 * 
	 * @f$\mathbf{X}(t_j)@f$.
	 */
	public void getCurrentObservation(double[] obs) {
		for (int i = 0; i < c; i++)
			obs[i] = path[c * observationIndex + i];
	}

	/**
	 * Returns in `x0` the initial value @f$\mathbf{X}(t_0)@f$ for this process.
	 */
	public double[] getX0(double[] x0) {
		for (int i = 0; i < c; i++)
			x0[i] = this.x0[i];
		return x0;
	}

	protected void init() {
		if (observationTimesSet)
			createPath();
		if (path != null) // If observation times are not defined, do nothing.
			for (int i = 0; i < c; i++)
				path[i] = x0[i];
	}

	/*** Called by 'init' to create new path ***/
	protected void createPath() {
		path = new double[c * (d + 1)];
	}

	/**
	 * Returns the dimension of @f$\mathbf{X}@f$.
	 */
	public int getDimension() {
		return c;
	}

}