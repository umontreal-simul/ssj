/*
 * Class:        StochasticProcess
 * Description:  Base class for all stochastic processes
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
package umontreal.ssj.stochprocess;

import umontreal.ssj.rng.RandomStream;

/**
 * Abstract base class for a stochastic process @f$\{X(t) : t \geq 0 \}@f$ sampled (or observed) at
 * a finite number of time points, @f$0 = t_0 < t_1 < \cdots< t_d@f$. The observation times are
 * usually all specified before generating a sample path. This can be done via
 * `setObservationTimes`. The method `generatePath` generates @f$X(t_1),\dots,X(t_d)@f$ and
 * memorizes them in a vector, which can be recovered by `getPath`.
 *
 * Alternatively, for some types of processes, the observations @f$X(t_j)@f$ can be generated
 * sequentially, one at a time, by invoking `resetStartProcess` first, and then `nextObservation`
 * repeatedly. For some types of processes, the observation times can be specified one by one as
 * well, when generating the path. This may be convenient or even necessary if the observation times
 * are random, for example.
 *
 * <strong>WARNING:</strong> After having called the constructor for one of the subclass, it is
 * important to set the observation times of the process, usually by calling `setObservationTimes`.
 *
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */

public abstract class StochasticProcess {

	protected boolean observationTimesSet = false;  // Used in some subclasses
	// to make sure 'setObservationTimes' has been invoked before calling 'init'.

	protected double x0 = 0.0;       // Default Initial Value of the process
	protected int d = -1;            // Number of observation times
	protected int observationIndex = 0; // Index of last generated obs. time
	protected int observationCounter = 0; // Counts how many observations have
	// been generated so far. Useful when they are not generated in
	// chronological order.
	protected double[] t;            // Observation times
	protected double[] path;         // Observations of the process (skeleton).
	// protected RandomStream stream; // Random stream used to generate the process
	protected int[] observationIndexFromCounter; // Links counter# to index#

	/**
	 * Sets the observation times of the process to a copy of `T`, with
	 * 
	 * @f$t_0 =@f$ `T[0]` and @f$t_d =@f$ `T[d]`. The size of `T` must be @f$d+1@f$.
	 */
	public void setObservationTimes(double[] T, int d) {
		if (d <= 0)
			throw new IllegalArgumentException("Number of observation times d <= 0");
		this.d = d;
		observationTimesSet = true;

		// Copy of the observation times
		this.t = new double[d + 1];
		System.arraycopy(T, 0, this.t, 0, d + 1);

		// Test chronological order
		for (int i = 0; i < d; i++) {
			if (T[i + 1] < T[i])
				throw new IllegalArgumentException("Observation times T[] are not time-ordered");
		}

		// Construction of 'path' object
		// We do not do it in 'init()' because we do not want to change the
		// path object if the user only calls 'setParams'
		path = new double[d + 1];

		// Process-specific initialization; usually precomputes quantities
		// that depend on the observation times.
		init();
	}

	/**
	 * Sets equidistant observation times at @f$t_j = j\delta@f$, for
	 * 
	 * @f$j=0,\dots,d@f$, and `delta` = @f$\delta@f$.
	 */
	public void setObservationTimes(double delta, int d) {
		t = new double[d + 1];
		for (int i = 0; i <= d; i++)
			t[i] = i * delta;
		setObservationTimes(t, d);
	}

	/**
	 * Returns a reference to the array that contains the observation times
	 * 
	 * @f$(t_0,…,t_d)@f$. *Warning*: This method should only be used to *read* the observation
	 *                    times. Changing the values in the array directly may have unexpected
	 *                    consequences. The method `setObservationTimes` should be used to modify
	 *                    the observation times.
	 */
	public double[] getObservationTimes() {
		return t;
	}

	/**
	 * Returns the number @f$d@f$ of observation times, *excluding* the time @f$t_0@f$.
	 */
	public int getNumObservationTimes() {
		return d;
	}

	/**
	 * Generates, returns, and saves the sample path @f$\{X(t_0), X(t_1), \dots, X(t_d)\}@f$. It can
	 * then be accessed via `getPath`, `getSubpath`, or `getObservation`. The generation method
	 * depends on the process type.
	 */
	public abstract double[] generatePath();

	/**
	 * Same as `generatePath()`, but first resets the stream to `stream`.
	 */
	public double[] generatePath(RandomStream stream) {
		setStream(stream);
		return generatePath();
	}

	/**
	 * Returns a *reference* to the last generated sample path @f$\{X(t_0), ... , X(t_d)\}@f$.
	 * *Warning*: The returned array and its size should not be modified, because this is the one
	 * that memorizes the observations (not a copy of it). To obtain a copy, use `getSubpath`
	 * instead.
	 */
	public double[] getPath() {
		return path;
	}

	/**
	 * Returns in `subpath` the values of the process at a subset of the observation times,
	 * specified as the times @f$t_j@f$ whose indices
	 * 
	 * @f$j@f$ are in the array `pathIndices`. The size of `pathIndices` should be at least as much
	 *         as that of `subpath`.
	 */
	public void getSubpath(double[] subpath, int[] pathIndices) {
		for (int j = 0; j < subpath.length; j++) {
			subpath[j] = path[pathIndices[j]];
		}
	}

	/**
	 * Returns @f$X(t_j)@f$ from the current sample path. *Warning*: If the observation @f$X(t_j)@f$
	 * for the current path has not yet been generated, then the value returned is unpredictable.
	 */
	public double getObservation(int j) {
		return path[j];
	}

	/**
	 * Resets the observation counter to its initial value @f$j=0@f$, so that the current
	 * observation @f$X(t_j)@f$ becomes @f$X(t_0)@f$. This method should be invoked before
	 * generating observations sequentially one by one via #nextObservation, for a new sample path.
	 */
	public void resetStartProcess() {
		observationIndex = 0;
		observationCounter = 0;
	}

	/**
	 * Returns `true` if @f$j<d@f$, where @f$j@f$ is the number of observations of the current
	 * sample path generated since the last call to #resetStartProcess. Otherwise returns `false`.
	 */
	public boolean hasNextObservation() {
		if (observationCounter < d)
			return true;
		else
			return false;
	}

	/**
	 * Generates and returns the next observation @f$X(t_j)@f$ of the stochastic process. The
	 * processes are usually sampled *sequentially*, i.e. if the last observation generated was for
	 * time
	 * 
	 * @f$t_{j-1}@f$, the next observation returned will be for time
	 * @f$t_j@f$. In some cases, subclasses extending this abstract class may use non-sequential
	 *            sampling algorithms (such as bridge sampling). The order of generation of
	 *            the @f$t_j@f$’s is then specified by the subclass. All the processes generated
	 *            using principal components analysis (PCA) do not have this method.
	 */
	public double nextObservation() {
		throw new UnsupportedOperationException("Method not defined in this class");
	}

	/**
	 * Returns the value of the index @f$j@f$ corresponding to the time
	 * 
	 * @f$t_j@f$ of the last generated observation.
	 */
	public int getCurrentObservationIndex() {
		return observationIndex;
	}

	/**
	 * Returns the value of the last generated observation @f$X(t_j)@f$.
	 */
	public double getCurrentObservation() {
		return path[observationIndex];
	}

	/**
	 * Returns the initial value @f$X(t_0)@f$ for this process.
	 */
	public double getX0() {
		return x0;
	}

	/**
	 * Sets the initial value @f$X(t_0)@f$ for this process to `s0`, and reinitializes.
	 */
	public void setX0(double s0) {
		x0 = s0;
		init();
	}

	/**
	 * Resets the random stream of the underlying generator to `stream`.
	 */
	public abstract void setStream(RandomStream stream);

	/**
	 * Returns the random stream of the underlying generator.
	 */
	public abstract RandomStream getStream();

	/***
	 * Called by 'setObservationTimes' to initialize arrays and precompute constants to speed up
	 * execution. See overriding method 'init' in subclasses for details
	 ***/
	protected void init() {
		if (observationTimesSet) // If observation times are not defined, do nothing.
			path[0] = x0;
		// We do this here because the s0 parameter may have changed through
		// a call to the 'setParams' method.
	}

	/**
	 * Returns a reference to an array that maps an integer @f$k@f$ to @f$i_k@f$, the index of the
	 * observation @f$S(t_{i_k})@f$ corresponding to the
	 * @f$k@f$-th observation to be generated for a sample path of this process. If this process is
	 *            sampled sequentially, then this map is trivial (i.e.
	 * @f$i_k = k@f$). But it can be useful in a more general setting where the process is not
	 *        sampled sequentially (for example, by a Brownian or gamma bridge) and one wants to
	 *        know which observations of the current sample path were previously generated or will
	 *        be generated next.
	 */
	public int[] getArrayMappingCounterToIndex() {
		return observationIndexFromCounter;
	}

}