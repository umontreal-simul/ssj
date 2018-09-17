/*
 * Class:        GeometricLevyProcess
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
package umontreal.ssj.stochprocess;

import umontreal.ssj.rng.*;

/**
 * Abstract class used as a parent class for the exponentiation of a Lévy process @f$X(t)@f$: @f[
 * S(t) = S(0) \exp\left(X(t) + (r - \omega_{RN}) t\right). @f] The interest rate is denoted @f$r@f$
 * and is referred to as `muGeom` in the class below. The risk neutral correction is given
 * by @f$\omega_{RN}@f$ and takes into account risk aversion in the pricing of assets; its value
 * depends on the specific Lévy process that is used. *From Pierre:* Maybe we should remove these
 * parameters and make it more general by just putting a linear trend in @f$X@f$.
 *
 * @ref GeometricNormalInverseGaussianProcess is implemented as a child of this class and so
 *      could @ref GeometricVarianceGammaProcess and
 * @ref GeometricBrownianMotion.
 *
 *      <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public abstract class GeometricLevyProcess extends StochasticProcess {

	protected StochasticProcess levyProcess;
	protected double omegaRiskNeutralCorrection;
	protected double muGeom;  // usually the interest rate
	protected double[] muGeomRNdt; // risk neutral corrected
	protected double[] muGeomRNdT; // risk neutral corrected, from time t0.

	protected void init() {
		super.init();
		if (observationTimesSet) {
			// Telling the g\'eom\'etric process about the observation times
			levyProcess.setObservationTimes(t, d);

			// We need to know in which order the observations are generated
			this.observationIndexFromCounter = levyProcess.getArrayMappingCounterToIndex();

			muGeomRNdt = new double[d];
			for (int i = 0; i < d; i++) {
				muGeomRNdt[i] = (muGeom - omegaRiskNeutralCorrection) * (t[i + 1] - t[i]);
			}
			muGeomRNdT = new double[d + 1];
			for (int i = 0; i <= d; i++) {
				muGeomRNdT[i] = (muGeom - omegaRiskNeutralCorrection) * (t[i] - t[0]);
			}
		}
	}

	/**
	 * Generates a path.
	 */
	public double[] generatePath() {
		double s = x0;
		resetStartProcess();
		double[] arithmPath = levyProcess.generatePath();
		for (int i = 0; i < d; i++) {
			s *= Math.exp(muGeomRNdt[i] + arithmPath[i + 1] - arithmPath[i]);
			path[i + 1] = s;
		}
		observationIndex = d;
		return path;
	}

	/**
	 * Returns the next observation. It will also work on a Lévy process which is sampled using the
	 * bridge order, but it will return the observations in the bridge order. If the underlying Lévy
	 * process is of the PCA type, this method is not usable.
	 */
	public double nextObservation() {
		double levy = levyProcess.nextObservation();
		observationIndex = levyProcess.getCurrentObservationIndex();
		path[observationIndex] = x0 * Math.exp(muGeomRNdT[observationIndex] + levy);
		return path[observationIndex];
	}

	/**
	 * Resets the step counter of the geometric process and the underlying Lévy process to the start
	 * value.
	 */
	public void resetStartProcess() {
		super.init();
		levyProcess.resetStartProcess();
	}

	/**
	 * Sets the observation times on the geometric process and the underlying Lévy process.
	 */
	public void setObservationTimes(double[] time, int d) {
		super.setObservationTimes(time, d);
		levyProcess.setObservationTimes(time, d);
	}

	/**
	 * Returns the risk neutral correction.
	 */
	public double getOmega() {
		return omegaRiskNeutralCorrection;
	}

	/**
	 * Returns the geometric drift parameter, which is usually the interest rate, @f$r@f$.
	 */
	public double getMuGeom() {
		return muGeom;
	}

	/**
	 * Sets the drift parameter (interest rate) of the geometric term.
	 */
	public void setMuGeom(double muGeom) {
		this.muGeom = muGeom;
	}

	/**
	 * Returns the Lévy process.
	 */
	public StochasticProcess getLevyProcess() {
		return levyProcess;
	}

	/**
	 * Changes the value of @f$\omega_{RN}@f$. There should usually be no need to redefine the risk
	 * neutral correction from the value set by the constructor. However it is sometimes not unique,
	 * e.g. in
	 * 
	 * @ref GeometricNormalInverseGaussianProcess @cite fALB04a&thinsp;.
	 */
	public void resetRiskNeutralCorrection(double omegaRN) {
		omegaRiskNeutralCorrection = omegaRN;
		init();
	}

	/**
	 * Returns the stream from the underlying Lévy process. If the underlying Lévy process has
	 * multiple streams, it returns what the `getStream()` method of that process was made to
	 * return.
	 */
	public RandomStream getStream() {
		return levyProcess.getStream();
	}

	/**
	 * Resets the stream in the underlying Lévy process. If the underlying Lévy process has multiple
	 * streams, it sets the streams on this process in the same way as `setStream()` for that
	 * process.
	 */
	public void setStream(RandomStream stream) {
		levyProcess.setStream(stream);
	}

}