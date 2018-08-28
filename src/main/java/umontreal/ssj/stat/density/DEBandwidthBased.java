package umontreal.ssj.stat.density;

import java.util.ArrayList;

import umontreal.ssj.mcqmctools.MonteCarloModelDouble;

/**
 * Implements an abstract class for density estimators that rely on the
 * selection of a bandwidth \f$h\f$ such as, e.g., \ref DEHistogram or \ref
 * DEKernelDensity. Furthermore, it provides basic methods such as to estimate
 * the integrated variance (IV), for instance.
 * 
 * The IV of such estimators usually locally follow the model
 * \f$\textrm{IV}\approx C n^{-\beta}h^{-\delta}\f$ and their integrated square
 * bias (ISB) is of the form \f$\textrm{ISB}\approx B h^{\alpha}\f$, where
 * \f$n\f$ denotes the number of data points used to construct the density.
 * 
 * @author puchhamf
 *
 */
public abstract class DEBandwidthBased extends DensityEstimator {

	protected double h;
	/**<bandwidth */
	protected double[] theHs;
	/**<array containing various bandwidths */

	protected double alpha;
	/**<model parameters for ISB */
	protected double beta;
	/**<model parameters for IV */
	protected double delta;
	/**<model parameters for IV */
	protected double C;
	/**<model parameters for IV */
	protected double B;
	/**<model parameters for ISB */

	protected double a;
	/**<left boundary of the interval over which we estimate. */
	protected double b;

	/**<right boundary of the interval over which we estimate. */

	/**
	 * Gives the bandwidth.
	 * 
	 * @return the bandwidth.
	 */
	public double getH() {
		return h;
	}

	/**
	 * Sets the bandwidth to the value of \f$h\f$.
	 * 
	 * @param h
	 *            the desired bandwidth
	 */
	public void setH(double h) {
		this.h = h;
	}

	/**
	 * Gives the array {@link #theHs} containing various bandwidths.
	 * 
	 * @return array with various bandwidths.
	 */
	public double[] getTheHs() {
		return theHs;
	}

	/**
	 * Sets the optional array containing various bandwidths {@link #theHs} to the
	 * value of \a tH.
	 * 
	 * @param tH
	 *            the desired array of bandwidths
	 */
	public void setTheHs(double[] tH) {
		this.theHs = new double[tH.length];
		this.theHs = tH;

	}

	/**
	 * Gives the current value of {@link #alpha}.
	 * 
	 * @return {@link #alpha}
	 */
	public double getAlpha() {
		return alpha;
	}

	/**
	 * Sets the current value of {@link #alpha} to \a alpha.
	 * 
	 * @param alpha
	 *            the desired value for {@link #alpha}.
	 */
	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}

	/**
	 * Gives the current value of {@link #B}.
	 * 
	 * @return {@link #B}
	 */
	public double getB() {
		return B;
	}

	/**
	 * Sets the current value of {@link #B} to \a B.
	 * 
	 * @param B
	 *            the desired value for {@link #B}
	 */
	public void setB(double B) {
		this.B = B;
	}

	/**
	 * Sets the current value of {@link #beta} to \a beta$.
	 * 
	 * @param beta
	 *            the desired value for {@link #beta}.
	 */
	public void setBeta(double beta) {
		this.beta = beta;
	}

	/**
	 * Gives the current value of {@link #beta}.
	 * 
	 * @return {@link #beta}
	 */
	public double getBeta() {
		return beta;
	}

	/**
	 * Sets the current value of {@link #C} to \a C.
	 * 
	 * @param C
	 *            the desired value for {@link #C}.
	 */
	public void setC(double C) {
		this.C = C;
	}

	/**
	 * Gives the current value of {@link #C}.
	 * 
	 * @return {@link #C}
	 */
	public double getC() {
		return C;
	}

	/**
	 * Gives the current value of {@link #delta}.
	 * 
	 * @return {@link #delta}
	 */
	public double getDelta() {
		return delta;
	}

	/**
	 * Sets the current value of {@link #delta} to \a delta.
	 * 
	 * @param delta
	 *            the desired value for {@link #delta}.
	 */
	public void setDelta(double delta) {
		this.delta = delta;
	}

	/**
	 * Gives the estimated IV based on the (local) model assumption
	 * \f$Cn^{-\beta}h^{-\delta} \f$. Note that this requires the parameters \a C,
	 * \a beta, and \a delta to be set for this estimator.
	 * 
	 * @param n
	 *            the number of observations.
	 * @return the estimated IV.
	 */
	public double computeDensityEstimatedIV(int n) {
		return C * Math.pow(n, -beta) * Math.pow(h, -delta);
	}

	/**
	 * Gives the estimated ISB for when the exact density of the underlying model is
	 * not known, based on its asymptotic value \f$Bh^{alpha}\f$. Note that this
	 * requires the parameter \a B to be set.
	 * 
	 * @return the estimated ISB.
	 */
	public double computeDensityISB() {
		return B * Math.pow(h, alpha);
	}

	/**
	 * Same as {@link #computeDensityISB()} with dummy arguments to overload
	 * {@link #computeDensityISB(umontreal.ssj.mcqmctools.MonteCarloModelDensityKnown, double[][], double[])}
	 * for situations, where the true density is not known.
	 * 
	 * @return the estimated ISB.
	 */
	public double computeDensityISB(MonteCarloModelDouble model, double[][] data, double[] evalPoints) {
		return computeDensityISB();
	}

	/**
	 * Same as {@link #computeDensityISB()} with dummy arguments to overload
	 * {@link #computeDensityISB(umontreal.ssj.mcqmctools.MonteCarloModelDensityKnown, double[][], int)}
	 * for situations, where the true density is not known.
	 * 
	 * @return the estimated ISB.
	 */
	public double computeDensityISB(MonteCarloModelDouble model, double[][] data, int numEvalPoints) {
		return computeDensityISB();
	}

	/**
	 * Computes the estimated MISE, i.e. the sum of
	 * {@link #computeDensityEstimatedIV(int)} and {@link #computeDensityISB()}.
	 * 
	 * @param n
	 *            the number of observations.
	 * @return the estimated MISE
	 */
	public double computeDensityEstimatedMISE(int n) {
		return computeDensityISB() + computeDensityEstimatedIV(n);
	}

	/**
	 * Computes the semi-empirical MISE for situations, where the true density is
	 * not known. I.e., it takes the sum of the empirical IV
	 * {@link #computeDensityIV(double[][], double[])} and the estimated ISB
	 * {@link #computeDensityISB()}.
	 * 
	 * @param data
	 *            the observations to construct the density.
	 * @param evalPoints
	 *            the points used to compute the empirical IV.
	 * @return the semi-empirical MISE
	 */
	public double computeDensityMISE(double[][] data, double[] evalPoints) {
		double iv = computeDensityIV(data, evalPoints);
		return iv + computeDensityISB();
	}

	public double computeDensityMISE(double[][] data, int numEvalPoints) {
		return computeDensityMISE(data, equidistantPoints(numEvalPoints));
	}

	/**
	 * Same as {@link #computeDensityMISE(double[][], double[])} but with a dummy
	 * argument \a model to overload
	 * {@link #computeDensityMISE(umontreal.ssj.mcqmctools.MonteCarloModelDensityKnown, double[][], double[])}
	 * for cases where the true density is not known.
	 * 
	 * @param model
	 * @param data
	 *            the observations to construct the density.
	 * @param evalPoints
	 *            the points used to compute the empirical IV.
	 * @return the semi-empirical MISE
	 */
	public double computeDensityMISE(MonteCarloModelDouble model, double[][] data, double[] evalPoints) {
		return computeDensityMISE(data, evalPoints);
	}

	/**
	 * Same as {@link #computeDensityIV(double[][], int)} but with a dummy argument
	 * \a model to overload
	 * {@link #computeDensityMISE(umontreal.ssj.mcqmctools.MonteCarloModelDensityKnown, double[][], double[])}
	 * for cases where the true density is not known.
	 * 
	 * @param model
	 * @param data
	 *            the observations to construct the density.
	 * @param evalPoints
	 *            the points used to compute the empirical IV.
	 * @return the semi-empirical MISE.
	 */
	public double computeDensityMISE(MonteCarloModelDouble model, double[][] data, int numEvalPoints) {
		return computeDensityMISE(data, numEvalPoints);
	}

	/**
	 * Computes all traits for the bandwidth based density estimator that are
	 * specified in the list \a traitList. The possible choices of traits are "isb",
	 * "iv", "estimatediv", "mise" and "estimatedmise", with deviations due to
	 * letter-capitalization considered.
	 * 
	 * Note that in order to be able to compute the traits, several parameters of
	 * the density estimator need to be known and set in advance. Conversely, if
	 * only the ISB is to be computed, \a data and \a evalPoints can be empty.
	 * 
	 * @param traitsList
	 *            the traits which shall be computed.
	 * 
	 *            The computed traits are returned in an array in the same order as
	 *            they appear in \a traitsList.
	 * @param data
	 *            matrix containing the observations of \f$m\f$ independent
	 *            repetitions.
	 * @param evalPoints
	 *            the points at which the density estimator is evaluated to compute
	 *            the traits.
	 * @return the values of the specified traits.
	 */
	public double[] computeDensityTraits(ArrayList<String> traitsList, double[][] data, double[] evalPoints) {
		int t = traitsList.size();
		double[] traitsVals = new double[t];
		String traitName;
		//
		double iv = 0.0;
		double isb = 0.0;

		// check if "mise" is in the list. Not done with "contains", because "mise"
		// might be written in uppercase, etc.
		boolean containsMise = false;
		for (int s = 0; s < t; s++) {
			traitName = traitsList.get(s).toLowerCase();
			if (traitName == "mise")
				containsMise = true;
		}

		// if "mise" is in the list compute iv and isb automatically
		if (containsMise) {
			iv = computeDensityIV(data, evalPoints);
			isb = computeDensityISB();
		}

		// go through all traits of the list and add value to output array. If "mise" is
		// not included in the list, also compute the traits.
		for (int s = 0; s < t; s++) {
			traitName = traitsList.get(s).toLowerCase();
			switch (traitName) {
			case "iv":
				if (!containsMise) {
					iv = computeDensityIV(data, evalPoints);
				}
				traitsVals[s] = iv;
				break;
			case "estimatediv":
				traitsVals[s] = computeDensityEstimatedIV(data[0].length);
				break;
			case "isb":
				if (!containsMise)
					isb = computeDensityISB();
				traitsVals[s] = isb;
				break;
			case "mise":
				traitsVals[s] = iv + isb;
			case "estimatedmise":
				traitsVals[s] = computeDensityEstimatedMISE(data[0].length);
			default:
				System.out.println("The trait " + traitName
						+ " cannot be interpreted. Supported traits are 'isb', 'iv', and 'mise'.");
			}
		}

		return traitsVals;
	}

	/**
	 * Same as {@link #computeDensityTraits(ArrayList, double[][], double[])} but
	 * with \a numEvalPoints equidistant evaluation points.
	 * 
	 * @param traitsList
	 *            the traits which shall be computed.
	 * 
	 *            The computed traits are returned in an array in the same order as
	 *            they appear in \a traitsList.
	 * @param data
	 *            matrix containing the observations of \f$m\f$ independent
	 *            repetitions.
	 * @param numEvalPoints
	 *            the number of equidistant points at which the density estimator is
	 *            evaluated to compute the traits.
	 * @return the values of the specified traits.
	 */
	public double[] computeDensityTraits(ArrayList<String> traitsList, double[][] data, int numEvalPoints) {
		return computeDensityTraits(traitsList, data, equidistantPoints(numEvalPoints));
	}

	/**
	 * Same as {@link #computeDensityTraits(ArrayList, double[][], double[])} but
	 * with the dummy argument \a model to overload
	 * {@link #computeDensityTraits(ArrayList, MonteCarloModelDouble, double[][], double[])}
	 * for cases where the true density is not known.
	 * 
	 * @param traitsList
	 *            the traits which shall be computed.
	 * 
	 *            The computed traits are returned in an array in the same order as
	 *            they appear in \a traitsList.
	 * @param model
	 * @param data
	 *            matrix containing the observations of \f$m\f$ independent
	 *            repetitions.
	 * @param evalPoints
	 *            the points at which the density estimator is evaluated to compute
	 *            the traits.
	 * @return the values of the specified traits.
	 */

	public double[] computeDensityTraits(ArrayList<String> traitsList, MonteCarloModelDouble model, double[][] data,
			double[] evalPoints) {
		return computeDensityTraits(traitsList, data, evalPoints);
	}

	/**
	 * Same as {@link #computeDensityTraits(ArrayList, double[][], int)} but with
	 * the dummy argument \a model to overload
	 * {@link #computeDensityTraits(ArrayList, MonteCarloModelDouble, double[][], int)}
	 * for cases where the true density is not known.
	 * 
	 * @param traitsList
	 *            the traits which shall be computed.
	 * 
	 *            The computed traits are returned in an array in the same order as
	 *            they appear in \a traitsList.
	 * @param model
	 * @param data
	 *            matrix containing the observations of \f$m\f$ independent
	 *            repetitions.
	 * @param numEvalPoints
	 *            the number of equidistant points at which the density estimator is
	 *            evaluated to compute the traits.
	 * @return the values of the specified traits.
	 */
	public double[] computeDensityTraits(ArrayList<String> traitsList, MonteCarloModelDouble model, double[][] data,
			int numEvalPoints) {
		return computeDensityTraits(traitsList, data, numEvalPoints);
	}

	/**
	 * {@inheritDoc}
	 */
	public abstract void constructDensity(double[] data);

	/**
	 * {@inheritDoc}
	 */
	public abstract double evalDensity(double x);

	/**
	 * {@inheritDoc}
	 */
	public abstract String toString();

}
