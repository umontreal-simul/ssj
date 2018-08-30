package umontreal.ssj.stat.density;

import java.util.ArrayList;

import umontreal.ssj.mcqmctools.MonteCarloModelDensityKnown;
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
	public double computeEstimatedIV(int n) {
		return C * Math.pow(n, -beta) * Math.pow(h, -delta);
	}

	/**
	 * Gives the estimated ISB for when the exact density of the underlying model is
	 * not known, based on its asymptotic value \f$Bh^{\alpha}\f$. Note that this
	 * requires the parameter \a B to be set.
	 * 
	 * @return the estimated ISB.
	 */
	public double computeEstimatedISB() {
		return B * Math.pow(h, alpha);
	}

	/**
	 * Same as {@link #computeEstimatedISB()} with dummy arguments to overload
	 * {@link #computeISB(MonteCarloModelDensityKnown, double[][], double[])}
	 * for situations, where the true density is not known.
	 * 
	 * @return the estimated ISB.
	 */
	public double computeISB(MonteCarloModelDouble model, double[][] estDensities, double[] evalPoints) {
		return computeEstimatedISB();
	}


	/**
	 * Computes the estimated MISE, i.e. the sum of
	 * {@link #computeEstimatedIV(int)} and {@link #computeDensityISB()}.
	 * 
	 * @param n
	 *            the number of observations.
	 * @return the estimated MISE
	 */
	public double computeEstimatedMISE(int n) {
		return computeEstimatedISB() + computeEstimatedIV(n);
	}

	/**
	 * Computes the semi-empirical MISE for situations, where the true density is
	 * not known. I.e., it takes the sum of the empirical IV
	 * {@link #computeIV(double[][])} and the estimated ISB
	 * {@link #computeEstimatedISB()}. 
	 * 
	 * The estimate of the empirical IV is based on \f$m\f$ realizations of the density estimator,
	 * which have previously been evaluated at the \f$k\f$ points stored in \a
	 * evalPoints. The matrix \a estDensities has
	 * dimensions \f$m\times k \f$, i.e. each row contains the evaluations of one
	 * density estimator.
	 * 
	 * @param the \f$m\times k\f$matrix containing the results of evaluating
	 *            \f$m\f$ densities at \f$k\f$ evaluation points.
	 *            the observations to construct the density.
	 * @return the semi-empirical MISE
	 */
	public double computeMISE(double[][] estDensities) {
		double iv = computeIV(estDensities);
		return iv + computeEstimatedISB();
	}


	/**
	 * Same as {@link #computeMISE(double[][], double[])} but with a dummy
	 * argument \a model to overload
	 * {@link #computeMISE(MonteCarloModelDensityKnown, double[][], double[])}
	 * for cases where the true density is not known.
	 * 
	 * @param model
	 * @param estDensities
	 *           the \f$m\times k\f$matrix containing the results of evaluating
	 *            \f$m\f$ densities at \f$k\f$ evaluation points.
	 *            the observations to construct the density.
	 * @return an estimate for the MISE
	 */
	public double computeMISE(MonteCarloModelDouble model, double[][] estDensities, double[] evalPoints) {
		return computeMISE(estDensities);
	}

	/**
	 * This method estimates the 
	 * IV, ISB, and MISE based on \f$m\f$ realizations of the density estimator,
	 * which have previously been evaluated at the \f$k\f$ points stored in \a
	 * evalPoints, and returns them in an array in this order. The matrix \a estDensities has
	 * dimensions \f$m\times k \f$, i.e. each row contains the evaluations of one
	 * density estimator.
	 * 
	 * The estimate for the IV is obtained from an estimate of the empirical IV {@link #computeIV(double[][])},
	 * the ISB is estimated by {@link #computeEstimatedISB()}, and the MISE-estimate is obtained by their sum.
	 * 
	 * 
	 * @param estDensities
	 *            the \f$m\times k\f$matrix containing the results of evaluating
	 *            \f$m\f$ densities at \f$k\f$ evaluation points.
	 * @return estimates for the IV, the ISB, and the MISE.
	 */
	public double[] computeIVandISBandMISE(double[][] estDensities) {
		double iv = computeIV(estDensities);
		double isb = computeEstimatedISB();
		double[] res = {iv,isb,iv+isb};
		return res;
	}
	
	/**
	 * Same as {@link #computeIVandISBandMISE(double[][])} but with dummy arguments to overload
	 * {@link #computeIVandISBandMISE(MonteCarloModelDensityKnown, double[][], double[])} when the
	 * true density is not known.
	 * 
	 * @param model
	 * @param estDensities
	 * @param evalPoints
	 * @return
	 */
	public double[] computeIVandISBandMISE(MonteCarloModelDouble model, double[][] estDensities,
			double[] evalPoints) {
		return computeIVandISBandMISE(estDensities);
	}


}
