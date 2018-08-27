/**
 *  @author Pierre L'Ecuyer
 */
package umontreal.ssj.stat.density;

import java.util.ArrayList;

import umontreal.ssj.mcqmctools.MonteCarloModelDensityKnown;

/**
 * This abstract class implements a univariate density estimator over a finite
 * interval \f$[a,b]\f$. The density is constructed from a set of \f$n\f$
 * observations. Some of these observations can lie outside the interval
 * \f$[a,b]\f$.
 * 
 * The density can be evaluated at a single point \f$x\f$ or on a grid of
 * evaluation points.
 * 
 * Furthermore, this class implements basic routines for the computation of the
 * integrated variance (IV) as well as for calculating basic error measures,
 * such as the integrated square bias (ISB) and the mean integrated square error
 * (MISE) in case the true density is known.
 */
public abstract class DensityEstimator {

	double a, b; // The density can be evaluated over [a, b].

	/**
	 * Sets the interval \f$[a,b]\f$ over which we estimate.
	 * 
	 * @param a
	 *            left boundary of the interval.
	 * @param b
	 *            right boundary of the interval.
	 */
	public void setRange(double a, double b) {
		this.a = a;
		this.b = b;
	}

	/**
	 * Constructs a density estimator from the data points in vector #data.
	 * 
	 * @param data
	 *            the data points.
	 */
	public abstract void constructDensity(double[] data);

	/**
	 * Returns the value of the density evaluated at point x.
	 * 
	 * @param x
	 *            the point at which the density is to be evaluated.
	 * @return the value of the estimated density at x.
	 */
	public abstract double evalDensity(double x);

	/**
	 * Returns in array #density the value of the density at the evaluation points
	 * in #evalPoints. These two arrays must have the same size.
	 * 
	 * By default, it calls <tt>evalDensity(double )</tt> for each element of
	 * #evalpoints.
	 * 
	 * @param evalPoints
	 *            the evaluation points
	 * @param density
	 *            values of the density at these points
	 */
	public void evalDensity(double[] evalPoints, double[] density) {
		for (int i = 0; i < evalPoints.length; i++)
			density[i] = evalDensity(evalPoints[i]);
	}

	/**
	 * Gives the left boundary of the interval over which we estimate.
	 * 
	 * @return the left boundary of the interval
	 */
	public double getMin() {
		return a;
	}

	/**
	 * Gives the right boundary of the interval over which we estimate.
	 * 
	 * @return the right boundary of the interval
	 */
	public double getMax() {
		return b;
	}

	/**
	 * Gives a short description of the estimator.
	 * 
	 * @return a short description.
	 */
	public abstract String toString();

	/**
	 * This method estimates the density and computes the empirical IV based on
	 * \f$m\f$ independently generated simulations stored in \a data. This is done
	 * by numerically integrating the variance with integration nodes \a evalPoints
	 * over the interval \f$[a,b]\f$ on which the density estimator is defined.
	 * 
	 * 
	 * @param data
	 *            matrix in which the realizations of the model are stored.
	 * @param evalPoints
	 *            the integration nodes.
	 * @return the empirical IV
	 */
	public double computeDensityIV(double[][] data, double[] evalPoints) {
		int m = data.length;
		int numEvalPoints = evalPoints.length;

		double x, y;
		// TO DO:
		// If the density estimator is a histogram, here we may reset numEvalPoints to
		// the number of bins of the histogram.
		// int numEvalPoints = evalPoints.length();
		double estimDens[] = new double[numEvalPoints]; // Value of the density at those points
		double meanDens[] = new double[numEvalPoints]; // Average value over the rep replicates
		double varDens[] = new double[numEvalPoints]; // Variance at each evaluation point
		// Arrays.fill(meanDens, 0.0);
		// Arrays.fill(varDens, 0.0);
		for (int rep = 0; rep < m; rep++) {
			// Estimate the density for this rep and evaluate it at the evaluation points
			constructDensity(data[rep]);
			evalDensity(evalPoints, estimDens);
			// Update the empirical mean and sum of squares of centered observations at each
			// evaluation point.
			for (int j = 0; j < numEvalPoints; j++) {
				x = estimDens[j];
				y = x - meanDens[j];
				meanDens[j] += y / (double) (rep + 1);
				varDens[j] += y * (x - meanDens[j]);
			}
		}
		double sumVar = 0.0;
		for (int i = 0; i < numEvalPoints; ++i)
			sumVar += varDens[i];
		return sumVar * (b - a) / (double) (numEvalPoints * (m - 1)); // Empirical integrated variance.
	}

	/**
	 * Same as {@link #computeDensityIV(double[][], double[])} but now the
	 * evaluation points are taken as \a numEvalPoints equidistant points over
	 * \f$[a,b]\f$. More precisely, we fix the distance between two points as
	 * \f$\delta = (b-a)/\textrm{numEvalPoints}\f$ and set the first point as \f$a +
	 * \delta/2\f$.
	 * 
	 * @param data
	 *            matrix in which the realizations of the model are stored.
	 * @param numEvalPoints
	 *            the number of evaluation points.
	 * @return the empirical IV.
	 */
	public double computeDensityIV(double[][] data, int numEvalPoints) {

		double evalPoints[] = new double[numEvalPoints]; // Points at which the density will be evaluated
		double delta = (b - a) / (numEvalPoints);
		for (int j = 0; j < numEvalPoints; j++)
			evalPoints[j] = a + delta * (0.5 + j);
		return computeDensityIV(data, evalPoints);
	}
	
	/**
	 * For a \ref MonteCarloModelDensityKnown this method estimates the density and
	 * computes the ISB based on \f$m\f$ independently generated simulations stored
	 * in \a data. This is done by numerically integrating the squared bias with
	 * integration nodes \a evalPoints over the interval on which the density estimator is defined.
	 * 
	 * 
	 * This is particularly useful for testing the density estimator with a toy model.
	 * 
	 * @param model
	 *            the underlying model.
	 * @param data
	 *            matrix in which the realizations of the model are stored.
	 * @param evalPoints
	 *            the integration nodes.
	 * @return the empirical ISB
	 */
	public double computeDensityISB(MonteCarloModelDensityKnown model, double[][] data,
			double[] evalPoints) {

		int numEvalPoints = evalPoints.length;
		int m = data.length;


		double meanDens[] = new double[numEvalPoints];
		double sqBiasDens[] = new double[numEvalPoints]; // squared bias

		double estimDens[] = new double[numEvalPoints]; // Value of the density at evalPoints

		double x, y, z;

		for (int rep = 0; rep < m; rep++) {
			// Estimate the density and evaluate it at evalPoints
			constructDensity(data[rep]);
			evalDensity(evalPoints, estimDens);

			// Update the empirical mean and sum of squares of centered
			// observations at each evaluation point.
			for (int j = 0; j < numEvalPoints; j++) {
				x = estimDens[j];
				y = x - meanDens[j];
				meanDens[j] += y / (double) (rep + 1);
			}

		}

		for (int j = 0; j < numEvalPoints; j++) {
			z = (meanDens[j] - model.density(evalPoints[j]));
			sqBiasDens[j] = z * z;
		}

		double isb = 0.0;
		for (int i = 0; i < numEvalPoints; ++i)
			isb += sqBiasDens[i];

		return isb * (b - a) / ((double) numEvalPoints);

	}

	
	/**
	 * For a \ref MonteCarloModelDensityKnown this method estimates the density and
	 * computes the MISE based on \f$m\f$ independently generated simulations stored
	 * in \a data. This is done by numerically integrating the mean square error with
	 * integration nodes \a evalPoints over the interval on which the density estimator is defined.
	 * 
	 * 
	 * This is particularly useful for testing the density estimator with a toy model.
	 * 
	 * @param model
	 *            the underlying model.
	 * @param data
	 *            matrix in which the realizations of the model are stored.
	 * @param evalPoints
	 *            the integration nodes.
	 * @return the empirical MISE
	 */
	
	public double computeDensityMISE(MonteCarloModelDensityKnown model, double[][] data, 
			double[] evalPoints) {
		
		int m = data.length;
		int numEvalPoints = evalPoints.length;
	
		double estimDens[] = new double[numEvalPoints]; // Value of the density
														// at evalPoints
		double mseDens[] = new double[numEvalPoints]; // mse at evalPoints


		double y;
		for (int rep = 0; rep < m; rep++) {
			// Estimate the density and evaluate it at eval points
			constructDensity(data[rep]);
			evalDensity(evalPoints, estimDens);

			
			for (int j = 0; j < numEvalPoints; j++) {
				y = estimDens[j] - model.density(evalPoints[j]); // model.density(x);
				mseDens[j] += y * y;
			}
		}

		double sumMISE = 0.0;
		for (int i = 0; i < numEvalPoints; ++i)
			sumMISE += mseDens[i];
		double fact = (b - a) / ((double) (numEvalPoints * m));
		return sumMISE * fact;
	}
}
