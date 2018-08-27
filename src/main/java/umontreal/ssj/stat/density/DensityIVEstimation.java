package umontreal.ssj.stat.density;

import java.util.ArrayList;

import umontreal.ssj.mcqmctools.MonteCarloModelDensityKnown;
import umontreal.ssj.mcqmctools.MonteCarloModelDouble;

/**
 * This class provides methods to compute density estimators and to estimate their integrated variance (IV),
 * based on data from m independent replications, with n data points each.  
 * The data may have been obtained in any way, e.g., by Monte Carlo or by RQMC.
 * The density estimator is evaluated at #numEvalPoints points either equally spaced or specified explicitly
 * over a given interval [a,b], for each of the m replicates. 
 * By computing the empirical average and variance at each evaluation point, one can estimate the density and the IV.
 */

/**
 * @author Pierre L'Ecuyer
 * 
 */
public class DensityIVEstimation {

	/**
	 * Takes data from previous simulations (m independent replications, with n data
	 * points each), a density estimator de, and an array of evaluation points for
	 * this density over the interval [a,b]. Computes and returns the density
	 * estimator and the empirical integrated variance (IV) for this density
	 * estimator, obtained by estimating the variance at the numEvalPoints
	 * evaluation points provided and taking the average multiplied by (b-a).
	 * 
	 * \remark{<b>Florian</b>: \a m and \a numEvalPoints don't need to be passed.
	 * They can be extracted from \a data and \a evalPoints.}
	 */

	public static double computeDensityVariance(int m, double[][] data, DensityEstimator de, double a, double b,
			int numEvalPoints, double[] evalPoints) {
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
			de.constructDensity(data[rep]);
			de.evalDensity(evalPoints, estimDens);
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
		return sumVar * (b - a) / (numEvalPoints * (m - 1)); // Empirical integrated variance.
	}

	/**
	 * Similar to the previous method, except that the numEvalPoints evaluation
	 * points are equidistant. We partition [a,b] in numEvalPoints equal intervals
	 * and put one point in the middle of each interval. With delta = (b-a) /
	 * numEvalPoints, we put the first point at a + delta/2, and then the successive
	 * points are spaced by delta.
	 */
	public static double computeDensityVariance(int m, double[][] data, DensityEstimator de, double a, double b,
			int numEvalPoints) {
		double evalPoints[] = new double[numEvalPoints]; // Points at which the density will be evaluated
		double delta = (b - a) / (numEvalPoints);
		for (int j = 0; j < numEvalPoints; j++)
			evalPoints[j] = a + delta * (0.5 + j);
		return computeDensityVariance(m, data, de, a, b, numEvalPoints, evalPoints);
	}

	/**
	 * Similar to computeDensityVariance but does it for a list of density
	 * estimators, and returns the results (integrated variance for each DE) in
	 * array #integVariance. NOTE: In case we would like to provide more results, we
	 * may want to return them in a list like listDE.
	 * 
	 * We may also want a method somewhere that plots the densities in a listDE ...
	 */
	public static void computeDensityVarianceListDE(int m, double[][] data, ArrayList<DensityEstimator> listDE,
			double a, double b, int numEvalPoints, double[] evalPoints, double[] integVariance) {
		integVariance = new double[listDE.size()];
		int deNumber = 0;
		for (DensityEstimator de : listDE) {
			integVariance[deNumber] = computeDensityVariance(m, data, de, a, b, numEvalPoints, evalPoints);
			deNumber++;
		}
	}

	/**
	 * Similar to computeDensityVariance but does it for a list of density
	 * estimators, and returns the results (integrated variance for each DE) in
	 * array integVariance. NOTE: In case we would like to provide more results, we
	 * may want to return them in a list like listDE, but for the results (?).
	 */
	public static void computeDensityVarianceListDE(int m, double[][] data, ArrayList<DensityEstimator> listDE,
			double a, double b, int numEvalPoints, double[] integVariance) {
		double evalPoints[] = new double[numEvalPoints]; // Points at which the density will be evaluated
		double delta = (b - a) / (numEvalPoints);
		for (int j = 0; j < numEvalPoints; j++)
			evalPoints[j] = a + delta * (0.5 + j);
		computeDensityVarianceListDE(m, data, listDE, a, b, numEvalPoints, evalPoints, integVariance);
	}

	/**
	 * For a \ref MonteCarloModelDensityKnown this method estimates the density and
	 * computes the ISB based on \f$m\f$ independently generated simulations stored
	 * in \a data. This is done by numerically integrating the squared bias with
	 * integration nodes \a evalPoints over the interval on which \a de is defined.
	 * 
	 * 
	 * This is particularly useful for testing a density estimator with a toy model.
	 * 
	 * @param model
	 *            the underlying model.
	 * @param data
	 *            matrix in which the realizations of the model are stored.
	 * @param de
	 *            density estimator which is used.
	 * @param evalPoints
	 *            the integration nodes.
	 * @return the empirical ISB
	 */
	public static double computeDensityISB(MonteCarloModelDensityKnown model, double[][] data, DensityEstimator de,
			double[] evalPoints) {

		int numEvalPoints = evalPoints.length;
		int m = data.length;

		double a = de.getMin();
		double b = de.getMax();

		double meanDens[] = new double[numEvalPoints];
		double sqBiasDens[] = new double[numEvalPoints]; // squared bias

		double estimDens[] = new double[numEvalPoints]; // Value of the density at evalPoints

		double x, y, z;

		for (int rep = 0; rep < m; rep++) {
			// Estimate the density and evaluate it at evalPoints
			de.constructDensity(data[rep]);
			de.evalDensity(evalPoints, estimDens);

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
	 * Computes the estimated ISB for the \ref DensityEstimatorBandwidthBased \a de,
	 * i.e. \f$\textrm{ISB}\approx B h^{\alpha}\f$ Note that \f$B,h \f$, and \f$
	 * \alpha\f$ need to be set for \a de before the ISB can be calculated.
	 * 
	 * @param de
	 *           the density estimator which is used.
	 * @return the estimated ISB
	 */
	public static double computeDensityISB(DensityEstimatorBandwidthBased de) {
		return de.getB() * Math.pow(de.getH(), de.getAlpha());
	}

	/**
	 * Same as {@link #computeDensityISB(DensityEstimatorBandwidthBased)} but with additional dummy-arguments to match
	 * the structure of {@link #computeDensityISB(MonteCarloModelDensityKnown, double[][], DensityEstimator, double[])}
	  * @param model
	 *            the underlying model.
	 * @param data
	 *            matrix in which the realizations of the model are stored.
	 * @param de
	 *            density estimator which is used.
	 * @param evalPoints
	 *            the integration nodes.
	 * @return the estimated ISB
	 */
	public static double computeDensityISB(MonteCarloModelDouble model, double[][] data, DensityEstimatorBandwidthBased de,
			double[] evalPoints) {
		return computeDensityISB(de);
		
	}
	
	/**
	 * For a \ref MonteCarloModelDensityKnown this method estimates the density and
	 * computes the MISE based on \f$m\f$ independently generated simulations stored
	 * in \a data. This is done by numerically integrating the mean square error with
	 * integration nodes \a evalPoints over the interval on which \a de is defined.
	 * 
	 * 
	 * This is particularly useful for testing a density estimator with a toy model.
	 * 
	 * @param model
	 *            the underlying model.
	 * @param data
	 *            matrix in which the realizations of the model are stored.
	 * @param de
	 *            density estimator which is used.
	 * @param evalPoints
	 *            the integration nodes.
	 * @return the empirical MISE
	 */
	
	public static double computeDensityMISE(MonteCarloModelDensityKnown model, double[][] data, DensityEstimator de,
			double[] evalPoints) {
		
		int m = data.length;
		int numEvalPoints = evalPoints.length;
	
		double estimDens[] = new double[numEvalPoints]; // Value of the density
														// at evalPoints
		double mseDens[] = new double[numEvalPoints]; // mse at evalPoints


		double a = de.getMin();
		double b = de.getMax();

		double y;
		for (int rep = 0; rep < m; rep++) {
			// Estimate the density and evaluate it at eval points
			de.constructDensity(data[rep]);
			de.evalDensity(evalPoints, estimDens);

			
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
	
	/**
	 * Computes the estimated MISE for the \ref DensityEstimatorBandwidthBased \a de.
	 * More precisely, it computes the empirical IV via {@link #computeDensityVariance(int, double[][], DensityEstimator, double, double, int, double[])}
	 * and an estimate for the ISB via {@link #computeDensityISB(DensityEstimatorBandwidthBased)} and adds them to 
	 * obtain an estimate for the MISE.
	 * @param data matrix in which the realizations of the model are stored.
	 * @param de density estimator which is used.
	 * @param evalPoints
	 *            the integration nodes.
	 * @return the estimated MISE.
	 */
	public static double computeDensityMISE(double[][] data, DensityEstimatorBandwidthBased de,
			double[] evalPoints) {
		double iv = computeDensityVariance(data.length,data,de,de.getMin(),de.getMax(),evalPoints.length,evalPoints);
		double isb = computeDensityISB(de);
		return iv + isb;
	}
	
	/**
	 * Same as #
	 * @param model
	 * @param data
	 * @param de
	 * @param evalPoints
	 * @return
	 */
	public static double computeDensityMISE(MonteCarloModelDouble model, double[][] data, DensityEstimatorBandwidthBased de,double[] evalPoints) {
		return computeDensityMISE(data,de,evalPoints);
	}
	
	/**
	 * Computes all traits for the \ref DensityEstimatorBandwidthBased \a de that
	 * are specified in the list \a traitList. The possible choices of traits are
	 * "isb", "iv", and "mise", with deviations due to letter-capitalization
	 * considered.
	 * 
	 * The IV is computed from the observations given in \a data, whereas the ISB
	 * only gives the estimate \f$\textrm{ISB}\approx B h^{\alpha}\f$. The MISE is
	 * computed by adding the values of the IV and the ISB. Note that in order to be
	 * able to compute the ISB and the MISE, several parameters of \a de need to be
	 * known and set in advance. Conversely, if only the ISB is to be computed, \a
	 * data and \a evalPoints can be empty.
	 * 
	 * @param traitsList
	 *            the traits which shall be computed.
	 * 
	 *            The computed traits are returned in an array in the same order as
	 *            they appear in \a traitsList.
	 * @param data
	 *            matrix containing the observations of \f$m\f$ independent
	 *            repetitions.
	 * @param de
	 *            the density estimator.
	 * @param evalPoints
	 *            the points at which the density estimator is evaluated to compute
	 *            the trait.
	 * @return the values of the desired traits.
	 */
	public static double[] computeDensityTraits(ArrayList<String> traitsList, double[][] data,
			DensityEstimatorBandwidthBased de, double[] evalPoints) {
		int t = traitsList.size();
		double[] traitsVals = new double[t];
		String traitName;
		// number of indep. repetitions. Not initialized here,
		// because for isb "data" can be empty
		int m;
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
			m = data.length; // number of indep. repetitions
			iv = computeDensityVariance(m, data, de, de.getMin(), de.getMax(), evalPoints.length, evalPoints);
			isb = de.getB() * Math.pow(de.getH(), de.getAlpha());
		}

		// go through all traits of the list and add value to output array. If "mise" is
		// not included in the list, also compute the traits.
		for (int s = 0; s < t; s++) {
			traitName = traitsList.get(s);
			switch (traitName.toLowerCase()) {
			case "iv":
				if (!containsMise) {
					m = data.length; // number of indep. repetitions
					iv = computeDensityVariance(m, data, de, de.getMin(), de.getMax(), evalPoints.length, evalPoints);
				}
				traitsVals[s] = iv;
				break;
			case "isb":
				if (!containsMise)
					isb = de.getB() * Math.pow(de.getH(), de.getAlpha());
				traitsVals[s] = isb;
				break;
			case "mise":
				traitsVals[s] = iv + isb;
			default:
				System.out.println("The trait " + traitName
						+ " cannot be interpreted. Supported traits are 'isb', 'iv', and 'mise'.");
			}
		}

		return traitsVals;
	}

	/**
	 * Same as
	 * {@link #computeDensityTraits(ArrayList, double[][], DensityEstimatorBandwidthBased, double[])}.
	 * The additional argument \a model clarifies that the exact sought density is
	 * not known in order to distinguish this method from
	 * {@link #computeDensityTraits(ArrayList, MonteCarloModelDensityKnown, double[][], DensityEstimatorBandwidthBased, double[])}.
	 * Hence, the values returned for the ISB and the MISE are only estimates.
	 * 
	 * @param traitsList
	 *            the traits which shall be computed.
	 * 
	 *            The computed traits are returned in an array in the same order as
	 *            they appear in \a traitsList.
	 * @param model
	 *            the model from which the observations were generated.
	 * @param data
	 *            matrix containing the observations of \f$m\f$ independent
	 *            repetitions.
	 * @param de
	 *            the density estimator.
	 * @param evalPoints
	 *            the points at which the density estimator is evaluated to compute
	 *            the trait.
	 * @return the values of the desired traits.
	 */
	public static double[] computeDensityTraits(ArrayList<String> traitsList, MonteCarloModelDouble model,
			double[][] data, DensityEstimatorBandwidthBased de, double[] evalPoints) {
		return computeDensityTraits(traitsList, data, de, evalPoints);
	}

	/**
	 * 
	 * Computes all traits for the \ref DensityEstimator \a de that are specified in
	 * the list \a traitList. The possible choices of traits are "isb", "iv", and
	 * "mise", with deviations due to letter-capitalization considered.
	 * 
	 * Since this method takes a \a model of type \ref MonteCarloModelDensityKnown,
	 * all the traits are computed empirically instead of merely estimated. This
	 * makes it particularly useful for testing estimators with toy models.
	 * 
	 * If \a traitsList contains only 1 element, then this method calls the usual
	 * method to compute the respective trait. Otherwise, it computes the traits in
	 * parallel.
	 * 
	 * The computed traits are returned in an array in the same order as they appear
	 * in \a traitsList.
	 *
	 * @param traitsList
	 *            the traits which shall be computed.
	 * 
	 * @param model
	 *            the model from which the observations were generated.
	 * @param data
	 *            matrix containing the observations of \f$m\f$ independent
	 *            repetitions.
	 * @param de
	 *            the density estimator.
	 * @param evalPoints
	 *            the points at which the density estimator is evaluated to compute
	 *            the trait.
	 * @return the values of the desired traits.
	 */
	public static double[] computeDensityTraits(ArrayList<String> traitsList, MonteCarloModelDensityKnown model,
			double[][] data, DensityEstimatorBandwidthBased de, double[] evalPoints) {

		double x, y, z;
		double dens;
		int numEvalPoints = evalPoints.length;
		int m = data.length;

		double estimDens[] = new double[numEvalPoints]; // Value of the density
														// at those points
		double meanDens[] = new double[numEvalPoints]; // Average value over rep
														// replicates
		double varDens[] = new double[numEvalPoints]; // Variance at each
														// evaluation point
		double mseDens[] = new double[numEvalPoints]; // MSE at each evaluation
														// point

		double a = de.getMin();
		double b = de.getMax();

		for (int rep = 0; rep < m; rep++) {
			// Estimate the density and evaluate it at evalPoints
			de.constructDensity(data[rep]);
			de.evalDensity(evalPoints, estimDens);

			// Update the empirical mean, sum of squares, and mse of
			// observations at each evaluation point.
			for (int j = 0; j < numEvalPoints; j++) {
				x = estimDens[j];
				y = x - meanDens[j];
				dens = model.density(evalPoints[j]);
				z = x - dens;

				meanDens[j] += y / (double) (rep + 1);
				varDens[j] += y * (x - meanDens[j]);
				mseDens[j] += z * z;
			}

		}

		double iv = 0.0;
		double mise = 0.0;
		// double isb = 0.0;
		for (int i = 0; i < numEvalPoints; ++i) {
			iv += varDens[i];
			mise += mseDens[i];
			// isb += biasDens[i];
		}

		double fact = (b - a) / ((double) (numEvalPoints * (m + 1.0)));
		iv *= fact;
		mise *= (b - a) / ((double) (numEvalPoints * m));

		double[] result = new double[3];
		result[0] = iv;
		result[1] = mise - iv;
		result[2] = mise;
	}

}
