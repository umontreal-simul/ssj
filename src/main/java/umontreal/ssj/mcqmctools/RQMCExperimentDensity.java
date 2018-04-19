package umontreal.ssj.mcqmctools;

// import java.util.Arrays;
import java.util.ArrayList;

import umontreal.ssj.stat.density.DensityEstimator;
// import umontreal.ssj.stat.list.ListOfTallies;


/**
 * NOTE:  Since this class does not perform simulations, but just constructs density estimators
 * from an array of previously collected data, it should perhaps rather go in  @ref stat.density.DensityEstimator.
 * Also, we do not need a model here, only a and b.
 * Maybe these methods could go in  DensityEstimator.
 * 
 * This class extends @ref RQMCExperimentDouble by providing additional methods 
 * to estimate the integrated variance (IV) of a density estimator over a bounded interval [a,b], 
 * for a model that implements the MonteCarloModelBounded interface,
 * using the observations from m independent replicates of a simulation experiment with 
 * sample size n.  The latter may have been obtained in any way (e.g., by Monte Carlo or RQMC).
 * The density estimator is evaluated over a grid of #numEvalPoints points over the interval
 * [a,b] for each of the m replicates.  By computing the empirical average and variance at each evaluation
 * point, this permit one to estimate the density and the IV.
 */

/**
 * @author Pierre L'Ecuyer
 * 
 */
public class RQMCExperimentDensity {  //  extends RQMCExperimentDouble {


	/**
	 * Takes data from previous simulations (m independent replications, with n data points each).
	 * and a density estimator de.
	 * Computes and returns the empirical integrated variance (IV) 
	 * for this density estimator, obtained by estimating the variance 
	 * at numEvalPoints equidistant points over [a,b] and summing up.
	 * 
	 * TO DO: Maybe pass an array of arbitrary evaluation points instead?
	 *    Maybe ModelBounded is not essential:  could pass a and b directly perhaps?
	 */
	public static double computeDensityVariance (int n, int m,
			double[][] data, DensityEstimator de, double a, double b, int numEvalPoints) {

		double x, y;
		// TO DO:
		// If the density estimator is a histogram, here we may reset numEvalPoints to 
		// the number of bins of the histogram.
		double evalPoints[] = new double[numEvalPoints];  // Points at which the density will be evaluated
		double estimDens[] = new double[numEvalPoints];   // Value of the density at those points
		double meanDens[] = new double[numEvalPoints];    // Average value over the rep replicates
		double varDens[] = new double[numEvalPoints];     // Variance at each evaluation point
		// Arrays.fill(meanDens, 0.0);
		// Arrays.fill(varDens, 0.0);
		for (int rep = 0; rep < m; rep++) {
			// Estimate the density for this rep and evaluate it at the evaluation points
			de.constructDensity(data[rep]);
			de.evalDensity(evalPoints, estimDens);
	        // Update the empirical mean and sum of squares of centered observations at each evaluation point.
			for (int j = 0; j < numEvalPoints; j++) {
				x = estimDens[j];
				y = x - meanDens[j];
				meanDens[j] += y / (double) (rep+1);
				varDens[j] += y * (x - meanDens[j]);
			}
		}
		double sumVar = 0.0;
		for (int i = 0; i < numEvalPoints; ++i)
			sumVar += varDens[i];
		return sumVar * (b - a) / (numEvalPoints * (m - 1));   // Empirical integrated variance.
	}
	

		
	/**
	 * Similar to computeDensityVariance but does it for a list of density estimators,
	 * and returns the results (integrated variance for each DE) in array integVariance.
	 * NOTE: In case we would like to provide more results, we may want to return
	 * them in a list like listDE, but for the results (?).
	 * 
	 * We may also want a method that plots the densities?   
	 */
	public static void computeDensityVarianceListDE (int n, int m,
			double[][] data, ArrayList<DensityEstimator> listDE, double a, double b, int numEvalPoints, 
            double[] integVariance) {
        integVariance = new double[listDE.size()];
        int deNumber = 0;
		for (DensityEstimator de : listDE) {
            integVariance[deNumber] = computeDensityVariance (n, m,
        			data, de, a, b, numEvalPoints);
            deNumber++;
		}
	}
		
}
