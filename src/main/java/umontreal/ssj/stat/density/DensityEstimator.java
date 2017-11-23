/**
 *  @author Pierre L'Ecuyer
 */
package umontreal.ssj.stat.density;

import java.util.ArrayList;

/**
 * This abstract class implements a univariate density estimator over a finite interval [a,b].
 * The density is constructed from a set of n data points.
 * Some of these data points could be outside the interval [a,b].
 * The density can be evaluated at a single point x or on a grid of evaluation points.
 */
public abstract class DensityEstimator {
	
    double a, b;  // The density can be evaluated over [a, b].
    
	/**
	 * From now on, the density will be estimated over the interval [a,b].
	 */
	public abstract void setRange (double a, double b);
	 
	 
    /**	
     * Constructs a density estimator from the data points in vector #data.
     * @param data  the data points.
     */
	public abstract void constructDensity (double[] data);
	

	/**
	 * Returns the value of the density evaluated at point x.
	 */
	public abstract double evalDensity (double x);

	
    /**	
     * Returns in array #density the value of the density at the evaluation points in #evalPoints.
     * These two arrays must have the same size. 
     * @param evalPoints  the evaluation points
     * @param density     values of the density at these points
     */
	public abstract void evalDensity (double[] evalPoints, double[] density);
	
	/**
	 * Takes data from previous simulations (m independent replications, with n data points each).
	 * Computes and returns the empirical integrated variance (IV) 
	 * for this density estimator, obtained by estimating the variance 
	 * at numEvalPoints equidistant points over [a,b] and summing up.
	 */

	public double computeDensityVariance (int n, int m,
			double[][] data, double a, double b, int numEvalPoints) {

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
			constructDensity(data[rep]);
			evalDensity(evalPoints, estimDens);
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
	 * This method is static, so we cannot declare it in the @ref DensityEstimator interface.
	 * 
	 * We may also want a method that plots the densities for pgfplots?   
	 */
	public static void computeDensityVarianceListDE (int n, int m,
			double[][] data, ArrayList<DensityEstimator> listDE, double a, double b, int numEvalPoints, 
            double[] integVariance) {
        integVariance = new double[listDE.size()];
        int deNumber = 0;
		for (DensityEstimator de : listDE) {
            integVariance[deNumber] = de.computeDensityVariance (n, m,
        			data, a, b, numEvalPoints);
            deNumber++;
		}
	}
		
}
