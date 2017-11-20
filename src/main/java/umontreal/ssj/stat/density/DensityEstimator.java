/**
 *  @author Pierre L'Ecuyer
 */
package umontreal.ssj.stat.density;

import java.util.ArrayList;

/**
 * An interface for a univariate density estimator over a finite interval [a,b].
 * The density is constructed from a set of n data points.
 * Some of these data points could be outside the interval [a,b].
 * The density can be evaluated at a single point x or on a grid of evaluation points.
 */
public interface DensityEstimator {
	
	
	/**
	 * From now on, the density will be estimated over the interval [a,b].
	 */
	public void setRange (double a, double b);
	 
	 
    /**	
     * Constructs a density estimator from the data points in vector #data.
     * @param data  the data points.
     */
	public void constructDensity (double[] data);
	

	/**
	 * Returns the value of the density evaluated at point x.
	 */
	public double evalDensity (double x);

	
    /**	
     * Returns in array #density the value of the density at the evaluation points in #evalPoints.
     * These two arrays must have the same size. 
     * @param evalPoints  the evaluation points
     * @param density     values of the density at these points
     */
	public void evalDensity (double[] evalPoints, double[] density);
	
	/**
	 * Takes data from previous simulations (m independent replications, with n data points each).
	 * Computes and returns the empirical integrated variance (IV) 
	 * for this density estimator, obtained by estimating the variance 
	 * at numEvalPoints equidistant points over [a,b] and summing up.
	 */
	public double computeDensityVariance (int n, int m,
			double[][] data, double a, double b, int numEvalPoints);
		
}
