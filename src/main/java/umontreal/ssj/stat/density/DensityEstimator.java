/**
 * 
 */
package umontreal.ssj.stat.density;

/**
 * @author Lecuyer
 * 
 * An interface for a density estimator. 
 *
 */
public interface DensityEstimator {
	
	
	/**
	 * From now on, the density will be estimated over the interval [a,b].
	 * 
	 */
	public void setRange (double a, double b);
	 
	 
    /**	
     * Constructs a density estimator from 
     * the data points in vector x.
     * @param x  data points.
     */
	public void constructDensity (double[] x);
	

	/**
	 * Returns the value of the density evaluated at x.
	 */
	public double evalDensity (double x);

	
    /**	
     * Returns in array f the value of the density at the evaluation points in x.
     * These two arrays must have the same size. 
     * @param x  evaluation points
     * @param f  values of the density
     */
	public void evalDensity (double[] x, double[] f);
		

}
