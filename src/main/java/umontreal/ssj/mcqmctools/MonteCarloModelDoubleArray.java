package umontreal.ssj.mcqmctools;
import umontreal.ssj.rng.*;

/**
 * Similar to @ref MonteCarloModelDouble except that the returned performance is an array of real numbers.  
 * The dimension of that array must be returned by @ref getPerformanceDim().
*/

public interface MonteCarloModelDoubleArray {

	// Optional
	// public void simulate ();

	// Simulates the model for one run
	public void simulate (RandomStream stream);

	// Recovers the realization of the vector of performance measures for the last run.
	public double[] getPerformance();

	// Returns the dimension of the array of performance measures.
	public int getPerformanceDim();

	// Returns the max number of uniforms required to simulate the model.
	// Can be useful for RQMC simulation.
    // NOTE: Cannot always implement this, because often the model itself cannot know its dimension,
	// for example if the dimension depends on the stochastic process that is used... !!!!!
	public int getDimension ();

	// Returns a short description of the model and its parameters.
	public String toString();

}
