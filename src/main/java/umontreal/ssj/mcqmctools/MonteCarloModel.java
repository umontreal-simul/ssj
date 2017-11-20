package umontreal.ssj.mcqmctools;

import umontreal.ssj.rng.*;

/**
 * An interface for a simulation model for which Monte Carlo (MC) and RQMC
 * experiments are to be performed. It generalizes @ref MonteCarloModelDouble.
 * This interface allows the output (performance) from the model to be of arbitrary type.
 * 
 */

public interface MonteCarloModel<E> {

	// Optional
	// public void simulate ();

	// Simulates the model for one run
	public void simulate (RandomStream stream);

	// Recovers the realization of the performance measure $X$.
	public E getPerformance();

	// Returns the max number of uniforms required to simulate the model.
	// Can be useful for RQMC simulation.
    // NOTE: Cannot always implement this exactly, because often the model itself cannot know its dimension,
	// for example if the dimension depends on the stochastic process that is used!
	public int getDimension ();

	// Returns a short description of the model and its parameters.
	public String toString();

}
