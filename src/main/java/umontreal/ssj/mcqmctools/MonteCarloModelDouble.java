package umontreal.ssj.mcqmctools;
import umontreal.ssj.rng.*;

/**
 * An interface for a simulation model for which Monte Carlo (MC) and RQMC
 * experiments are to be performed. This interface is used by the classes
 * MCExperiment and RQMCExperimentDensity, among others, to run the model.  
 * The interface assumes that the output from the model is a random
 * variable X.  The reason we separate simulate and getValue  is to be able
 * to call getValue  several times without having to simulate the model each time.
 * In some models, there may also be different kinds of measures accessible by
 * other functions than getValue, or by subclasses that just redefine getValue.
 * Is this really useful or better ???
 */

public interface MonteCarloModelDouble {

	// Optional
	// public void simulate ();

	// Simulates the model for one run
	public void simulate (RandomStream stream);

	// Recovers the realization of the performance measure $X$.
	public double getPerformance();

	// Returns the max number of uniforms required to simulate the model.
	// Can be useful for RQMC simulation.
    // NOTE: Cannot always implement this, because often the model itself cannot know its dimension,
	// for example if the dimension depends on the stochastic process that is used... !!!!!
	public int getDimension ();

	// Returns a short description of the model and its parameters.
	public String toString();

}
