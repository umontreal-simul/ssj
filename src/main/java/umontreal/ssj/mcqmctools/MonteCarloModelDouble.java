package umontreal.ssj.mcqmctools;
import umontreal.ssj.rng.*;

/**
 * An interface for a simulation model for which Monte Carlo (MC) and RQMC
 * experiments are to be performed. This interface is used by the classes
 * MCExperiment and RQMCExperimentDensity, among others, to run the model.  
 * The interface assumes that the output (sample performance) from the model 
 * is a real-valued random variable X. 
 * We could have asked @ref simulate to return the performance to avoid a 
 * separate call to @ref getPerformance, but we preferred not,
 * because in some situations one may not always need the performance computed
 * by @ref getPerformance, which might be costly to compute for nothing, 
 * but sometimes only some other output information.
 * This also applies to @ref MonteCarloModel.
*/

public interface MonteCarloModelDouble {

	// Optional
	// public void simulate ();

	// Simulates the model for one run
	public void simulate (RandomStream stream);

	// Recovers the realization of the performance measure $X$ for the last run.
	public double getPerformance();

	// Returns the max number of uniforms required to simulate the model.
	// Can be useful for RQMC simulation.
    // NOTE: Cannot always implement this, because often the model itself cannot know its dimension,
	// for example if the dimension depends on the stochastic process that is used... !!!!!
	public int getDimension ();

	// Returns a short description of the model and its parameters.
	public String toString();

}
