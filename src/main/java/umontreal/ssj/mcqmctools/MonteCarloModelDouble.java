package umontreal.ssj.mcqmctools;
import umontreal.ssj.rng.*;

/**
 * An interface for a very simple simulation model for which Monte Carlo (MC) and RQMC
 * experiments are to be performed. This interface is used by the classes
 * MonteCarloExperiment and RQMCExperiment, among others, to run the model.  
 * This interface assumes that the simulation requires a single `RandomStream` and that
 * the output (sample performance) from the model is a real-valued random variable X. 
 * We could have asked @ref simulate to return the performance to avoid a 
 * separate call to @ref getPerformance, but we preferred not,
 * because in some situations one may not always need the performance computed
 * by @ref getPerformance, which might be costly to compute for nothing, 
 * but only some other output information.
 * This also applies to @ref MonteCarloModel.
*/

public interface MonteCarloModelDouble {

	// Optional
	// public void simulate ();

	/**
	 * Simulates the model for one run
	 */
	public void simulate (RandomStream stream);

	/** 
	 * Recovers and returns the realization of the performance measure, of type E.
	 */
	public double getPerformance();

	/** 
	 * Returns a short description of the model and its parameters.
	 */
	public String toString();

}
