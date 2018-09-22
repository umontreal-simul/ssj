package umontreal.ssj.mcqmctools;
import umontreal.ssj.rng.*;

/**
 * An interface for a very simple simulation model for which Monte Carlo (MC) and RQMC
 * experiments are to be performed. This interface is used by the classes
 * @ref MonteCarloExperiment and @ref RQMCExperiment, among others, to run the model.  
 * This interface assumes that the simulation requires a single `RandomStream` and that
 * the output (sample performance) from the model is a real-valued random variable @f$X@f$. 
 * We could have asked `simulate` to return the performance to avoid a 
 * separate call to `getPerformance`, but we decided not,
 * because in some situations one may not always need the performance computed
 * by `getPerformance`, which might be costly to compute for nothing, 
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
