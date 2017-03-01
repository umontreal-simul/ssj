package umontreal.ssj.mcqmctools;

/**
 * An interface for a simulation model for which Monte Carlo (MC) and RQMC
 * experiments are to be performed. This interface is used by the classes
 * MCExperimentDensity and RQMCExperimentDensity, among others, 
 * to run the model.  The interface assumes that the output from the model is a random
 * variable X with known density and pdf (for the purpose of the experiment).
 */

public interface MonteCarloModelDensityU01 extends MonteCarloModelDensityKnown {

	// Recovers the cdf of the last realization of the performance measure $X$.
	// Returns F(X).
	public double getValueU01();

	// Recovers the cdf of $X$ evaluated at x.
	public double cdf(double x);
	
}
