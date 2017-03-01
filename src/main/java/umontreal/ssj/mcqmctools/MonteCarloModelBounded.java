package umontreal.ssj.mcqmctools;

/**
 * An interface for a simulation model for which Monte Carlo (MC) and RQMC
 * experiments are to be performed. This interface is used by the classes
 * MCExperiment and RQMCExperimentDensity, among others, to run the model.  
 * The interface assumes that the output from the model is a random
 * variable X with some density, mostly over the interval [a,b], i.e., 
 * it can be neglected outside of that interval when estimating the density.
 */

public interface MonteCarloModelBounded extends MonteCarloModelDouble {

	// Returns the min value a for density estimation for this model.
	public double getMin();

	// Returns the max value b for density estimation for this model.
	public double getMax();

}
