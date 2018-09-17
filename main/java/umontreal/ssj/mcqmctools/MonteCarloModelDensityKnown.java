package umontreal.ssj.mcqmctools;


/**
 * An interface for a simulation model for which Monte Carlo (MC) and RQMC
 * experiments are to be performed. This interface is used by the classes
 * MCExperimentDensity and RQMCExperimentDensity, among others, 
 * to run the model.  The interface assumes that the output from the model is a random
 * variable X with known density and pdf (for the purpose of the experiment).
 */

public interface MonteCarloModelDensityKnown extends MonteCarloModelDouble{

	// Recovers the density of X evaluated at x.
	public double density(double x);
	
	/**
	 * Recovers the cumulative density function (cdf) of \f$X\f$ evaluated at \f$x\f$
	 * @param x the point at which the cdf shall be evaluated.
	 * @return the cdf at \f$x\f$.
	 */
	public double cdf(double x);

}
