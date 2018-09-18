package umontreal.ssj.mcqmctools.anova;

import umontreal.ssj.mcqmctools.*;

/**
 * An interface for a simple simulation model for which Monte Carlo (MC) or RQMC
 * experiments are to be performed. It generalizes @ref MonteCarloModelDouble.
 * This interface allows the output (performance) from the model to be of arbitrary type @ref E.
 * It could be a scalar, an array, etc.
 */

public interface MonteCarloModelDoubleRQMC extends MonteCarloModelDouble {

	// Optional
	// public void simulate ();

	/**
	 * Returns the number of uniforms required to simulate this model.
	 */
	 public int getDimension ();

}
