package umontreal.ssj.mcqmctools;

/**
 * An MonteCarloModelDouble that also implements a vector of control variates.
 */

public interface MonteCarloModelCV extends MonteCarloModelDouble {

	// Recovers the realizations of the control variates for the the last run.
	public double[] getValuesCV();

	// Returns the number of control variates.
	public int getNumberCV();

}
