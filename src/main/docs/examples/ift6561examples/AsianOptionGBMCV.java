package ift6561examples;

import umontreal.ssj.stochprocess.*;
import umontreal.ssj.mcqmctools.*;
import umontreal.ssj.probdist.NormalDist;

/**
 * This class represents an Asian Option based on a geometric Brownian motion,
 * for which we use the payoff under a geometric average (instead of the arithmetic
 * average) as a control variate.  It is a subclass of `AsianOption`.
 */

public class AsianOptionGBMCV extends AsianOption implements MonteCarloModelCV {
	
	double expectedGeo; // Expected value of the CV (payoff under geometric average).

	/**
	 * Array <TT>obsTimes[0..d+1]</TT> must contain <TT>obsTimes[0] = 0</TT>,
	 * plus the <SPAN CLASS="MATH"><I>d</I></SPAN> observation times.
	 * 
	 */
	public AsianOptionGBMCV (double r, int d, double[] obsTimes, double strike) {
		super (r, d, obsTimes, strike);
	}

	public AsianOptionGBMCV (GeometricBrownianMotion sp, double r, int d,
			double[] obsTimes, double strike) {
		super (r, d, obsTimes, strike);
		setProcess (sp);  // Must make sure we compute expectedGeo.
	}

	public AsianOptionGBMCV (double r, int d, double T1, double T, double strike) {
		super (r, d, T1, T, strike);
	}
	
	/**
	 * Set (or reset) the process to `sp', which much be a GBM in this class.
	 * Also resets the observation times of `sp` to those of this `AsianOption` object.
	 */
	public void setProcess(GeometricBrownianMotion sp) {
		// Reset the process to sp. Assumes that obsTimes have been set.
		priceProcess = sp;
		sp.setObservationTimes(obsTimes, d);
		computeExpectedGeo();
	}
	
	/**
	 * Computes the expected value of geometric CV.
	 */
	public double computeExpectedGeo() {
		// First, get the underlying BM and its parameters.
		BrownianMotion bm = ((GeometricBrownianMotion) priceProcess).getBrownianMotion();
		double s0 = priceProcess.getX0();  // Initial value of the GBM.
        double mu = bm.getMu();
        double sigma = bm.getSigma();
       	double my = 0;
		double s2y = 0;
		for (int j = 1; j <= d; j++) {
			my += obsTimes[j];
			s2y += (obsTimes[j] - obsTimes[j - 1]) * (d - j + 1) * (d - j + 1);
		}
		my = Math.log(s0) + my * mu / d;
		s2y *= sigma * sigma / (d * d);
		double dd = (-Math.log(strike) + my) / Math.sqrt(s2y);
		return expectedGeo = discount
		        * (Math.exp(my + 0.5 * s2y) * NormalDist.cdf01(dd + Math.sqrt(s2y))
		                - strike * NormalDist.cdf01(dd));
	}

	/**
	 * Returns the expected value of geometric CV (must have been computed before).
	 * 
	 */
	public double getExpectedGeo() {
		return expectedGeo;
	}

	/**
	 * Computes and returns discounted payoff for geometric average. Assumes that the path has been
	 * generated.
	 * 
	 */
	public double getPayoffGeo() {
		// Computes and returns discounted payoff for geometric average.
		// Assumes that the path has been generated.
		double[] pathBM = ((GeometricBrownianMotion) priceProcess).getBrownianMotion().getPath();
		double average = 0.0;      // Average over BM sample path.
		for (int j = 1; j <= d; j++)
			average += pathBM[j];
		average /= d;
		average = path[0] * Math.exp(average);
		if (average > strike)
			return discount * (average - strike);
		else
			return 0.0;
	}
	
	// Recovers the realizations of the control variates for the the last run.
	// Returns the geometric average minus its expectation.  Assumes sp is a GBM process!!!
	public double[] getValuesCV() {
		double[] cv = new double[1];
        cv[0] = (getPayoffGeo() - expectedGeo); 
        // System.out.println(expectedGeo + "   " + cv[0]);
        return cv;
	}

	// Returns the number of control variates.
	public int getNumberCV() {
		return 1;
	}

	public String toString() {
		return "Asian option model with " + d + " observation times,"
				+ "under a GBM process,\n with payoff for geometric average as a control variate.";
	}
}
