package umontreal.ssj.mcqmctools.florian.examples;

import umontreal.ssj.mcqmctools.MonteCarloModelDoubleArray;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stochprocess.BrownianMotion;
import umontreal.ssj.stochprocess.BrownianMotionBridge;

/**
 * Class to generate several random variates involved in the estimation of the
 * current value of a Lookback call option under a geometric Brownian motion,
 * which can be used, e.g., for density estimation with a conditional density
 * estimator. The variates that are generated are the Brownian Motion at several
 * given observation times, including the initial value. Here, the observation times are automatically
 * selected as equidistant. The Brownian motion is constructed via bridge
 * sampling.
 * 
 * @remark **Florian:** What this class does, is basically computing the path of
 *         a standard Brownian motion with bridge sampling.
 * 
 * @author puchhamf
 *
 */
public class LookbackOptionGBMVars implements MonteCarloModelDoubleArray {

	int dim;
	double[] path;
	double h; // stepsize
	BrownianMotionBridge bm;

	/**
	 * Constructor.
	 * @param dim the number of observation times
	 * @param bm the Brownian motion
	 */
	public LookbackOptionGBMVars(int dim, BrownianMotionBridge bm) {
		this.dim = dim;
		path = new double[dim + 1];
		this.h = 1.0 / (double) dim;
		this.bm = bm;
		this.bm.setObservationTimes(h, dim);
	}

	/**
	 * Constructor that automatically selects a Brownian bridge with initial value 0, mean 0 and standard deviation 1.
	 * @param dim the number of observation times.
	 */
	public LookbackOptionGBMVars(int dim) {
		this(dim, new BrownianMotionBridge(0.0, 0.0, 1.0, new MRG32k3a()));
	}

	@Override
	public void simulate(RandomStream stream) {

		bm.setStream(stream);
		path = bm.generatePath();

	}

	@Override
	public double[] getPerformance() {
		return path;
	}

	@Override
	public int getPerformanceDim() {
		return (dim + 1);
	}

	public String toString() {
		return "Vars for lookback option with " + (dim) + " observations";
	}
}
