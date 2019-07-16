package umontreal.ssj.mcqmctools.florian.examples;

import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stochprocess.BrownianMotionBridge;

/**
 * Implements the current value of a Lookback call option under a geometric
 * Brownian motion. This value is defined as the maximum of the
 * (discounted)values of the asset at several observation times. Here, the observation times are automatically selected as equidistant. The
 * Brownian motion is constructed via bridge sampling.
 * 
 * @author florian
 *
 */
public class LookbackOptionGBM implements MonteCarloModelDouble {

	int dim;
	double[] path;

	double s0, K, r, sigma, h;

	BrownianMotionBridge bm;

	/**
	 * Constructor
	 * @param dim number of observation times.
	 * @param s0 initial value of the asset.
	 * @param K strike price.
	 * @param r discount rate.
	 * @param sigma volatility.
	 * @param bm the brownian bridge.
	 */
	public LookbackOptionGBM(int dim, double s0, double K, double r, double sigma, BrownianMotionBridge bm) {
		this.dim = dim;
		this.s0 = s0;
		this.K = K;
		this.r = r;
		this.sigma = sigma;

		this.h = 1.0 / (double) dim;
		this.path = new double[dim + 1];
		this.bm = bm;
		this.bm.setObservationTimes(h, dim);
	}

	/**
	 * Constructor that automatically selects a brownian bridge with initial value 0, mean 0 and standard deviation 1.
	 * @param dim
	 * @param s0
	 * @param K
	 * @param r
	 * @param sigma
	 */
	public LookbackOptionGBM(int dim, double s0, double K, double r, double sigma) {
		this(dim, s0, K, r, sigma, new BrownianMotionBridge(0.0, 0.0, 1.0, new MRG32k3a()));
	}

	@Override
	public void simulate(RandomStream stream) {
		bm.setStream(stream);
		path = bm.generatePath();

	}

	@Override
	public double getPerformance() {

		double res = sigma * path[0];
		double temp;
		for (int j = 1; j <= dim; j++) {
			temp = (r - sigma * sigma * 0.5) * h * j + sigma * path[j];
			if (res < temp)
				res = temp;
		}

		return s0 * Math.exp(res);
	}

	public String toString() {
		return "Lookback Option with " + dim + "observations.";
	}
}
