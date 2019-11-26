package umontreal.ssj.stat.density.florian.examples;

import umontreal.ssj.stat.density.ConditionalDensityEstimator;
import umontreal.ssj.stochprocess.BrownianMotionBridge;

/**
 * This class implements a likelihood-ratio estimator for the current value of a Lookback call option under a geometric
 * Brownian motion. This value is defined as the maximum of the
 * (discounted)values of the asset at several observation times. 
 * 
 * Realizations of the random variates involved can for instance be obtained with the class #LookbackOptionGBMVars.
 * @author florian
 *
 */
public class LRLookbackOptionGBM extends ConditionalDensityEstimator {

	int dim; 

	double s0,  r,sigma, a, b, h, K;


	/**
	 * Constructor
	 * @param dim number of observation times.
	 * @param s0 initial value of the asset.
	 * @param K strike price.
	 * @param r discount rate.
	 * @param sigma the volatility.
	 */
	public LRLookbackOptionGBM(int dim, double s0, double K, double r, double sigma) {
		this.dim = dim;
		this.s0 = s0;
		this.K = K;
		this.sigma = sigma;
		this.r = r;
		this.h = 1.0 / (double) dim;

	}
	
	//note: data includes X(0)
	double g(double[] data) {
		double res = sigma * data[0];
		double temp;
		for(int j = 1; j <= dim; j++) {
			temp = (r- sigma*sigma*0.5) * h * j + sigma * data[j];
			if(res < temp)
				res = temp;
		}
		return s0*Math.exp(res);
	}
	
	@Override
	public double evalEstimator(double x, double[] data) {
		double res = 0.0;
		if (g(data) > x)
			res = -data[1]/(sigma * x * h);
		
		return res;
	}

}
