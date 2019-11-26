package umontreal.ssj.mcqmctools.florian.examples;

import umontreal.ssj.mcqmctools.MonteCarloModelDoubleArray;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.rng.RandomStream;

public class BucklingStrengthVars implements MonteCarloModelDoubleArray {

	/**
	 * This class is intended to create and return several random variables related
	 * to the buckling strength @f$X@f$ of a steel plate as given by Schields and
	 * Zhang '16. This is used for estimating the density of @f$X@f$ with a
	 * conditional density estimator, for instance. The model is given by
	 * 
	 * \f[ X=\left( \frac{2.1}{\lambda}-\frac{0.9}{\lambda^2} \right) \left(
	 * 1-\frac{0.75\delta_0}{\lambda}\right) \left(1-\frac{2\eta t}{b} \right), \f]
	 * 
	 * 
	 * where \lambda = (b/t) \sqrt{\sigma_0/eta} and where the
	 * parameters @f$b,t,\sigma_0,E,\delta,\eta@f$ are random variables. More
	 * precisely, $t$ and $b_0$ follow a lognormal distribution and all the others
	 * are normally distributed.
	 * 
	 * The random variates that are returned by this class (in the stated order)
	 * are @f$2t/b, \lambda, \delta_0, \eta@f$.
	 * 
	 * @author florian
	 *
	 */
	double[] mu;
	double[] sigma;
	double[] performance; // 2t/b, lambda , delta0, eta

	/**
	 * Constructor to pass the parameters of the distributions for all random
	 * variables involved. The order in which they should be passed is @f$b, t,
	 * sigma0, E, delta0, eta@f$. Note that for the lognormal variables also the
	 * means and standard deviations are required, not the distribution parameters!
	 * 
	 * @param mu    array with the means.
	 * @param sigma array with the standard deviations.
	 */
	public BucklingStrengthVars(double[] mu, double[] sigma) {
		this.mu = mu;
		this.sigma = sigma;
		performance = new double[4];
	}

	private double transformMu(double mu, double sigma) {
		return (Math.log(mu) - 0.5 * Math.log(sigma * sigma / (mu * mu) + 1.0));
	}

	private double transformSigma(double mu, double sigma) {
		return (Math.sqrt(Math.log(1.0 + sigma * sigma / (mu * mu))));
	}

	/*
	 * b,t,delta0,eta,lambda
	 */
	@Override
	public void simulate(RandomStream stream) {
		// 2t/b
		performance[0] = 2.0 * Math.exp(
				NormalDist.inverseF(transformMu(mu[1], sigma[1]), transformSigma(mu[1], sigma[1]), stream.nextDouble()))
				/ NormalDist.inverseF(mu[0], sigma[0], stream.nextDouble()); // b
		// lambda
		performance[1] = 2.0
				* Math.sqrt(Math.exp(NormalDist.inverseF(transformMu(mu[2], sigma[2]), transformSigma(mu[2], sigma[2]),
						stream.nextDouble())) / NormalDist.inverseF(mu[3], sigma[3], stream.nextDouble()))
				/ performance[0];
		// delta0
		performance[2] = NormalDist.inverseF(mu[4], sigma[4], stream.nextDouble());
		// eta
		performance[3] = NormalDist.inverseF(mu[5], sigma[5], stream.nextDouble());// eta
	}

	@Override
	public double[] getPerformance() {

		return performance;
	}

	@Override
	public int getPerformanceDim() {
		return performance.length;
	}

}
