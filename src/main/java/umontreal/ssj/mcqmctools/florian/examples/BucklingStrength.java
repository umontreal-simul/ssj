package umontreal.ssj.mcqmctools.florian.examples;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.LMScrambleShift;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.RQMCPointSet;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.mcqmctools.RQMCExperiment;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;

/**
 * This class implements the buckling strength @f$X@f$ of a steel plate as given
 * by Schields and Zhang '16, i.e.,
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
 * @author florian
 *
 */

public class BucklingStrength implements MonteCarloModelDouble {

	/*
	 * Order of parameters: b, t, sigma0, E, delta0, eta
	 */
	double[] mu;
	double[] sigma;
	double performance;

	/**
	 * Constructor to pass the parameters of the distributions for all random
	 * variables involved. The order in which they should be passed is @f$b, t,
	 * sigma0, E, delta0, eta@f$. Note that for the lognormal variables also the
	 * means and standard deviations are required, not the distribution parameters!
	 * 
	 * @param mu    array with the means.
	 * @param sigma array with the standard deviations.
	 */
	public BucklingStrength(double[] mu, double[] sigma) {
		this.mu = mu;
		this.sigma = sigma;
	}

	private double transformMu(double mu, double sigma) {
		return (Math.log(mu) - 0.5 * Math.log(sigma * sigma / (mu * mu) + 1.0));
	}

	private double transformSigma(double mu, double sigma) {
		return (Math.sqrt(Math.log(1.0 + sigma * sigma / (mu * mu))));
	}

	public void simulate(RandomStream stream) {
		double b = NormalDist.inverseF(mu[0], sigma[0], stream.nextDouble());
		double t = Math.exp(NormalDist.inverseF(transformMu(mu[1], sigma[1]), transformSigma(mu[1], sigma[1]),
				stream.nextDouble()));
		double lambda = b / t * Math.sqrt(Math.exp(
				NormalDist.inverseF(transformMu(mu[2], sigma[2]), transformSigma(mu[2], sigma[2]), stream.nextDouble()))
				/ NormalDist.inverseF(mu[3], sigma[3], stream.nextDouble()));

		performance = (2.1 / lambda - 0.9 / (lambda * lambda))
				* (1.0 - 0.75 * NormalDist.inverseF(mu[4], sigma[4], stream.nextDouble()) / lambda)
				* (1.0 - 2.0 * t * NormalDist.inverseF(mu[5], sigma[5], stream.nextDouble()) / b);

	}

	public double getPerformance() {
		return performance;
	}

	public String toString() {
		return "BucklingStrength";
	}

	
}
