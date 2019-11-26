package umontreal.ssj.mcqmctools.florian.examples;

import umontreal.ssj.mcqmctools.MonteCarloModelDoubleArray;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.rng.RandomStream;

/**
 * Class, that generates an of independent normal variables.This is frequently
 * used for density estimation with, e.g., a conditional density estimator,
 * since many simulation models make use of normal variates.
 * 
 * @author florian
 *
 */
public class MultiNormalIndependent implements MonteCarloModelDoubleArray {
	int dim;
	double[] mus;
	double[] sigmas;
	double[] performance;

	/**
	 * Constructor which passes the means and standard deviations for each
	 * coordinate of the random vector.
	 * 
	 * @param mus the means.
	 * @param sigmas the standard deviations.
	 */
	public MultiNormalIndependent(double[] mus, double[] sigmas) {
		this.mus = mus;
		this.sigmas = sigmas;
		dim = this.mus.length;
		performance = new double[dim];

	}

	@Override
	public void simulate(RandomStream stream) {
		for (int j = 0; j < dim; j++) {
			performance[j] = NormalDist.inverseF(mus[j], sigmas[j], stream.nextDouble());
		}
	}

	@Override
	public double[] getPerformance() {
		return performance;
	}

	@Override
	public int getPerformanceDim() {
		return dim;
	}

	@Override
	public String toString() {
		return "MultiNormalInd";
	}
}
