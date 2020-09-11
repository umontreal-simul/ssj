package umontreal.ssj.stat.density.florian.examples;

import java.util.Arrays;

import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.stat.density.ConditionalDensityEstimator;

/**
 * This class implements the CDE for a sum of @f$d@f$ normals with mean 0, which
 * is normalized to have variance 1. The realizations of the random variables at
 * which this CDE  can be obtained with the class
 * #MultiNormalIndependent.
 * 
 * This estimator consists of @f$d@f$ individual
 * estimators @f$\hat{f}_{j}@f$, @f$1 \leq j \leq d@f$, where the index
 * indicates which summand is hidden. For the final estimator we take a convex
 * combination of the above estimators, i.e.,
 * 
 * \f[\hat{f} = \sum_{j=1}^d \omega_j\hat{f}_{j},\f]
 * 
 * where the @f$\omega_j@f$ are positive weight whose sum is one.
 * 
 * If one wants to omit a certain estimator, one can set the
 * corresponding @f$\omega_j@f$ to a negative number. In this case the weight
 * will be treated as zero. In this case, the remaining weights should sum up to
 * one.
 * 
 * @author florian
 *
 */
public class CDENormalizedSumOfNormals extends ConditionalDensityEstimator {

	/*
	 * Array containing d weights, where the d is the number of estimators.
	 */
	private double[] weights;
	private int dimension;
	private double sigma;
	private double[] sigmas;

	/**
	 * Constructor, passing the number of normals in the sum, the weights, and the
	 * standard deviations of the normal variates.
	 * 
	 * @param dim     number of summands.
	 * @param weights the weights.
	 * @param sigmas  the standard deviations of the normals.
	 */
	public CDENormalizedSumOfNormals(int dim, double[] weights, double[] sigmas) {
		this.weights = new double[dim];
		this.weights = weights;
		this.setDimension(dim);
		this.sigmas = sigmas;
		this.sigma = 0.0;
		for (double s : sigmas)
			sigma += s * s;
		sigma = Math.sqrt(sigma);
	}

	/**
	 * Same as above but the standard deviations are not yet set.
	 * 
	 * @param dim
	 * @param weights
	 */
	public CDENormalizedSumOfNormals(int dim, double[] weights) {
		this.setWeights(weights);
		this.setDimension(dim);
	}

	/**
	 * Same as above, but with @f$\omega_j=1/d@f$ and the standard deviations are
	 * not yet set.
	 * 
	 * @param dim
	 */
	public CDENormalizedSumOfNormals(int dim) {
		weights = new double[dim];
		Arrays.fill(weights, 1.0 / (double) dim);
		this.setDimension(dim);

	}

	/**
	 * Constructor, for which, so far, only the weights are set.
	 * @param weights
	 */
	public CDENormalizedSumOfNormals(double[] weights) {
		this(weights.length, weights);
	}

	/**
	 * Setter for the weights.
	 * 
	 * @param weights
	 */
	public void setWeights(double[] weights) {
		for (int i = 0; i < dimension; i++)
			this.weights[i] = weights[i];
//		normalizeWeights();
	}

	/**
	 * @return the dimension
	 */
	public int getDimension() {
		return dimension;
	}

	/**
	 * @param dimension the number of summands
	 */
	public void setDimension(int dimension) {
		this.dimension = dimension;
	}

	@Override
	public double evalEstimator(double x, double[] data) {
//		for(double d : weights)
//		System.out.println("w = " + d);
		double val = 0.0;
		double sum;
		for (int leave = 0; leave < getDimension(); ++leave) {
			sum = 0.0;
			if (weights[leave] > 0) {

				for (int j = 0; j < getDimension(); ++j)
					if (j != leave)
						sum += data[j];
				val += weights[leave] * NormalDist.density01(x * sigma / sigmas[leave] - sum / sigmas[leave]) * sigma
						/ sigmas[leave];
			} // endif
		}
		return val;
	}

	public String toString() {
		return "CDENormalizedSumOf" + getDimension() + "Normals";
	}

}
