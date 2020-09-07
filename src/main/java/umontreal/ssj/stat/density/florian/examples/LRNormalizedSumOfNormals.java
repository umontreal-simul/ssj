package umontreal.ssj.stat.density.florian.examples;

import java.util.Arrays;

import umontreal.ssj.stat.density.ConditionalDensityEstimator;

public class LRNormalizedSumOfNormals extends ConditionalDensityEstimator {

	/**
	 * This class implements a likelihood-ratio density estimator (LR) for a sum
	 * of @f$d@f$ normals with mean 0, which is normalized to have variance 1. The
	 * realizations of the random variables at which this LR  can
	 * be obtained with the class #MultiNormalIndependent.
	 * 
	 * 
	 * @author florian
	 *
	 */

	private int dimension;

	/**
	 * Constructor setting the number of summands @f$d@f$.
	 * @param dim
	 */
	public LRNormalizedSumOfNormals(int dim) {
		this.setDimension(dim);
	}

	/**
	 * @return the dimension
	 */
	public int getDimension() {
		return dimension;
	}

	/**
	 * @param dimension the dimension to set
	 */
	public void setDimension(int dimension) {
		this.dimension = dimension;
	}

	@Override
	public double evalEstimator(double x, double[] data) {
		double sum = 0.0;
		double sumSq = 0.0;
		for (int j = 0; j < getDimension(); j++) {
			sum += data[j];
		}
		if (sum > x)
			return 0;
		else {
			for (int j = 0; j < getDimension(); j++)
				sumSq += -data[j] * data[j];
			return ((sumSq + getDimension()) / x);
		}

	}

	public String toString() {
		return "LRNormalizedSumOf" + getDimension() + "normals";
	}

}
