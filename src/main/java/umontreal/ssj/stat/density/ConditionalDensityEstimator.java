package umontreal.ssj.stat.density;

/**
 * This is an abstract class that implements an @ref
 * DensityEstimatorDoubleArray. It is designed for situations, where the
 * estimator @f$\hat{f}@f$ can be written in the form
 * 
 * \f[ \hat{f}(x) = \frac{1}{n}\sum_{i=0}^{n-1}g(x;\mathbf{X}_i), \f]
 * 
 * where @f$\mathbf{X}_i@f$ is a realization of a @f$t@f$-dimensional random
 * variable @f$\mathbf{X}_i@f$ and where @f$g@f$ is a real function. This is the
 * case for a conditional density estimator, or a likelihood-ratio estimator,
 * for instance.
 * 
 * @remark **Florian:** Maybe the name is misleading, because this estimator
 *         applies in many different settings.
 * @remark **TODO:** Add citations when papers for these estimators are published. 
 * 
 *         The evaluation of @f$g(x;\mathbf{X})@f$ at @f$x@f$ given @f$\mathbf{X}@f$
 *         is handled by the abstract function #evalEstimator.
 * 
 * @author florian
 *
 */

public abstract class ConditionalDensityEstimator extends DensityEstimatorDoubleArray {

	@Override
	public void setData(double[][] data) {
		this.data = new double[data.length][];
		for (int i = 0; i < data.length; i++) {
			this.data[i] = new double[data[i].length];
			for (int j = 0; j < data[i].length; j++)
				this.data[i][j] = data[i][j];
		}

	}

	@Override
	public double evalDensity(double x) {
		double dens = 0.0;
		int N = data.length;
		double Ninv = 1.0 / (double) N;
		for (int i = 0; i < N; i++) {
			dens += evalEstimator(x, data[i]);
			dens *= Ninv;
		}
		return dens;
	}

	@Override
	public double[] evalDensity(double[] x) {
		int k = x.length;
		double[] dens = new double[k];
		int N = data.length;
		double Ninv = 1.0 / (double) N;
		for (int j = 0; j < k; j++) {
			dens[j] = 0.0;
			for (int i = 0; i < N; i++) {
				dens[j] += evalEstimator(x[j], data[i]);

			}
			dens[j] *= Ninv;
		}
		return dens;
	}

	@Override
	public double[] evalDensity(double[] evalPoints, double[][] data) {
		setData(data);
		return evalDensity(evalPoints);
	}

	/**
	 * Evaluates the function @f$g@f$ at the point \a x and the realization
	 * of @f$\mathbf{X}@f$ given in the @f$t@f$-dimensional array \a data
	 * 
	 * @param x    the evaluation point.
	 * @param data the realization of @f$\mathbf{X}@f$.
	 * @return the function @f$g(x,\mathbf{X})@f$.
	 */
	public abstract double evalEstimator(double x, double[] data);

	@Override
	public String toString() {
		return "Conditional Density Estimator";
	}

}
