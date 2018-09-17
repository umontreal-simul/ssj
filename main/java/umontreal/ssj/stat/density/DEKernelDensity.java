package umontreal.ssj.stat.density;

import umontreal.ssj.probdist.ContinuousDistribution;

/**
 * This class provides methods to construct a kernel density estimator (KDE) for
 * univariate  densities from a set of \f$n\f$ individual observations
 * \f$X_0, …, X_{n-1}\f$, and to evaluate it at a single point or at a set of
 * selected evaluation points. The observations can be realizations of a \ref
 * umontreal.ssj.mcqmctools.MonteCarloModelDouble, for instance.
 * 
 * The KDE takes a fixed bandwidth \f$ h>0\f$ as well as a kernel function
 * \f$K\f$, which is also referred to as kernel density. The kernel density
 * should be non-negative and integrate to one. For \f$x\in[a,b]\f$, the KDE
 * itself is defined as
 * 
 * @anchor REF_stat_density_DEKernelDensity_KDE 
 *    \f[ \hat{f}_{n}(x) =
 *         \hat{f}_{n,h}(x) = \frac{1}{nh} \sum_{i = 0}^{n-1} K\left( \frac{x -
 *         X_i}{h} \right). \tag{KDE} \f]
 * 
 * 
 * @author puchhamf
 *
 */

public class DEKernelDensity extends DensityEstimator {

	/** 
	 * The kernel function \f$K\f$.
    */
	private ContinuousDistribution kernel;
	
	/**
	 * The bandwidth \f$h\f$.
	 */
	private double h;

	/**
	 * threshold value for density-evaluation. Summands w.r.t. \f$i\f$ in (
	 * {@link REF_stat_density_DEKernelDensity_KDE KDE} ) smaller than this value
	 * are considered negligible. It is set to \f$ 10^{-10} by default\f$.
	 */
	private double eps = 1.0E-10;

	/**
	 * Constructs a KDE with the kernel function \a kernel.
	 * 
	 * @param kernel
	 *            the kernel density function.
	 */
	public DEKernelDensity(ContinuousDistribution kernel) {
		setKernel(kernel);
	}

	/**
	 * Constructs a KDE with bandwidth \f$h\f$.
	 * 
	 * @param h
	 *            the bandwidth.
	 */
	public DEKernelDensity(double h) {
		setH(h);
	}

	/**
	 * Constructs a KDE with the kernel function \a kernel and bandwidth \a h.
	 * @param kernel
	 *            the kernel density function
	 * @param h
	 *            the bandwidth
	 */
	public DEKernelDensity(ContinuousDistribution kernel, double h) {
		this(kernel);
		setH(h);
	}

	/**
	 * Constructs a KDE with the kernel function \a kernel and bandwidth \a h.
	 * @param kernel
	 *            the kernel density function
	 * @param h
	 *            the bandwidth
	 */
	public DEKernelDensity(ContinuousDistribution kernel, double h, double[] data) {
		this(kernel);
		setH(h);
		this.data = data;
	}

    /**
	 * Sets the bandwidth to \a h.
	 * 
	 * @param h
	 *            the desired bandwidth.
	 */
	public void setH(double h) {
		this.h = h;
	}

	/**
	 * Gives the bandwidth \f$h\f$.
	 * 
	 * @return the bandwidth.
	 */
	public double getH() {
		return h;
	}

	/**
	 * Gives the kernel density function \f$K\f$.
	 * 
	 * @return the kernel density function.
	 */
	public ContinuousDistribution getKernel() {
		return kernel;
	}

	/**
	 * Sets the kernel density function to \a kernel.
	 * 
	 * @param kernel
	 *            the kernel density function to be used.
	 */
	public void setKernel(ContinuousDistribution kernel) {
		this.kernel = kernel;
	}

	/**
	 * Gives the threshold-level \f$\varepsilon\f$ for the evaluation of the
	 * density. Summands w.r.t. \f$i\f$ in (
	 * {@link REF_stat_density_DEKernelDensity_KDE KDE} ) smaller than this value
	 * are considered negligible.
	 * 
	 * @return the threshold-level for evaluation.
	 */
	public double getEps() {
		return eps;
	}

	/**
	 * Sets the threshold-level \f$\varepsilon\f$  for the evaluation of the density to \a eps.
	 * Summands  w.r.t. \f$i\f$ in (
	 * {@link REF_stat_density_DEKernelDensity_KDE KDE} ) smaller than this value are considered
	 * negligible.
	 * 
	 * @param eps
	 *            the threshold-level for evaluation.
	 */
	public void setEps(double eps) {
		this.eps = eps;
	}

	/**
	 * Evaluates the KDE defined by the observations \a data at the evaluation point
	 * \a x. Each summand w.r.t. \f$i\f$ in (
	 * {@link REF_stat_density_DEKernelDensity_KDE KDE} ) is only considered if it
	 * is larger than \f$\varepsilon\f$ . For this method the kernel function \f$K\f$ is assumed
	 * to be unimodal, i.e. increasing and then decreasing.
	 * 
	 * @param x
	 *            the evaluation point.
	 * @param data
	 *            the observations of the underlying model.
	 * @return the value of the KDE at \a x.
	 */
	public double evalDensity(double x, double[] data) {
		double density;
		int n = data.length;
		double invh = 1.0 / h;
		double invhn = invh / (double) n;
		double sum = 0.0;
		double term; // A term to be added to the sum that defines the density
						// estimate.
		int imin = 0; // We know that the terms for i < imin do not contribute
						// significantly.

		term = kernel.density((x - data[imin]) * invh);
		while ((term < eps) && (imin < n - 1) && (data[imin] < x))
			term = kernel.density((x - data[++imin]) * invh);
		// System.out.println(imin);
		// pmin=imin;
		sum = term; // The first significant term.
		for (int i = imin + 1; (i < n) && ((term > eps) || (data[i] < x)); i++) // Data
																				// indexed
			// by i.
			sum += (term = kernel.density((x - data[i]) * invh));
		density = sum * invhn;
		// System.out.println(density);
		return density;
	}

	/**
	 * Same as #evalDensity(double, double[]) but with two placeholder arguments to
	 * implement the corresponding function demanded by the abstract superclass \ref
	 * umontreal.ssj.stat.density.DensityEstimator.
	 */
	@Override
	public double evalDensity(double x, double[] data, double a, double b) {
		return evalDensity(x, data);
	}

	/**
	 * Evaluates the KDE defined by the observations \a data at each of the
	 * evaluation points \a evalPoints and returns the results in an array. Each
	 * summand w.r.t. \f$i\f$ in ( {@link REF_stat_density_DEKernelDensity_KDE KDE}
	 * ) is only considered if it is larger than \f$\varepsilon\f$ .
	 * 
	 * For this method the kernel function \f$K\f$  is assumed to be unimodal, i.e.
	 * increasing and then decreasing, and that the points in \a evalPoints are
	 * sorted in increasing order. This allows this method to avoid looping over all
	 * \f$i\f$ for each evaluation point by remembering that some summands have
	 * already been deemed too small.
	 * 
	 * @param evalPoints
	 *            the evaluation points.
	 * @param data
	 *            the observations of the underlying model.
	 * @return the value of the KDE at \a evalPoints.
	 */

	public double[] evalDensity(double[] evalPoints, double[] data) {
		int k = evalPoints.length;
		double[] density = new double[k];
		int n = data.length;
		double invh = 1.0 / h;
		double invhn = invh / (double) n;
		double y;
		double sum = 0.0;
		double term; // A term to be added to the sum that defines the density
						// estimate.
		int imin = 0; // We know that the terms for i < imin do not contribute
						// significantly.
		for (int j = 0; j < k; j++) { // Evaluation points are indexed by j.
			y = evalPoints[j];
			term = kernel.density((y - data[imin]) * invh);
			while ((term < eps) && (imin < n - 1) && (data[imin] < y))
				term = kernel.density((y - data[++imin]) * invh);
			sum = term; // The first significant term.
			for (int i = imin + 1; (i < n) && ((term > eps) || (data[i] < y)); i++)
				// Data indexed by i.
				sum += (term = kernel.density((y - data[i]) * invh));
			density[j] = sum * invhn;
		}

		return density;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String toString() {
		String str = "KDE [h = " + h + ", Kernel: " + kernel.toString() + "]";
		return str;
	}

}
