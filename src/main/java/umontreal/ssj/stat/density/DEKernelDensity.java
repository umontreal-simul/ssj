package umontreal.ssj.stat.density;

import umontreal.ssj.probdist.ContinuousDistribution;

/**
 * This class provides methods to construct a kernel density estimator (KDE) for
 * univariate densities from a set of \f$n\f$ individual observations \f$x_0, …,
 * x_{n-1}\f$, and to evaluate it at a single point or at a set of selected
 * evaluation points. The observations can be collected data or realizations of a \ref
 * umontreal.ssj.mcqmctools.MonteCarloModelDouble, for instance.
 * 
 * The KDE takes a fixed bandwidth \f$ h>0\f$ as well as a kernel function
 * \f$K\f$, which is also referred to as kernel density. The kernel density
 * should be non-negative and integrate to one. For \f$x\in[a,b]\f$, the KDE
 * itself is defined as
 * 
 * @anchor REF_stat_density_DEKernelDensity_KDE \f[ \hat{f}_{n}(x) =
 *         \hat{f}_{n,h}(x) = \frac{1}{nh} \sum_{i = 0}^{n-1} K\left( \frac{x -
 *         x_i}{h} \right). \tag{KDE} \f]
 * 
 * This class also offers static methods so that the user can simply evaluate the density based on
 * a set of observations without having to construct a KDE.
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
	 * Constructs a KDE from the observations \a data.
	 * 
	 * @param data
	 *            the observations.
	 */
	public DEKernelDensity(double[] data) {
		this.data = data;
	}

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
	 * Constructs a KDE with the kernel function \a kernel from the observations in
	 * \a data.
	 * 
	 * @param kernel
	 *            the kernel density function.
	 * @param data
	 *            the observations.
	 */
	public DEKernelDensity(ContinuousDistribution kernel, double[] data) {
		setKernel(kernel);
		this.data = data;
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
	 * Constructs a KDE with bandwidth \a h from the observations \a data.
	 * 
	 * @param h
	 *            the bandwidth
	 * @param data
	 *            the observations
	 */
	public DEKernelDensity(double h, double[] data) {
		this(h);
		this.data = data;
	}

	/**
	 * Constructs a KDE with the kernel function \a kernel and bandwidth \a h.
	 * 
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
	 * Constructs a KDE with the kernel function \a kernel and bandwidth \a h from
	 * the observations \a data
	 * 
	 * @param kernel
	 *            the kernel density function
	 * @param h
	 *            the bandwidth
	 * @param data
	 *            the observations
	 */
	public DEKernelDensity(ContinuousDistribution kernel, double h, double[] data) {
		this(kernel, data);
		setH(h);
	}

	/**
	 * Sets a new set of observations.
	 * 
	 * @param data
	 *            the desired observations.
	 */
	public void setData(double[] data) {
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
	 * Sets the threshold-level \f$\varepsilon\f$ for the evaluation of the density
	 * to \a eps. Summands w.r.t. \f$i\f$ in (
	 * {@link REF_stat_density_DEKernelDensity_KDE KDE} ) smaller than this value
	 * are considered negligible.
	 * 
	 * @param eps
	 *            the threshold-level for evaluation.
	 */
	public void setEps(double eps) {
		this.eps = eps;
	}

	/**
	 * Evaluates the KDE at the evaluation point \a x. Each summand w.r.t. \f$i\f$
	 * in ( {@link REF_stat_density_DEKernelDensity_KDE KDE} ) is only considered if
	 * it is larger than \f$\varepsilon\f$ . For this method the kernel function
	 * \f$K\f$ is assumed to be unimodal, i.e. increasing and then decreasing.
	 * 
	 * @param x
	 *            the evaluation point.
	 * @return the value of the KDE at \a x.
	 */
	public double evalDensity(double x) {
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
		sum = term; // The first significant term.
		for (int i = imin + 1; (i < n) && ((term > eps) || (data[i] < x)); i++) // Data
																				// indexed
			// by i.
			sum += (term = kernel.density((x - data[i]) * invh));
		density = sum * invhn;
		return density;
	}

	/**
	 * Evaluates the KDE at each of the evaluation points \a evalPoints and returns
	 * the results in an array. Each summand w.r.t. \f$i\f$ in (
	 * {@link REF_stat_density_DEKernelDensity_KDE KDE} ) is only considered if it
	 * is larger than \f$\varepsilon\f$ .
	 * 
	 * For this method the kernel function \f$K\f$ is assumed to be unimodal, i.e.
	 * increasing and then decreasing, and that the points in \a evalPoints are
	 * sorted in increasing order. This allows this method to avoid looping over all
	 * \f$i\f$ for each evaluation point by remembering that some summands have
	 * already been deemed too small.
	 * 
	 * @param evalPoints
	 *            the evaluation points.
	 * 
	 * @return the value of the KDE at \a evalPoints.
	 */
	@Override
	public double[] evalDensity(double[] evalPoints) {
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

	// STATIC METHODS:

	/**
	 * Evaluates the KDE with kernel density function\a kernel, bandwidth \a h which
	 * is defined by the observations \a data at the evaluation point \a x. Each
	 * summand w.r.t. \f$i\f$ in ( {@link REF_stat_density_DEKernelDensity_KDE KDE}
	 * ) is only considered if it is larger than \a eps . For this method the kernel
	 * function \f$K\f$ is assumed to be unimodal, i.e. increasing and then
	 * decreasing.
	 * 
	 * @param x
	 *            the evaluation point.
	 * @param kernel
	 *            the kernel density function.
	 * @param h
	 *            the bandwidth.
	 * @param data
	 *            the observations.
	 * @param eps
	 *            the threshold level.
	 * @return the KDE defined by the above parameters evaluated at \a x.
	 */
	public static double evalDensity(double x, ContinuousDistribution kernel, double h, double[] data, double eps) {
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
	 * Evaluates the KDE with kernel density function\a kernel, bandwidth \a h which
	 * is defined by the observations \a data at each of the evaluation points \a
	 * evalPoints and returns the results in an array. Each summand w.r.t. \f$i\f$
	 * in ( {@link REF_stat_density_DEKernelDensity_KDE KDE} ) is only considered if
	 * it is larger than \f$\varepsilon\f$ .
	 * 
	 * For this method the kernel function \f$K\f$ is assumed to be unimodal, i.e.
	 * increasing and then decreasing, and that the points in \a evalPoints are
	 * sorted in increasing order. This allows this method to avoid looping over all
	 * \f$i\f$ for each evaluation point by remembering that some summands have
	 * already been deemed too small.
	 * 
	 * @param evalPoints
	 *            the evaluation points.
	 * @param kernel
	 *            the kernel density function.
	 * @param h
	 *            the bandwidth.
	 * @param data
	 *            the observations.
	 * @param eps
	 *            the threshold level.
	 * @return the KDE defined by the above parameters evaluated at the evaluation
	 *         points \a evalPoints.
	 */
	public static double[] evalDensity(double[] evalPoints, ContinuousDistribution kernel, double h, double[] data,
			double eps) {
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
	 * Assume that we have \f$m\f$ independent realizations of the underlying model.
	 * For each such realization this method evaluates a KDE with kernel density
	 * function\a kernel, bandwidth \a h at the points from \a evalPoints. The
	 * independent realizations are passed via the 2-dimensional \f$m\times n\f$
	 * array \a data, where \f$n\f$ denotes the number of observations per
	 * realization. Hence, its first index identifies the independent realization
	 * while its second index identifies a specific observation of this realization.
	 * 
	 * The result is returned as a \f$m\times k\f$ matrix, where \f$k \f$ is the
	 * number of evaluation points, i.e., the length of \a evalPoints. The first
	 * index, again, identifies the independent realization whereas the second index
	 * corresponds to the point of \a evalPoints at which the KDE was
	 * evaluated.
	 * 
	 * For this method the kernel function \f$K\f$ is assumed to be unimodal, i.e.
	 * increasing and then decreasing, and that the points in \a evalPoints are
	 * sorted in increasing order. For a specific observation, this allows this
	 * method to avoid looping over all \f$i\f$ for each evaluation point by
	 * remembering that some summands have already been deemed too small.
	 * 
	 * @param evalPoints
	 *            the evaluation points.
	 * @param kernel
	 *            the kernel density function.
	 * @param h
	 *            the bandwidth.
	 * @param data
	 *            the observations.
	 * @param eps
	 *            the threshold level.
	 * @return the KDE for each realization evaluated at \a evalPoints.
	 */
	public static double[][] evalDensity(double[] evalPoints, ContinuousDistribution kernel, double h, double[][] data,
			double eps) {
		int m = data.length;
		int k = evalPoints.length;
		double[][] density = new double[m][k];
		for (int rep = 0; rep < m; rep++)
			density[rep] = evalDensity(evalPoints, kernel, h, data[rep], eps);
		return density;

	}

}
