package umontreal.ssj.stat.density;

import umontreal.ssj.probdist.ContinuousDistribution;
import umontreal.ssj.probdist.EmpiricalDist;

/**
 * This class provides methods to compute a kernel density estimator (KDE) for
 * univariate unimodal densities over an interval \f$[a,b]\f$ from a set of
 * \f$n\f$ individual observations \f$x_0, …, x_{n-1}\f$, and returns its value
 * at a set of selected points.
 * 
 * This estimator takes a fixed bandwidth \f$ h>0\f$ as well as a kernel function \f$K\f$, which is also referred to
 * as kernel density. The kernel density should be non-negative and integrate to one. For \f$x\in[a,b]\f$, the KDE itself is defined as
 * \f[ \hat{f}_{n}(x) = \hat{f}_{n,h}(x) = \frac{1}{nh} \sum_{i = 0}^{n-1} k\left( \frac{x - x_i}{h} \right). \f]
 * 
 * @author puchhamf
 *
 */

public class DEKernelDensity extends DEBandwidthBased {

	private EmpiricalDist dist;
	/** < container for the observations */
	private ContinuousDistribution kernel;
	/** < the kernel density function */
	private double eps = 1.0E-10; /**< threshold value for density-evaluation. Contributions smaller than this
									 * value are considered negligible.
									 */
	

	/**
	 * Constructs a KDE over the interval \f$[a,b]\f$.
	 * 
	 * @param a
	 *            the left boundary of the interval
	 * @param b
	 *            the right boundary of the interval
	 */
	public DEKernelDensity(double a, double b) {
		setRange(a, b);
		setAlpha(4.0);
	}

	/**
	 * Constructs a KDE over the interval \f$[a,b]\f$ with the kernel function
	 * <tt>kernel</tt>.
	 * 
	 * @param a
	 *            the left boundary of the interval
	 * @param b
	 *            the right boundary of the interval
	 * @param kernel
	 *            the kernel density function
	 */
	public DEKernelDensity(double a, double b, ContinuousDistribution kernel) {

		this(a, b);
		this.setKernel(kernel);

	}

	/**
	 * Constructs a KDE over the interval \f$[a,b]\f$ with bandwidth \f$h\f$ .
	 * 
	 * @param a
	 *            the left boundary of the interval
	 * @param b
	 *            the right boundary of the interval
	 * @param h
	 *            the bandwidth
	 */
	public DEKernelDensity(double a, double b, double h) {

		this(a, b);
		setH(h);

	}

	/**
	 * Constructs a KDE over the interval \f$[a,b]\f$ with the kernel function
	 * \p kernel and bandwidth \f$h\f$.
	 * 
	 * @param a
	 *            the left boundary of the interval
	 * @param b
	 *            the right boundary of the interval
	 * @param h
	 *            the bandwidth
	 * @param kernel
	 *            the kernel density function
	 */
	public DEKernelDensity(double a, double b, double h, ContinuousDistribution kernel) {

		this(a, b,kernel);
		setH(h);
	}

	/**
	 * Constructs a KDE over the interval \f$[a,b]\f$ with the kernel function
	 * \p kernel and an array of bandwidths \p theHs.
	 * 
	 * @param a
	 *            the left boundary of the interval
	 * @param b
	 *            the right boundary of the interval
	 * @param theHs
	 *            an array of bandwidths
	 * @param kernel
	 *            the kernel density function
	 */
	public DEKernelDensity(double a, double b, double[] theHs, ContinuousDistribution kernel) {

		this(a, b, kernel);
		this.theHs = new double[theHs.length];
		this.theHs = theHs;
	}

	/**
	 * Constructs a KDE over the interval \f$[a,b]\f$ with an array of bandwidths
	 * \p theHs.
	 * 
	 * @param a
	 *            the left boundary of the interval
	 * @param b
	 *            the right boundary of the interval
	 * @param theHs
	 *            an array of bandwidths
	 */
	public DEKernelDensity(double a, double b, double[] theHs) {

		this(a, b);
		this.theHs = new double[theHs.length];
		this.theHs = theHs;
	}

	/**
	 * Gives the kernel density function.
	 * 
	 * @return the kernel density function.
	 */
	public ContinuousDistribution getKernel() {
		return kernel;
	}

	/**
	 * Sets the kernel density function to \p kernel.
	 * 
	 * @param kernel
	 *            the kernel density function to be used.
	 */
	public void setKernel(ContinuousDistribution kernel) {
		this.kernel = kernel;
	}

	/**
	 * Gives the threshold-level for the evaluation of the density. Contributions of
	 * observations smaller than this value are considered negligible.
	 * 
	 * @return the threshold-level for evaluation.
	 */
	public double getEps() {
		return eps;
	}

	/**
	 * Sets the threshold-level for the evaluation of the density to \p eps.
	 * Contributions of observations smaller than this value are considered
	 * negligible.
	 * 
	 * @param eps
	 *            the threshold-level for evaluation.
	 */
	public void setEps(double eps) {
		this.eps = eps;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void constructDensity(double[] data) {
		dist = new EmpiricalDist(data);

	}

	/**
	 * Returns the value of the density evaluated at point \f$x\f$. For this
	 * evaluation method to work we assume that {@link #kernel} is unimodal, i.e.
	 * increasing then decreasing. This algorithm works by summing up all
	 * contributions of all observations \f$x_0, \ldots, x_{n-1}\f$ which are larger
	 * than {@link #eps}.
	 * 
	 * @param x
	 *            the point at which the density is to be evaluated.
	 * @return the value of the estimated density at \f$x\f$.
	 */
	@Override
	public double evalDensity(double x) {
		double density;
		int n = dist.getN();
		double invhn = 1.0 / (h * (double) n);
		double invh = 1.0 / h;
		double sum = 0.0;
		double term; // A term to be added to the sum that defines the density
						// estimate.
		int imin = 0; // We know that the terms for i < imin do not contribute
						// significantly.
		// Evaluation points are indexed by j.
		// int imin= pmin;
		term = kernel.density((x - dist.getObs(imin)) * invh);
		while ((term < eps) && (imin < n - 1) && (dist.getObs(imin) < x))
			term = kernel.density((x - dist.getObs(++imin)) * invh);
		// System.out.println(imin);
		// pmin=imin;
		sum = term; // The first significant term.
		for (int i = imin + 1; (i < n) && ((term > eps) || (dist.getObs(i) < x)); i++) // Data
																						// indexed
			// by i.
			sum += (term = kernel.density((x - dist.getObs(i)) * invh));
		density = sum * invhn;
		// System.out.println(density);
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
