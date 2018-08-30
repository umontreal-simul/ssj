package umontreal.ssj.stat.density;

import umontreal.ssj.probdist.NormalDist;

/**
 * This class implements a density derivative estimator (DDE) based on a KDE
 * with a Gaussian Kernel.
 * 
 * @author puchhamf
 *
 */

public class DDEGaussian extends DensityDerivativeEstimator {

	/**
	 * Constructs a DDE of order \a order over the interval \f$[a,b]\f$.
	 * 
	 * @param a
	 *            the left boundary of the interval.
	 * @param b
	 *            the right boundary of the interval.
	 * @param order
	 *            the order of the derivative considered.
	 */
	public DDEGaussian(double a, double b, int order) {
		setRange(a, b);
		this.order = order;
	}

	/**
	 * Constructs a DDE of order \a order over the interval \f$[a,b]\f$ with
	 * bandwidth \a h.
	 * 
	 * @param a
	 *            the left boundary of the interval.
	 * @param b
	 *            the right boundary of the interval.
	 * @param h
	 *            the bandwidth.
	 * @param order
	 *            the order of the derivative considered.
	 */
	public DDEGaussian(double a, double b, int order, double h) {

		this(a, b, order);
		setH(h);

	}

	/**
	 * Constructs a DDE of order \a order over the interval \f$[a,b]\f$ with an
	 * array of bandwidths \a theHs.
	 * 
	 * @param a
	 *            the left boundary of the interval
	 * @param b
	 *            the right boundary of the interval
	 * @param order
	 *            the order of the derivative.
	 * @param theHs
	 *            an array of bandwidths
	 */
	public DDEGaussian(double a, double b, int order, double[] theHs) {

		this(a, b, order);
		this.theHs = new double[theHs.length];
		this.theHs = theHs;
	}

	/**
	 * {@inheritDoc}
	 */
	public double evalDensity(double x) {
		double z;
		double sign = order % 2 == 0 ? 1.0 : -1.0;
		double density = 0.0;
		int n = dist.getN();
		double nInv = 1.0 / (double) n;
		double hInv = 1.0 / h;

		double norma = sign * nInv * Math.pow(hInv, order + 1.0);

		for (int i = 0; i < n; i++) {
			z = (x - dist.getObs(i)) * hInv;

			density += NormalDist.density01(z) * hermitePoly(order, z);
		}
		density *= norma;
		return density;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String toString() {
		return "DDE_Gaussian with h = " + h;
	}

	/**
	 * Computes \f$\Phi_{2r}(p)\f$, i.e. the density functional of order @f$2r@f$ of
	 * a normal density \f$p\f$ with standard deviation \a sigma,
	 * 
	 * @f[\Phi_{2r}(p) = \int_{\mathbb{R}} p^{(2r)}(x) p(x)\mathrm{d}\!x@f]
	 * 
	 * and it can be comuted as
	 * 
	 * \f[\Phi_{2r}(p) = \frac{(-1)^r (2r)!}{(2\sigma)^{2r+1} r! \sqrt{\pi}}. \f]
	 * 
	 * It is worth mentioning that the second moment of the \f$r\f$-th derivative of
	 * any density \f$g\f$ is related to \f$\Phi_{2r}(g)\f$ via \f$\mu_2(f^{(r)}) =
	 * (-1)^r \Phi_{2r}(f)\f$.
	 * 
	 * @param r
	 *            the order of the functional. Has to be even.
	 * @param sigma
	 *            the standard deviation of the normal density considered.
	 * @return the density functional of order @f$2r@f$.
	 */
	// TODO: look for formulas with odd r.
	public static double densityFunctionalGaussian(int r, double sigma) {
		double sign = (r % 2 == 0 ? 1.0 : -1.0);
		double facTerm = 1.0;
		for (int i = r + 1; i <= 2 * r; i++) {
			facTerm *= (double) i;
		}
		double denom = Math.pow(2.0 * sigma, (double) (2.0 * r + 1.0));

		denom *= Math.sqrt(Math.PI);

		return sign * facTerm / denom;
	}

	/**
	 * Computes the asymptotically optimal multiplicative constant in the ISB for a
	 * histogram. This works exactly as
	 * {@link #computeB(DensityDerivativeEstimator, double[][], double, double[], double[], int, double, DEHistogram)}
	 * but with the kernel of the DDE being a Gaussian.
	 * 
	 * In this case, the second moment of the kernel is 1 and those of the
	 * derivatives can be obtained via {@link #densFunctionalGaussian(int, double)}.
	 * 
	 * Note that for a histogram estimator \a maxDerivative has to be odd!
	 * 
	 * @param dde
	 *            the DDE we use for determining the optimal bandwidth.
	 * @param data
	 *            the observations gained from @f$m@f$ independent simulations.
	 * @param evalPoints
	 *            the integration nodes.
	 * @param maxDerivative
	 *            the highest order of the density derivative considered.
	 * @param init
	 *            initial value of the roughness functional in the denominator.
	 * @param de
	 *            the histogram estimator for the sought density.
	 * @return the asymptotically optimal multiplicative constant in the ISB.
	 */
	public static double computeB(DDEGaussian dde, double[][] data, double[] evalPoints, int maxDerivative, double init,
			DEHistogram de) {
		int t = (maxDerivative - 1) / 2; // maxDerivative = 2t + 1
		double[] mu2Derivatives = new double[t - 1];
		for (int r = t - 1; r >= 1; r--)
			// the first order is 2t-1, this is init.
			// the first in mu2Derivatives is 2t-3, the second 2t-5,..., the last 2*1-1.
			mu2Derivatives[t - 1 - r] = densityFunctionalGaussian(2 * r - 1, 1.0);
		return computeB(dde, data, 1.0, mu2Derivatives, evalPoints, maxDerivative, init, de);

	}

	/**
	 * Same as
	 * {@link #computeB(DDEGaussian, double[][], double[], int, double, DEHistogram)}
	 * but with \a numEvalPoints equidistant integration nodes.
	 * 
	 * @param dde
	 *            the DDE we use for determining the optimal bandwidth.
	 * @param data
	 *            the observations gained from @f$m@f$ independent simulations.
	 * @param numEvalPoints
	 *            the number of equidistant integration nodes.
	 * @param maxDerivative
	 *            the highest order of the density derivative considered.
	 * @param init
	 *            initial value of the roughness functional in the denominator.
	 * @param de
	 *            the histogram estimator for the sought density.
	 * @return the asymptotically optimal multiplicative constant in the ISB.
	 */
	public static double computeB(DDEGaussian dde, double[][] data, int numEvalPoints, int maxDerivative, double init,
			DEHistogram de) {
		return computeB(dde, data, dde.getEquidistantPoints(numEvalPoints), maxDerivative, init, de);
	}

	/**
	 * Computes the asymptotically optimal multiplicative constant in the ISB for a
	 * KDE. This works exactly as
	 * {@link #computeB(DensityDerivativeEstimator, double[][], double, double[], double[], int, double, DEKernelDensity)}
	 * but with the kernel of the DDE being a Gaussian.
	 * 
	 * In this case, the second moment of the kernel is 1 and those of the
	 * derivatives can be obtained via {@link #densFunctionalGaussian(int, double)}.
	 * 
	 * Note that for a histogram estimator \a maxDerivative has to be odd!
	 * 
	 * @param dde
	 *            the DDE we use for determining the optimal bandwidth.
	 * @param data
	 *            the observations gained from @f$m@f$ independent simulations.
	 * @param evalPoints
	 *            the integration nodes.
	 * @param maxDerivative
	 *            the highest order of the density derivative considered.
	 * @param init
	 *            initial value of the roughness functional in the denominator.
	 * @param de
	 *            the KDE for the sought density.
	 * @return the asymptotically optimal multiplicative constant in the ISB.
	 */
	public static double computeB(DDEGaussian dde, double[][] data, double[] evalPoints, int maxDerivative, double init,
			DEKernelDensity de) {
		int t = maxDerivative / 2; // maxDerivative = 2t
		double[] mu2Derivatives = new double[t - 2];
		for (int r = t - 1; r >= 2; r--)
			// the first order is 2t-2, this is init.
			// the first in mu2Derivatives is 2t-4, the second 2t-6,..., the last 2*2-2.
			mu2Derivatives[t - 1 - r] = densityFunctionalGaussian(2 * r - 2, 1.0);
		return computeB(dde, data, 1.0, mu2Derivatives, evalPoints, maxDerivative, init, de);

	}

	/**
	 * Computes the probabilist's Hermite polynomial of order \a r at \a x, which is
	 * defined by the recursion \f[H_{r+1}(x)= x H_r(x) - r H_{r-1}(x) \f] with
	 * initial values \f$H_0(x) = 0\f$, \f$H_1(x) = x\f$.
	 * 
	 * @param r
	 *            the order of the Hermite polynomial.
	 * @param x
	 *            the evaluation point.
	 * @return the probabilist's Hermite polynomial.
	 */
	// TODO: should probably be located elsewhere!!!
	public static double hermitePoly(int r, double x) {
		if (r == 0)
			return 1.0;
		else if (r == 1)
			return x;
		else
			return hermitePoly(r - 1, x) * x - ((double) r - 1.0) * hermitePoly(r - 2, x);
	}

	/**
	 * Same as
	 * {@link #computeB(DDEGaussian, double[][], double[], int, double, DEKernelDensity)}
	 * but with \a numEvalPoints equidistant integration nodes.
	 * 
	 * @param dde
	 *            the DDE we use for determining the optimal bandwidth.
	 * @param data
	 *            the observations gained from @f$m@f$ independent simulations.
	 * @param numEvalPoints
	 *            the number of equidistant integration nodes.
	 * @param maxDerivative
	 *            the highest order of the density derivative considered.
	 * @param init
	 *            initial value of the roughness functional in the denominator.
	 * @param de
	 *            the KDE for the sought density.
	 * @return the asymptotically optimal multiplicative constant in the ISB.
	 */
	public static double computeB(DDEGaussian dde, double[][] data, int numEvalPoints, int maxDerivative, double init,
			DEKernelDensity de) {
		return computeB(dde, data, dde.getEquidistantPoints(numEvalPoints), maxDerivative, init, de);
	}
}
