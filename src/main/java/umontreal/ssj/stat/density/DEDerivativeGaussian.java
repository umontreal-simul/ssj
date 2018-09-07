package umontreal.ssj.stat.density;

import umontreal.ssj.probdist.NormalDist;

/**
 * This class implements a density derivative estimator (DDE) with a Gaussian (
 * i.e., standard normal) Kernel function. While the many general methods are
 * already handled by its superclass, this class is mainly concerned with the
 * construction and evaluation of this kind of DDE. The formula specific to a
 * Gaussian DDE is given by
 * 
 * \f[ \hat{f}^{(r)}_{n,h}(x) = \frac{(-1)^r}{n h^{r+1}}\sum_{i=0}^{n-1} \phi
 * \left( \frac{x-X_i}{h}\right) H_r\left( \frac{x-X_i}{h}\right), \f]
 * 
 * where \f$H_r\f$ denotes the probabilist's Hermite polynomial of order \f$r\f$ and 
 * \f$\phi\f$ denotes the standard normal density.
 * 
 * @author puchhamf
 *
 */

public class DEDerivativeGaussian extends DensityDerivativeEstimator {

	/**
	 * Constructs a DDE with a Gaussian kernel function of order \a order.
	 * 
	 * @param order the order of the derivative considered.
	 */
	public DEDerivativeGaussian(int order) {
		setOrder(order);
	}

	/**
	 * Constructs a DDE with a Gaussian kernel function of order \a order and
	 * bandwidth \a h.
	 * 
	 * @param h     the bandwidth.
	 * @param order the order of the derivative considered.
	 */
	public DEDerivativeGaussian(int order, double h) {

		this(order);
		setH(h);

	}

	/**
	 * {@inheritDoc}
	 * 
	 */
	@Override
	public double evalDensity(double x, double[] data, double a, double b) {
		double z;
		double sign = getOrder() % 2 == 0 ? 1.0 : -1.0;
		double density = 0.0;
		int n = data.length;
		double nInv = 1.0 / (double) n;
		double hInv = 1.0 / getH();

		double norma = sign * nInv * Math.pow(hInv, getOrder() + 1.0);

		for (int i = 0; i < n; i++) {
			z = (x - data[i]) * hInv;

			density += NormalDist.density01(z) * hermitePoly(getOrder(), z);
		}
		density *= norma;
		return density;
	}

	/**
	 * Constructs a DDE with a Gaussian kernel function from the observations \a
	 * data and the bandwith \f$h\f$ set to #getH(). Furthermore, the DDE is
	 * evaluated at the points in \a evalPoints and the resulting values are
	 * returned in an array.
	 * 
	 * @param evalPoints the evaluation points.
	 * @param data       the observations for constructing the density estimator.
	 * @param a          the left boundary of the interval.
	 * @param b          the right boundary of the interval
	 * @return the density estimator evaluated at the points \a evalPoints.
	 */
	@Override
	public double[] evalDensity(double[] evalPoints, double[] data, double a, double b) {
		double z = 0.0;
		double sign = getOrder() % 2 == 0 ? 1.0 : -1.0;
		int k = evalPoints.length;
		double[] density = new double[k];
		int n = data.length;
		double nInv = 1.0 / (double) n;
		double hInv = 1.0 / getH();

		double norma = sign * nInv * Math.pow(hInv, getOrder() + 1.0);

		for (int j = 0; j < k; j++) { // evalPoints indexed by j
			for (int i = 0; i < n; i++) { // data points indexed by i
				z = (evalPoints[j] - data[i]) * hInv;
			}
			density[k] += NormalDist.density01(z) * hermitePoly(getOrder(), z);
			density[k] *= norma;
		}

		return density;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String toString() {
		return "DDE [Gaussian kernel with h = " + getH() + "]";
	}

	/**
	 * Computes the probabilist's Hermite polynomial of order \a r at \a x, which is
	 * defined by the recursion \f[H_{r+1}(x)= x H_r(x) - r H_{r-1}(x) \f] with
	 * initial values \f$H_0(x) = 0\f$, \f$H_1(x) = x\f$.
	 * 
	 * @param r the order of the Hermite polynomial.
	 * @param x the evaluation point.
	 * @return the probabilist's Hermite polynomial.
	 * 
	 * @remark **Florian:** This should probably be located elsewhere.
	 */
	public static double hermitePoly(int r, double x) {
		if (r == 0)
			return 1.0;
		else if (r == 1)
			return x;
		else
			return hermitePoly(r - 1, x) * x - ((double) r - 1.0) * hermitePoly(r - 2, x);
	}

}
