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
	 * 			  the order of the derivative.
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
	
	

}
