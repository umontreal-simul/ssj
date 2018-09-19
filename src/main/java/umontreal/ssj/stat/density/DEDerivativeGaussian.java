package umontreal.ssj.stat.density;

import umontreal.ssj.probdist.NormalDist;

/**
 * This class implements a density derivative estimator (DDE) with a Gaussian (
 * i.e., standard normal) kernel function. While many general methods are
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
 * This class also offers static methods so that the user can simply evaluate the density based on
 * a set of observations without having to construct a KDE.
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
	 * Constructs a DDE with a Gaussian kernel function of order \a order based on the observations \a data.
	 * @param order order the order of the derivative considered.
	 * @param data the observations.
	 */
	public DEDerivativeGaussian(int order, double[] data) {
		this(order);
		this.data = data;
	}
	
	/**
	 * Constructs a DDE with a Gaussian kernel function with bandwidth \a h.
	 * @param h the bandwidth.
	 */
	public DEDerivativeGaussian(double h) {
		setH( h);
	}
	
	/**
	 * Constructs a DDE with a Gaussian kernel function with bandwidth \a h based on the observations \a data.
	 * @param h the bandwidth.
	 * @param data the observations.
	 */
	public DEDerivativeGaussian(double h, double[] data) {
		this(h);
		this.data = data;
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
	 * Constructs a DDE with a Gaussian kernel function of order \a order and
	 * bandwidth \a h based on the observations \a data.
	 * @param order the order of the derivative considered.
	 * @param h the bandwidth.
	 * @param data the observations.
	 */
	public DEDerivativeGaussian(int order, double h, double [] data) {
		this(order,h);
		this.data = data;

	}

	/**
	 *Evaluates the DDE at the point \a x.
	 *@param x the evaluation point.
	 *@return the DDE evaluated at \a x.
	 * 
	 */
	@Override
	public double evalDensity(double x) {
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
	 * Evaluates the DDE at all the points in \a evalPoints and the resulting values are
	 * returned in an array.
	 * 
	 * @param evalPoints the evaluation points.
	 * @return the DDE evaluated at the points \a evalPoints.
	 */
	@Override
	public double[] evalDensity(double[] evalPoints) {
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
	
	//STATIC METHODS
	/**
	 * Evaluates the DDE of order \a order with bandwidth \a h which
	 * is defined by the observations \a data at the evaluation point \a x. 
	 * @param x the evaluation point.
	 * @param order the order of the DDE.
	 * @param h the bandwidth.
	 * @param data the observations.
	 * @return The DDE defined by the above parameters evaluated at \a x
	 */
	public static double evalDensity(double x, int order, double h, double [] data) {
		double z;
		double sign = order % 2 == 0 ? 1.0 : -1.0;
		double density = 0.0;
		int n = data.length;
		double nInv = 1.0 / (double) n;
		double hInv = 1.0 / h;

		double norma = sign * nInv * Math.pow(hInv, order + 1.0);

		for (int i = 0; i < n; i++) {
			z = (x - data[i]) * hInv;

			density += NormalDist.density01(z) * hermitePoly(order, z);
		}
		density *= norma;
		return density;
	}
	
	/**
	 * Evaluates the DDE of order \a order with bandwidth \a h which
	 * is defined by the observations \a data at each of the evaluation points in \a evalPoints. 
	 * @param evalPoints the evaluation points.
	 * @param order the order of the DDE.
	 * @param h the bandwidth.
	 * @param data the observations.
	 * @return The DDE defined by the above parameters evaluated at each point in \a evalPoint.
	 */
	public static double [] evalDensity(double[] evalPoints , int order, double h, double [] data) {
		double z = 0.0;
		double sign = order == 0 ? 1.0 : -1.0;
		int k = evalPoints.length;
		double[] density = new double[k];
		int n = data.length;
		double nInv = 1.0 / (double) n;
		double hInv = 1.0 / h;

		double norma = sign * nInv * Math.pow(hInv, order + 1.0);

		for (int j = 0; j < k; j++) { // evalPoints indexed by j
			for (int i = 0; i < n; i++) { // data points indexed by i
				z = (evalPoints[j] - data[i]) * hInv;
			}
			density[k] += NormalDist.density01(z) * hermitePoly(order, z);
			density[k] *= norma;
		}

		return density;
	}
	
	/**
	 * Assume that we have \f$m\f$ independent realizations of the underlying model.
	 * For each such realization this method evaluates a DDE of order \a order with 
	 * bandwidth \a h at the points from \a evalPoints. The
	 * independent realizations are passed via the 2-dimensional \f$m\times n\f$
	 * array \a data, where \f$n\f$ denotes the number of observations per
	 * realization. Hence, its first index identifies the independent realization
	 * while its second index identifies a specific observation of this realization.
	 * 
	 * The result is returned as a \f$m\times k\f$ matrix, where \f$k \f$ is the
	 * number of evaluation points, i.e., the length of \a evalPoints. The first
	 * index, again, identifies the independent realization whereas the second index
	 * corresponds to the point of \a evalPoints at which the DDE was
	 * evaluated.
	 * @param evalPoints the evaluation points.
	 * @param order the order of the DDE.
	 * @param h the bandwidth.
	 * @param data the observations.
	 * @return the DDE for each realization evaluated at \a evalPoints.
	 */
	public static double [][] evalDensity(double[] evalPoints , int order, double h, double [][] data){
		int m = data.length;
		double[][] density= new double[m][];
		for(int r = 0; r < m; r++)
			density[r] = evalDensity(evalPoints,order,h,data[r]);
		return density;
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
