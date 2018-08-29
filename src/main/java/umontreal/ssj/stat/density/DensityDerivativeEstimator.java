package umontreal.ssj.stat.density;

import umontreal.ssj.probdist.EmpiricalDist;

/**
 * This class implements a density derivative estimator (DDE) over a finite
 * interval \f$[a,b] \f$ based on a kernel density estimator (KDE) with kernel
 * \f$k\f$. Such an estimator is of the form \f[ \hat{f}^{(r)}_n(x) = \frac{1}{n
 * h^{r + 1}} \sum_{i = 0}^{n - 1} k^{(r)}\left( \frac{x - X_i}{h} \right), \f]
 * where \f$X_0,X_1,\dots,X_{n-1} \f$ denote \f$n\f$ observations simulated from
 * the underlying model and \f$h\f$ the bandwidth.
 * 
 * This class is particularly useful for estimating optimal bandwidths and/or
 * the constant value in the asymptotic bias for histograms and KDEs. Recall
 * that asymptotically, in both cases, \f$\textrm{ISB} \approx B h^{\alpha}\f$.
 * Now, for MC and hence for RQMC too, it is known (see \cite tSCO15a) that \f$
 * B = R(f')/12 \f$ for histograms and \f$B = \mu_2(k)R(f'')/4 \f$ for KDEs,
 * where \f[ R(g) = \int_a^b g^2(x)\mathrm{d}\!x\f] denotes the roughness
 * functional and where \f$\mu_k \f$ denotes the \f$k\f$th moment. Thus, in
 * order to estimate \f$B\f$ one needs to estimate derivatives of the sought
 * density.
 * 
 * To select the asymptotically optimal (w.r.t. the MISE) bandwidth for a DDE
 * for the \f$r\f$-th derivative using \f$n\f$ observations one can use the
 * formula given in \cite tRAY06a \f[ h_{\text{AMISE}}^{(r)} = \left[
 * \frac{\mu_2(k^{(r)}) (2r+1)}{\mu^2_2(k) R(f^{(r+2)})n} \right]^{1/(2r+5)}.
 * \f]
 * 
 * So, in order to get a good estimate for \f$f^{(r)} \f$ we would need a good
 * estimate for \f$f^{(r+2)}\f$. When estimating \f$ B\f$, one usually resorts
 * to taking a reasonably good estimate for \f$f^{(r+2)}\f$ as initial value and
 * then proceeding iteratively up to the first or second derivative. To obtain
 * such a reasonable initial value, one can, for instance, assume \f$ f \f$ to
 * be a known density \f$p\f$ whose parameters can be estimated from the
 * observations \f$X_0,X_1,\dots,X_{n-1}\f$ and for which the roughness
 * functional \f$R(p^{(r+2)})\f$ can be estimated or even computed analytically.
 * 
 * It is worth mentioning that the second moment of \f$f^{(r)}\f$ can also be
 * expressed in terms of \f$\Phi_{2r}(f)\f$, the so-called density functional of
 * even order \f$2r \f$ , i.e. \f[\mu_2(f^(r)) = (-1)^r \int
 * f^{(2r)}(x)f(x)\mathrm{d}\!x = (-1)^r \mathbb{E}[f^{(2r)}] = \Phi_{2r}(f).\f]
 * 
 * 
 * @author puchhamf
 *
 */
public abstract class DensityDerivativeEstimator extends DEBandwidthBased {

	protected int order;
	/**<order of the derivative we want to estimate */
	protected EmpiricalDist dist;
	/**<contains the observations \f$X_0,\dots,X_{n-1}\f$. */

	
	/**
	 * Gives the {@link #order} of the DDE.
	 * 
	 * @return the order of the DDE.
	 */
	public int getOrder() {
		return order;
	}

	/**
	 * Sets the {@link #order} of the DDE to \a order
	 * 
	 * @param order
	 *            the desired order.
	 */
	public void setOrder(int order) {
		this.order = order;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void constructDensity(double[] data) {
		dist = new EmpiricalDist(data);

	}

	/**
	 * Given a value \a init for the roughness functional of \f$ f^{(r+2)}\f$, \a
	 * mu2 the second moment of \f$k\f$, and \a mu2Derivative the second moment of
	 * \f$k^{(r)}\f$ this method computes the asymptotically optimal bandwidth for
	 * the DDE based on @f$n@f$ observations for the @f$r@f$-th derivative of the
	 * sought density @f$f@f$.
	 * 
	 * @param r
	 *            the order of the sought derivative.
	 * @param mu2
	 *            the second moment of \f$k\f$.
	 * @param mu2Derivative
	 *            the second moment of \f$k^{(r)}\f$.
	 * @param init
	 *            estimate of the roughness functional of the @f$(r+2)@f$-th
	 *            derivative of the density
	 * @param n
	 *            the number of observations.
	 * @return the asymptotically optimal bandwidth.
	 */
	public static double hAmise(int r, double mu2, double mu2Derivative, double init, int n) {
		double invN = 1.0 / (double) n;
		double sign = r % 2 == 0 ? 1.0 : -1.0;
		double h = sign * mu2Derivative * (2.0 * (double) r + 1.0) * invN / (mu2 * init);
		// System.out.println("densFuncGaus = " + h);
		double exp = 1.0 / (2.0 * (double) r + 5.0);
		return Math.pow(h, exp);
	}

	/**
	 * Depending on whether \a maxDerivative is even or odd, this method computes
	 * \f$ h_{\text{AMISE}}^{(r_0)}\f$  for \f$r_0=2\f$ or \f$r_0=1\f$, respectively.
	 * This is done by recursively applying the formula from above. The recursion
	 * starts at @f$r=\text{maxDerivative} - 2@f$ with the initial value \a init
	 * for @f$R(f^{(r+2)})@f$. The second moment of \f$k\f$ is given by \a mu2 and the
	 * second moments of the derivatives \f$k^{(r)},k^{(r-2)},\dots,k^{(r_0)} \f$ are passed via
	 * \a mu2Derivatives. It is assumed that the elements of \a mu2Derivatives are sorted
	 * by order of derivative in a decreasing manner, the first one being \f$\mu_2(k^{(r)})\f$.
	 * 
	 * @param dde the DDE to use in the recursion.
	 * @param data
	 *            observations from one simulation.
	 * @param mu2
	 * 			the second moment of \f$k\f$.
	 * @param mu2Derivatives
	 * 			the second moments of the derivatives of \f$k\f$.
	 * @param evalPoints
	 *            the integration nodes.
	 * @param maxDerivative
	 *            the highest order of the density derivative considered.
	 * @param init
	 *            initial value of the roughness functional in the denominator.
	 * @return either \f$ h_{\text{AMISE}}^{(2)}\f$ or \f$
	 *         h_{\text{AMISE}}^{(1)}\f$, depending on whether \a maxDerivative is
	 *         even or odd.
	 */
	public static double hOptAmise(DensityDerivativeEstimator dde, double[] data, double mu2, double[] mu2Derivatives, double[] evalPoints,
			int maxDerivative, double init ) {

		int n = data.length;
		int order = maxDerivative - 2;
		//Index of the right derivative in mu2Derivatives
		int derIndex = 0;
		

		int numEvalPoints = evalPoints.length;

		double[] estDensDerivative = new double[numEvalPoints];
		// Arrays.fill(estDensDerivative, 0.0);

		double h = hAmise(order,mu2,mu2Derivatives[derIndex],init, n);
		// System.out.println("hAmiseGaussian = " + h);
		order -= 2;
		derIndex++;
		double roughnessFunctional = init;
		// only enters loop when the last computed h is not already h^(2) or
		// h^(1)
		while (order > 0) {

			dde.setOrder(order+2);
			dde.setH(h);
			dde.constructDensity(data);
			dde.evalDensity(evalPoints, estDensDerivative);
			roughnessFunctional = 0.0;
			for (int i = 0; i < numEvalPoints; i++)
				roughnessFunctional += estDensDerivative[i] * estDensDerivative[i];
			roughnessFunctional *= (dde.getMax() - dde.getMin()) / (double) numEvalPoints;

			h = hAmise(order, mu2,mu2Derivatives[derIndex],roughnessFunctional,n);
			order -= 2;
			derIndex++;

		}

		return h;
	}
	
	/**
	 * Same as {@link #hOptAmise(DensityDerivativeEstimator, double[], double, double[], double[], int, double)}
	 * but with \a numEvalPoints equidistant points as integration nodes.
	 * @param dde the DDE to use in the recursion.
	 * @param data
	 *            observations from one simulation.
	 * @param mu2
	 * 			the second moment of \f$k\f$.
	 * @param mu2Derivatives
	 * 			the second moments of the derivatives of \f$k\f$.
	 * @param numEvalPoints
	 *            the number of equidistant integration nodes.
	 * @param maxDerivative
	 *            the highest order of the density derivative considered.
	 * @param init
	 *            initial value of the roughness functional in the denominator.
	 * @return either \f$ h_{\text{AMISE}}^{(2)}\f$ or \f$
	 *         h_{\text{AMISE}}^{(1)}\f$, depending on whether \a maxDerivative is
	 *         even or odd.
	 * @return
	 */
	public static double hOptAmise(DensityDerivativeEstimator dde, double[] data, double mu2, double[] mu2Derivatives, int numEvalPoints,
			int maxDerivative, double init ) {
		return hOptAmise(dde,data,mu2,mu2Derivatives,dde.equidistantPoints(numEvalPoints),maxDerivative,init); 
	}
	
	/**
	 * Computes the asymptotically optimal multiplicative constant in the ISB for a
	 * histogram as described above. The bandwidth in these density derivate
	 * estimates is computed recursively via
	 * {@link #hOptAmise(DensityDerivativeEstimator, double[], double, double[], double[], int, double)}
	 * starting with initial value \a init for the roughness functional of the
	 * derivative of the density of order \a maxDerivative. Each occuring
	 * integration, e.g. for evaluating a roughness-functional, is carried out by a
	 * quadrature rule using the integration nodes \a evalPoints.
	 * 
	 * Analogously to
	 * {@link #hOptAmise(DensityDerivativeEstimator, double[], double, double[], double[], int, double)}
	 * \a mu2 is the second moment of the kernel and \a mu2Derivatives contains the
	 * second moments of the derivatives of \f$k\f$ decreasing w.r.t. the order of
	 * the derivatives, starting with order \a maxDerivative - 2.
	 * 
	 * Note that for a histogram estimator \a maxDerivative has to be odd!
	 * 
	 * @param dde
	 *            the DDE we use for determining the optimal bandwidth.
	 * @param data
	 *            the observations gained from @f$m@f$ independent simulations.
	 * @param mu2
	 *            the second moment of \f$k\f$.
	 * @param mu2Derivatives
	 *            the second moments of the derivatives of \f$k\f$.
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
	//TODO: add exception if maxDerivative is odd!
	public static double computeB(DensityDerivativeEstimator dde, double[][] data, double mu2, double[] mu2Derivatives, double[] evalPoints,
			int maxDerivative, double init,	DEHistogram de) {
		int m = data.length;
		double B = 0.0;
		double a = de.getMin();
		double b = de.getMax();

		int numEvalPoints = evalPoints.length;
		double[] estDensDerivative = new double[numEvalPoints];
//		Arrays.fill(estDensDerivative, 0.0);

		for (int r = 0; r < m; r++) {
			double h = hOptAmise(dde, data[r],mu2,mu2Derivatives, evalPoints, maxDerivative, init);
			// System.out.println("h = " + h);
			dde.setOrder(1);
			dde.setH(h);
			dde.constructDensity(data[r]);
			dde.evalDensity(evalPoints, estDensDerivative);
			double roughnessFunctional = 0.0;
			for (int i = 0; i < numEvalPoints; i++)
				roughnessFunctional += estDensDerivative[i] * estDensDerivative[i];
			roughnessFunctional *= (b - a) / (double) numEvalPoints;

			B += roughnessFunctional;
		}

		
			return B / (12.0 * (double) m);

		

	}
	/**
	 * Same as {@link #computeB(DensityDerivativeEstimator, double[][], double, double[], double[], int, double, DEHistogram)}
	 * but with \a numEvalPoints equidistant integration nodes.
	 * @param dde the DDE we use for determining the optimal bandwidth.
	 * @param data the observations gained from @f$m@f$ independent simulations. 
	 * @param mu2
	 * 			the second moment of \f$k\f$.
	 * @param mu2Derivatives
	 * 			the second moments of the derivatives of \f$k\f$.
	 * @param numEvalPoints
	 *            the number of equidistant integration nodes.
	 * @param maxDerivative
	 *            the highest order of the density derivative considered.
	 * @param init
	 *            initial value of the roughness functional in the denominator.
	 * @param de the histogram estimator for the sought density. 
	 * @return the asymptotically optimal multiplicative constant in the ISB. 
	 */

	public static double computeB(DensityDerivativeEstimator dde, double[][] data, double mu2, double[] mu2Derivatives, int numEvalPoints,
			int maxDerivative, double init,	DEHistogram de) {
		return computeB(dde,data,mu2,mu2Derivatives,dde.equidistantPoints(numEvalPoints),maxDerivative,init,de);
	}

	/**
	 * Computes the asymptotically optimal multiplicative constant in the ISB for a KDE as described above.
	 *  The bandwidth in these density derivate estimates is computed recursively via 
	 *  {@link #hOptAmise(DensityDerivativeEstimator, double[], double, double[], double[], int, double)}
	 * starting  with initial value \a init  for the roughness functional of the  derivative of the density
	 * of order \a maxDerivative.
	 * Each occuring integration, e.g. for evaluating a roughness-functional,
	 * is carried out by a quadrature rule using the integration nodes \a evalPoints.
	 * 
	 * Analogously to {@link #hOptAmise(DensityDerivativeEstimator, double[], double, double[], double[], int, double)}
	 * \a mu2 is the second moment of the kernel and \a mu2Derivatives contains the second moments of the derivatives
	 * of \f$k\f$ decreasing w.r.t. the order of the derivatives, starting with order \a maxDerivative - 2.
	 * 
	 * Note that for a KDE \a maxDerivative has to be even!
	 * 
	 * @param dde the DDE we use for determining the optimal bandwidth.
	 * @param data the observations gained from @f$m@f$ independent simulations. 
	 * @param * @param mu2
	 * 			the second moment of \f$k\f$.
	 * @param mu2Derivatives
	 * 			the second moments of the derivatives of \f$k\f$.
	 * @param evalPoints
	 *            the integration nodes.
	 * @param maxDerivative
	 *            the highest order of the density derivative considered.
	 * @param init
	 *            initial value of the roughness functional in the denominator.
	 * @param de the KDE for the sought density. 
	 * @return the asymptotically optimal multiplicative constant in the ISB. 
	 */
	//TODO: add exception that maxDerivative has to be even!
	public static double computeB(DensityDerivativeEstimator dde, double[][] data, double mu2, double[] mu2Derivatives, double[] evalPoints,
			int maxDerivative, double init,	DEKernelDensity de) {
		int m = data.length;
		double B = 0.0;
		double a = de.getMin();
		double b = de.getMax();

		int numEvalPoints = evalPoints.length;
		double[] estDensDerivative = new double[numEvalPoints];
//		Arrays.fill(estDensDerivative, 0.0);

		for (int r = 0; r < m; r++) {
			double h = hOptAmise(dde, data[r],mu2,mu2Derivatives, evalPoints, maxDerivative, init);
			// System.out.println("h = " + h);
			dde.setOrder(2);
			dde.setH(h);
			dde.constructDensity(data[r]);
			dde.evalDensity(evalPoints, estDensDerivative);
			double roughnessFunctional = 0.0;
			for (int i = 0; i < numEvalPoints; i++)
				roughnessFunctional += estDensDerivative[i] * estDensDerivative[i];
			roughnessFunctional *= (b - a) / (double) numEvalPoints;

			B += roughnessFunctional;
		}

		
			return 0.25 * mu2 * B / (double) m;	

	}
	
	/**
	 * Same as {@link #computeB(DensityDerivativeEstimator, double[][], double, double[], double[], int, double, DEKernelDensity)}
	 * but with \a numEvalPoints equidistant integration nodes.
	 * @param dde the DDE we use for determining the optimal bandwidth.
	 * @param data the observations gained from @f$m@f$ independent simulations. 
	 * @param mu2
	 * 			the second moment of \f$k\f$.
	 * @param mu2Derivatives
	 * 			the second moments of the derivatives of \f$k\f$.
	 * @param numEvalPoints
	 *            the number of equidistant integration nodes.
	 * @param maxDerivative
	 *            the highest order of the density derivative considered.
	 * @param init
	 *            initial value of the roughness functional in the denominator.
	 * @param de the KDE for the sought density. 
	 * @return the asymptotically optimal multiplicative constant in the ISB. 
	 */

	public static double computeB(DensityDerivativeEstimator dde, double[][] data, double mu2, double[] mu2Derivatives, int numEvalPoints,
			int maxDerivative, double init,	DEKernelDensity de) {
		return computeB(dde,data,mu2,mu2Derivatives,dde.equidistantPoints(numEvalPoints),maxDerivative,init,de);
	}
}
