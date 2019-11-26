package umontreal.ssj.stat.density;

/**
 * This class implements a density derivative estimator (DDE) based on a kernel
 * density estimator (KDE) with a sufficiently smooth kernel function \f$k\f$, 
 * see \ref umontreal.ssj.stat.density.DEKernelDensity. Such an estimator is used to find
 * the \f$r\f$-th derivative of an unknown density based on \f$n\f$ observations
 * \f$x_0,x_1,\dots,x_{n-1} \f$ of an underlying model. It is of the form
 * 
 * \f[ \hat{f}^{(r)}_n(x) = \hat{f}^{(r)}_{n,h}(x)=\frac{1}{n h^{r + 1}} \sum_{i
 * = 0}^{n - 1} k^{(r)}\left( \frac{x - x_i}{h} \right), \f]
 * 
 * where \f$h\f$ denotes the bandwidth. So, as a matter of fact, it is the
 * \f$r\f$-th derivative of a KDE. This class provides basic tools to construct
 * and evaluate DDEs. Note that sufficient smoothness of
 * the kernel function \f$K\f$ as well as of the unknown density \f$f\f$ is
 * required.
 * 
 * Moreover, for observations gained from a Monte Carlo simulation, the
 * asymptotically optimal (w.r.t. the mean integrated square error) bandwidth
 * \f$ h_{\text{AMISE}}^{(r)} \f$ can be computed explicitely via the formula
 * (see \cite tRAY06a)
 * 
 * @anchor REF_stat_density_DensityDerivativeEstimator_hopt
 * 
 * \f[ h_{\text{AMISE}}^{(r)} = \left[ \frac{\mu_2(k^{(r)}) (2r+1)}{\mu^2_2(k)
 * R(f^{(r+2)})n} \right]^{1/(2r+5)},\tag{hopt} \f]
 * 
 * So, in order to obtain a good estimate for \f$f^{(r)} \f$ one requires a good
 * estimate for \f$f^{(r+2)}\f$. A way to break out of this cyclic argument is
 * to resort to taking a reasonably good estimate for \f$f^{(r+2)}\f$ as initial
 * value. This is implemented as `hAmiseR(int, double, double, double, int)`. To
 * obtain such a reasonable initial value, one can, for instance, assume that
 * \f$ f \f$ belongs to a known family of distributions (e.g. normal
 * distributions), such that its defining parameters (e.g. mean, standard
 * deviation) can be estimated from the observations \f$x_0,x_1,\dots,x_{n-1}\f$
 * and for which the roughness functional \f$R(p^{(r+2)})\f$ (see
 * umontreal.ssj.stat.density.DensityEstimator#roughnessFunctional(double[],
 * double, double) ) can be easily estimated or even computed analytically.
 * 
 * Certainly, one can also try obtain a good approximation of \f$f^{(r+2t)} \f$
 * and iterate ({@link REF_stat_density_DensityDerivativeEstimator_hopt hopt})
 * \f$t\f$ times. The function `hAmiseR(int, int, double, double[], double,
 * DensityDerivativeEstimator, double[], double[], double, double)` implements
 * this recursion.
 * 
 * Since the above methods to compute \f$h_{\text{AMISE}}^{(r)}\f$ rely on initial
 * values for the second moments of derivatives of the unknown density, it is 
 * probably worth mentioning that the second moment of \f$f^{(r)}\f$ can also be
 * expressed in terms of \f$\Phi_{2r}(f)\f$, the so-called density functional of
 * even order \f$2r \f$ , i.e.
 * 
 * \f[\mu_2(f^{(r)}) = (-1)^r \int_{-\infty}^{\infty} f^{(2r)}(x)f(x)\mathrm{d}x = (-1)^r
 * \mathbb{E}[f^{(2r)}] = \Phi_{2r}(f). \f]
 * 
 * For such an initial value, one could, for instance, assume that the target
 * density is a normal distribution with standard deviation \f$\sigma\f$.
 * To this end, we include  the method `densityFunctionalGaussian(int, double)`.
 * 
 *
 *
 */
public abstract class DensityDerivativeEstimator extends DensityEstimator {

	/** order of the derivative we want to estimate. */
	private int order;

	/** the bandwith \f$h\f$ */
	private double h;

	/**
	 * Gives the order \f$r\f$ of the DDE.
	 * 
	 * @return the order of the DDE.
	 */
	public int getOrder() {
		return order;
	}

	/**
	 * Sets the order \f$r\f$ of the DDE to \a order
	 * 
	 * @param order the desired order.
	 */
	public void setOrder(int order) {
		this.order = order;
	}

	/**
	 * Gives the current bandwidth \f$h\f$.
	 * 
	 * @return the bandwidth.
	 */
	public double getH() {
		return h;
	}

	/**
	 * Sets the bandwidth \f$h\f$ to \a h
	 * 
	 * @param h the desired bandwidth.
	 */
	public void setH(double h) {
		this.h = h;
	}

	/**
	 * {@inheritDoc}
	 */
	public void setData(double[] data) {
		this.data = data;
	}
	/**
	 * Given a value \a init for the roughness functional of \f$ f^{(r+2)}\f$, \a
	 * mu2 the second moment of the kernel function \f$K\f$, and \a mu2Derivative
	 * the second moment of \f$K^{(r)}\f$, this method computes the asymptotically
	 * optimal bandwidth for the DDE based on @f$n@f$ observations simulated by
	 * Monte Carlo for the @f$r@f$-th derivative of the sought density @f$f@f$.
	 * 
	 * @param r             the order of the sought derivative.
	 * @param mu2           the second moment of \f$K\f$.
	 * @param mu2Derivative the second moment of \f$K^{(r)}\f$.
	 * @param init          estimate of the roughness functional of
	 *                      the @f$(r+2)@f$-th derivative of the density
	 * @param n             the number of observations.
	 * @return the asymptotically optimal bandwidth for the DDE of the \f$r\f$-th
	 *         derivative of the unknown density.
	 */
	public static double hAmiseR(int r, double mu2, double mu2Derivative, double init, int n) {
		double invN = 1.0 / (double) n;
		double sign = r % 2 == 0 ? 1.0 : -1.0;
		double h = sign * mu2Derivative * (2.0 * (double) r + 1.0) * invN / (mu2 * init);
		double exp = 1.0 / (2.0 * (double) r + 5.0);
		return Math.pow(h, exp);
	}

	/**
	 * Given an estimate of \f$R(f^{(r+2t)})\f$ via \a init as initial value, this
	 * function iterates over
	 * ({@link REF_stat_density_DensityDerivativeEstimator_hopt hopt}) \f$t\f$ times
	 * to obtain the asymptotically optimal bandwidth for the DDE based on @f$n@f$
	 * observations simulated by Monte Carlo for the @f$r@f$-th derivative of the
	 * sought density @f$f@f$.
	 * 
	 * Each recursion step calls #hAmiseR(int, double, double, double, int). The
	 * thereby obtained bandwidth can be used for estimating the density in the next
	 * step. The new initial value for the function call is computed by estimating
	 * the roughness functional of the corresponding derivative by \ref
	 * umontreal.ssj.stat.density.DensityEstimator#roughnessFunctional(double[],
	 * double, double) over \f$[a,b]\f$ with the quadrature points \a evalPoints.
	 * 
	 * The second moments of the derivatives of the ernel function
	 * \f$\mu_2\left(K^{(r+2(t-1))}\right),
	 * \mu_2\left(K^{r+2(t-2)}\right),\dots,\mu_2\left(K^{(r)}\right)\f$ are passed
	 * in the array \a mu2Derivative of length \f$t\f$ in exactly this order.
	 * 
	 * Note that the kernel function and the density have to be at least
	 * \f$(r+2(t-1))\f$ times and \f$(r+2t)\f$ times differentiable, respectively,
	 * to make use of this method.
	 * 
	 * @param r             the order of the sought derivative.
	 * @param t             the number of iteration steps.
	 * @param mu2           the second moment of the kernel function.
	 * @param mu2Derivative the second moments of the derivatives of the kernel
	 *                      function.
	 * @param init          the estimate of the roughness functional of the
	 *                      \f$(r+2t)\f$-th derivative of the sought density.
	 * @param dde           the DDE used for estimating the bandwidth.
	 * @param evalPoints    the quadrature points used for estimating the roughness
	 *                      functionals
	 * @param a             the left boundary of the interval considered.
	 * @param b             the right boundary of the interval considered.
	 * @return the asymptotically optimal bandwidth for the DDE of the \f$r\f$-th
	 *         derivative of the unknown.
	 */

	public static double hAmiseR(int r, int t, double mu2, double[] mu2Derivative, double init,
			DensityDerivativeEstimator dde, double[] evalPoints, double a, double b) {
		double h;
		int k = evalPoints.length;
		double[] estDensity = new double[k];
		int n = dde.data.length;

		for (int tau = t - 1; tau >= 1; tau--) {
			h = hAmiseR(r + 2 * tau, mu2, mu2Derivative[tau], init, n);
			dde.setH(h);
			dde.setOrder(r + 2 * tau);
			estDensity = dde.evalDensity(evalPoints);
			init = roughnessFunctional(estDensity, a, b);
		}

		return hAmiseR(r, mu2, mu2Derivative[0], init, n);
	}
	
	/**
	 * Computes \f$\Phi_{2r}(p)\f$, i.e. the density functional of order @f$2r@f$ of
	 * a normal density \f$p\f$ with standard deviation \a sigma,
	 * 
	 * @f[\Phi_{2r}(p) = \int_{-\infty}^{\infty} p^{(2r)}(x) p(x)\mathrm{d}x@f]
	 * 
	 * and it can be comuted as
	 * 
	 * \f[\Phi_{2r}(p) = \frac{(-1)^r (2r)!}{(2\sigma)^{2r+1} r! \sqrt{\pi}}. \f]
	 * 
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
	 * {@inheritDoc}
	 */
	public String toString() {
		return "DDE [h = " + h + "]";
	}
}
