package umontreal.ssj.stat.density;

import umontreal.ssj.probdist.EmpiricalDist;

/**
 * This class implements a density derivative estimator (DDE) based on a kernel
 * density estimator (KDE) with kernel \f$k\f$. Such an estimator is of the form
 * \f[ \hat{f}^{(r)}_n(x) = \frac{1}{n h^{r + 1}} \sum_{i = 0}^{n - 1}
 * k^{(r)}\left( \frac{x - X_i}{h} \right), \f] where \f$X_0,X_1,\dots,X_{n-1}
 * \f$ denote \f$n\f$ observations simulated from the underlying model and
 * \f$h\f$ the bandwidth.
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

}
