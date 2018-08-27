package umontreal.ssj.stat.density;

/**
 * Implements an abstract class for density estimators that rely on the
 * selection of a bandwidth \f$h\f$ such as, e.g., \ref DEHistogram or \ref
 * DEKernelDensity. Furthermore, it provides basic methods such as to estimate
 * the integrated variance (IV), for instance.
 * 
 * The IV of such estimators usually locally follow the model
 * \f$\textrm{IV}\approx C n^{-\beta}h^{-\delta}\f$ and their integrated square
 * bias (ISB) is of the form \f$\textrm{ISB}\approx B h^{\alpha}\f$, where
 * \f$n\f$ denotes the number of data points used to construct the density.
 * 
 * @author puchhamf
 *
 */
public abstract class DEBandwidthBased implements DensityEstimator {

	protected double h;
	/** < bandwidth */
	protected double[] theHs;
	/** < array containing various bandwidths */

	protected double alpha, beta, delta, C, B;
	/** < model parameters for IV and ISB */

	protected double a;
	/** < left boundary of the interval over which we estimate. */
	protected double b;

	/** < right boundary of the interval over which we estimate. */

	/**
	 * {@inheritDoc}
	 */
	public void setRange(double a, double b) {
		this.a = a;
		this.b = b;
	}

	/**
	 * Gives the bandwidth.
	 * 
	 * @return the bandwidth.
	 */
	public double getH() {
		return h;
	}

	/**
	 * Sets the bandwidth to the value of \f$h\f$.
	 * 
	 * @param h
	 *            the desired bandwidth
	 */
	public void setH(double h) {
		this.h = h;
	}

	/**
	 * Gives the array {@link #theHs} containing various bandwidths.
	 * 
	 * @return array with various bandwidths.
	 */
	public double[] getTheHs() {
		return theHs;
	}

	/**
	 * Sets the optional array containing various bandwidths {@link #theHs} to the
	 * value of \f$tH\f$.
	 * 
	 * @param tH
	 *            the desired array of bandwidths
	 */
	public void setTheHs(double[] tH) {
		this.theHs = new double[tH.length];
		this.theHs = tH;

	}

	/**
	 * Gives the current value of {@link #alpha}.
	 * 
	 * @return {@link #alpha}
	 */
	public double getAlpha() {
		return alpha;
	}

	/**
	 * Sets the current value of {@link #alpha} to \f$\alpha\f$.
	 * 
	 * @param alpha
	 *            the desired value for {@link #alpha}.
	 */
	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}

	/**
	 * Gives the current value of {@link #B}.
	 * 
	 * @return {@link #B}
	 */
	public double getB() {
		return B;
	}

	/**
	 * Sets the current value of {@link #B} to \f$B\f$.
	 * 
	 * @param B
	 *            the desired value for {@link #B}
	 */
	public void setB(double B) {
		this.B = B;
	}

	/**
	 * Sets the current value of {@link #beta} to \f$\beta\f$.
	 * 
	 * @param the
	 *            desired value for {@link #beta}.
	 */
	public void setBeta(double beta) {
		this.beta = beta;
	}

	/**
	 * Gives the current value of {@link #beta}.
	 * 
	 * @return {@link #beta}
	 */
	public double getBeta() {
		return beta;
	}

	/**
	 * Sets the current value of {@link #C} to \f$C\f$.
	 * 
	 * @param the
	 *            desired value for {@link #C}.
	 */
	public void setC(double C) {
		this.C = C;
	}

	/**
	 * Gives the current value of {@link #C}.
	 * 
	 * @return {@link #C}
	 */
	public double getC() {
		return C;
	}

	/**
	 * Gives the current value of {@link #delta}.
	 * 
	 * @return {@link #delta}
	 */
	public double getDelta() {
		return delta;
	}

	/**
	 * Sets the current value of {@link #delta} to \f$\delta\f$.
	 * 
	 * @param the
	 *            desired value for {@link #delta}.
	 */
	public void setDelta(double delta) {
		this.delta = delta;
	}

	
	
	

	/**
	 * Constructs a density estimator from the data points in the vector \p data.
	 * 
	 * @param data
	 *            the data points.
	 */
	public abstract void constructDensity(double[] data);

	/**
	 * Gives the value of the density estimator evaluated at the point \f$x\f$.
	 * 
	 * @param x
	 *            evaluation point
	 * @return the density estimator at \f$x\f$.
	 */
	public abstract double evalDensity(double x);

	/**
	 * Gives a short description of the estimator.
	 * 
	 * @return a short description.
	 */
	public abstract String toString();

}
