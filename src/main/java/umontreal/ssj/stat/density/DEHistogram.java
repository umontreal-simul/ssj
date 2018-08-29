package umontreal.ssj.stat.density;

import umontreal.ssj.stat.ScaledHistogram;
import umontreal.ssj.stat.TallyHistogram;

/**
 * This class provides methods to compute a histogram density estimator for
 * univariate densities over an interval \f$[a,b]\f$ from a set of \f$n\f$
 * individual observations \f$X_0, …, X_{n-1}\f$, and returns its value at a set
 * of selected points.
 * 
 * For a fixed number of bins \f$m>0\f$ we partition the interval \f$[a,b]\f$
 * into \f$m\f$ subintervals of equal lengths \f$h\f$. Observe that
 * \f$h=(b-a)/m\f$ and, hence, Histogram estimators can also be parametrized by a
 * bandwidth \f$h>0\f$, as long as the resulting number of bins is an integer.
 * The estimator itself is defined by 
 * \f[ \hat{f}_{n}(x) = \hat{f}_{n,h}(x) =
 * \frac{n_j}{nh},\quad\text{for } x\in[a+(j-1)h, a+jh), \quad j=1,\dots,m, \f] 
 * where
 * \f$n_j\f$ denotes the number of observations that fall in this interval.
 * 
 * Note that, due to the fact that the partition into \f$m\f$ bins relies on
 * half-open intervals, the boundary \f$b\f$ is not included in any of these
 * intervals. Since the probability of an observation being exactly equal to
 * \f$b\f$ is zero, we can effectively ignore this subtlety.
 * 
 * @author puchhamf
 *
 */

public class DEHistogram extends DEBandwidthBased {

	private int m;
	/**<the number of bins */
	private ScaledHistogram histDensity;
	private TallyHistogram hist;

	/**
	 * Constructs a histogram estimator over the interval \f$[a,b]\f$.
	 * 
	 * @param a
	 *            left boundary of the interval
	 * @param b
	 *            right boundary of the interval
	 */
	public DEHistogram(double a, double b) {
		setRange(a, b);
		setAlpha(2.0);
	}

	/**
	 * Constructs a histogram estimator with bandwidth \f$h\f$ over the interval
	 * \f$[a,b]\f$. Note that the actual bandwidth used might differ from \f$h\f$,
	 * since the number of bins {@link #m} has to be an integer.
	 * 
	 * @param a
	 *            left boundary of the interval
	 * @param b
	 *            right boundary of the interval
	 * @param h
	 *            the bandwidth
	 */
	public DEHistogram(double a, double b, double h) {
		this(a, b);
		setH(h);
	}

	/**
	 * Constructs a histogram estimator with @f$m@f$ bins over the
	 * interval @f$[a,b]@f$.
	 * 
	 * @param a
	 *            left boundary of the interval
	 * @param b
	 *            right boundary of the interval
	 * @param m
	 *            the desired number of bins
	 */
	public DEHistogram(double a, double b, int m) {
		this(a, b);
		setM(m);
	}

	/**
	 * 
	 * @param a
	 *            left boundary of the interval
	 * @param b
	 *            right boundary of the interval
	 * @param theHs
	 *            an array of bandwidths
	 * 
	 */
	public DEHistogram(double a, double b, double[] theHs) {
		this(a, b);
		this.theHs = new double[theHs.length];
		this.theHs = theHs;
	}

	/**
	 * Gives the number of bins.
	 * 
	 * @return the number of bins.
	 */
	public int getM() {
		return m;
	}

	/**
	 * Sets the number of bins to \f$m\f$.
	 * 
	 * @param m
	 *            the desired number of bins
	 */
	public void setM(int m) {
		this.m = m;
		setH((getMax() - getMin()) / (double) m);
	}

	/**
	 * Sets the bandwidth to \f$h\f$. Note that the actual bandwidth used might
	 * differ from \f$h\f$, since the number of bins {@link #m} has to be an
	 * integer.
	 */
	@Override
	public void setH(double h) {
		this.m = (int) ((getMax() - getMin()) / h);
		this.h = (getMax() - getMin()) / (double) m;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void constructDensity(double[] data) {
		hist = new TallyHistogram(getMin(), getMax(), m);
		hist.fillFromArray(data);
		histDensity = new ScaledHistogram(hist, 1.0);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public double evalDensity(double x) {
		double h = hist.getH();
		return histDensity.getHeights()[(int) ((x - getMin()) / h)];
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String toString() {
		return "Histogram estimator with " + m + " bins.";
	}

	/**
	 * Overrides the method in \ref DensityEstimator so that the number of points is
	 * reset to the minimum of \a numPoints and {@link #m}.
	 * 
	 * @param numPoints
	 *            the number of points to be returned.
	 * @return an array of equidistant points over \f$[a,b]\f$.
	 */
	// TODO: does this work???
	@Override
	protected double[] equidistantPoints(int numPoints) {
		int trueNumPoints = Math.max(numPoints, m);
		double evalPoints[] = new double[trueNumPoints];
		double delta = (b - a) / (trueNumPoints);
		for (int j = 0; j < trueNumPoints; j++)
			evalPoints[j] = a + delta * (0.5 + j);

		return evalPoints;
	}

}
