package umontreal.ssj.stat.density;

import umontreal.ssj.stat.ScaledHistogram;
import umontreal.ssj.stat.TallyHistogram;

/**
 * This class provides methods to construct, manipulate, and evaluate a histogram density estimator for
 * univariate densities. 
 *  
 * The construction of a histogram is based on a set of \f$n\f$
 * individual observations \f$X_0, …, X_{n-1}\f$ which, in turn, can be realizations of a
 *  \ref umontreal.ssj.mcqmctools.MonteCarloModelDouble, for instance. Note that it is essential 
 *  to confine oneself to a finite interval \f$[a,b]\f$.
 *
 * For a fixed number of bins \f$m>0\f$ we partition the interval \f$[a,b]\f$
 * into \f$m\f$ subintervals of equal lengths \f$h\f$. Observe that
 * \f$h=(b-a)/m\f$ and, hence, Histogram estimators can also be parameterized by a
 * binwidth \f$h>0\f$, as long as the resulting number of bins is an integer.
 * The estimator itself is defined by 
 * 
 *  \f[ 
 *  \hat{f}_{n}(x) = \hat{f}_{n,h}(x) =
 * \frac{n_j}{nh},\quad\text{for } x\in[a+(j-1)h, a+jh), \quad j=1,\dots,m,
 *  \f] 
 *  
 * where \f$n_j\f$ denotes the number of observations that fall in this interval.
 * 
 * It needs to be added that, due to the fact that the partition into \f$m\f$ bins relies on
 * half-open intervals, the boundary \f$b\f$ is not included in any of these
 * intervals. Since the probability of an observation being exactly equal to
 * \f$b\f$ is zero, we can effectively ignore this subtlety.
 * 
 * @author puchhamf
 *
 */

public class DEHistogram extends DensityEstimator {
	
	/**the number of bins */
	private int m;
	
	
	/**the actual histogram*/
	private ScaledHistogram histDensity;
//	private TallyHistogram hist;

	/**
	 * Constructor for a histogram estimator with \a m bins.
	 * 
	 * @param m the number of bins.
	 */
	public DEHistogram(int m) {
		this.m = m;
	}

	
	

	/**
	 * Gives the number of bins.
	 * 
	 * @return the number of bins.
	 */
	public int getM() {
		return m;
	}

	public double getH() {
		return (getB() - getA())/m;
	}
	
	public double getA() {
		return histDensity.getA();
	}
	
	public double getB() {
		return histDensity.getB();
	}

	

	
	public void constructDensity(double[] data, double a, double b) {
//		hist = new TallyHistogram(geta(), getb(), m);
		TallyHistogram hist = new TallyHistogram(a, b, m);
		constructDensity(hist);
	}
	
	/**
	 * Constructs a histogram density estimator from a \ref TallyHistogram. Furthermore, it automatically
	 * sets the parameters, {@link #a}, {@link #b} and {@link #m} to the values stored in \a tallyHist.
	 * @param tallyHist 
	 */
	public void constructDensity(TallyHistogram tallyHist) {
		m = tallyHist.getNumBins();
		histDensity = new ScaledHistogram(tallyHist, 1.0);
	}
	
	/**
	 * Constructs a histogram density estimator from a \ref ScaledHistogram. Furthermore, it automatically
	 * sets the parameters, {@link #a}, {@link #b} and {@link #m} to the values stored in \a scaledHist.
	 * @param hist
	 */
	public void constructDensity(ScaledHistogram scaledHist) {
		m = scaledHist.getNumBins();
		histDensity = scaledHist;
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
	protected double[] getEquidistantPoints(int numPoints) {
		int trueNumPoints = Math.max(numPoints, m);
		double evalPoints[] = new double[trueNumPoints];
		double delta = (b - a) / (trueNumPoints);
		for (int j = 0; j < trueNumPoints; j++)
			evalPoints[j] = a + delta * (0.5 + j);

		return evalPoints;
	}

	@Override
	public double evalDensity(double x, double[] data, double a, double b) {
		constructDensity(data,m,a,b);
		return 0;
	}

}
