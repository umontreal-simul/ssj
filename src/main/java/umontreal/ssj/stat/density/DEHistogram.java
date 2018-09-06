package umontreal.ssj.stat.density;

import umontreal.ssj.stat.ScaledHistogram;
import umontreal.ssj.stat.TallyHistogram;

/**
 * This class provides methods to construct, manipulate, and evaluate a
 * histogram density estimator for univariate densities.
 * 
 * The construction of a histogram is based on a set of \f$n\f$ individual
 * observations \f$X_0, …, X_{n-1}\f$ which, in turn, can be realizations of a
 * \ref umontreal.ssj.mcqmctools.MonteCarloModelDouble, for instance. In any
 * case, this type of density estimator can by built from pure data stored in an
 * array, form a \ref umontreal.ssj.stat.TallyHistogram, or even from a \ref
 * umontreal.ssj.stat.ScaledHistogram.
 * 
 * 
 * It is essential to confine oneself to a finite interval \f$[a,b]\f$. For a
 * fixed number of bins \f$s>0\f$ we partition the interval \f$[a,b]\f$ into
 * \f$s\f$ subintervals of equal lengths \f$h\f$. Observe that \f$h=(b-a)/s\f$
 * and, hence, histogram estimators can also be parameterized by a bin width
 * \f$h>0\f$, as long as the resulting number of bins is an integer. The
 * estimator itself is then defined by
 * 
 * \f[ \hat{f}_{n}(x) = \hat{f}_{n,h}(x) = \frac{n_j}{nh},\quad\text{for }
 * x\in[a+(j-1)h, a+jh), \quad j=1,\dots,s, \f]
 * 
 * where \f$n_j\f$ denotes the number of observations that fall in this
 * interval. Observe that changing \f$[a,b]\f$, \f$h\f$, or \f$s\f$ changes the
 * structure of the histogram completely. Thus, after any alteration of these
 * parameters the histogram needs to be constructed afresh.
 * 
 * As histograms are constant within one bin certain quantities, such as the
 * variance, are not affected by shifting an evaluation point within its bin.
 * Hence, selecting one evaluation point per bin is sufficient for such methods.
 * To this end, this class provides methods to evaluate the density estimator
 * once in each bin, without having to explicitly specify an array of evaluation
 * points.
 * 
 * It needs to be added that, due to the fact that the partition into \f$s\f$
 * bins relies on half-open intervals, the boundary \f$b\f$ is not included in
 * any of these intervals. Since the probability of an observation being exactly
 * equal to \f$b\f$ is zero, we can effectively ignore this subtlety.
 * 
 * @author puchhamf
 *
 */

public class DEHistogram extends DensityEstimator {

	/** the number of bins \f$s\f$ */
	private int numBins;

	/** the actual histogram */
	private ScaledHistogram histDensity;
	// private TallyHistogram hist;

	/**
	 * Constructor for a histogram estimator with \a numBins bins.
	 * 
	 * @param numBins
	 *            the number of bins.
	 */
	public DEHistogram(int numBins) {
		this.numBins = numBins;
	}

	/**
	 * Gives the number of bins \f$s\f$.
	 * 
	 * @return the number of bins.
	 */
	public int getNumBins() {
		return numBins;
	}

	/**
	 * Sets the number of bins to \a numBins. Note that this is only meaningful if
	 * the entire histogram will be reconstructed.
	 * 
	 * @param numBins
	 *            the number of bins.
	 */
	public void setNumBins(int numBins) {
		this.numBins = numBins;
	}

	/**
	 * Gives the bin width \f$h\f$.
	 * 
	 * @return the desired bin width.
	 */
	public double getH() {
		return (getB() - getA()) / numBins;
	}

	/**
	 * Sets the bin width to \a h. The actual value used may differ from the one
	 * being passed since \a numBins\f$=(b-a)/h\f$ has to be an integer. Note that
	 * setting the bin width is only meaningful if the entire histogram will be
	 * reconstructed.
	 * 
	 * @param h
	 *            the desired bin width.
	 */
	public void setH(double h) {
		this.numBins = (int) ((getB() - getA()) / h);
	}

	/**
	 * Gives the left boundary \f$a\f$ of the histogram.
	 * 
	 * @return the left boundary of the histogram.
	 */
	public double getA() {
		return histDensity.getA();
	}

	/**
	 * Gives the right boundary \f$b\f$ of the histogram.
	 * 
	 * @return the right boundary of the histogram.
	 */
	public double getB() {
		return histDensity.getB();
	}

	/**
	 * Constructs a histogram with #numBins bins over the interval \f$[a,b]\f$ from
	 * the observations passed in \a data. Note that the individual observations are
	 * not being stored.
	 * 
	 * @param data
	 *            the observations from the underlying model.
	 * @param a
	 *            the left boundary of the histogram.
	 * @param b
	 *            the right boundary of the histogram
	 */
	public void constructDensity(double[] data, double a, double b) {
		// hist = new TallyHistogram(geta(), getb(), m);
		TallyHistogram hist = new TallyHistogram(a, b, numBins);
		constructDensity(hist);
	}

	/**
	 * Constructs a histogram from a \ref umontreal.ssj.stat.TallyHistogram \a
	 * tallyHist. The method extracts all defining parameters such as the number of
	 * bins, the endpoints, etc. from \a tallyHist directly.
	 * 
	 * @param tallyHist
	 *            a \ref umontreal.ssj.stat.TallyHistogram from which the estimator
	 *            is constructed.
	 */
	public void constructDensity(TallyHistogram tallyHist) {
		numBins = tallyHist.getNumBins();
		histDensity = new ScaledHistogram(tallyHist, 1.0);
	}

	/**
	 * Constructs a histogram from a \ref umontreal.ssj.stat.ScaledHistogram \a
	 * scaledHist. The method extracts all defining parameters such as the number of
	 * bins, the endpoints, etc. from \a scaledHist directly.
	 * 
	 * @param scaledHist
	 *            a \ref umontreal.ssj.stat.ScaledHistogram from which the estimator
	 *            is constructed.
	 */
	public void constructDensity(ScaledHistogram scaledHist) {
		numBins = scaledHist.getNumBins();
		histDensity = scaledHist;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String toString() {
		return "Histogram estimator with " + numBins + " bins.";
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public double evalDensity(double x, double[] data, double a, double b) {
		constructDensity(data, a, b);
		return histDensity.getHeights()[(int) ((x - a) / getH())];
	}

	/**
	 * Constructs a histogram density estimator over the interval \f$[a,b]\f$ based
	 * on the observations \a data and evaluates it at the points in \a evalPoints.
	 * 
	 * @param evalPoints
	 *            the evaluation points.
	 * @param data
	 *            the observations for constructing the histogram.
	 * @param a
	 *            the left boundary of the interval.
	 * @param b
	 *            the right boundary of the interval
	 * @return the histogram density estimator evaluated at the points \a x.
	 */
	@Override
	public double[] evalDensity(double[] evalPoints, double[] data, double a, double b) {
		int k = evalPoints.length;
		double[] density = new double[k];

		constructDensity(data, a, b);

		for (int j = 0; j < k; j++)
			density[j] = histDensity.getHeights()[(int) ((evalPoints[j] - getA()) / getH())];

		return density;

	}

	/**
	 * Constructs a histogram density estimator over the interval \f$[a,b]\f$ based
	 * on the observations \a data and evaluates it at one point in each bin. Since
	 * the histogram is constant within its bins, it is of no concern for this
	 * method, which evaluation point is actually selected in each bin.
	 * 
	 * @param data
	 *            the observations for constructing the histogram.
	 * @param a
	 *            the left boundary of the interval.
	 * @param b
	 *            the right boundary of the interval
	 * @return the histogram density estimator evaluated at the points \a x.
	 */
	public double[] evalDensity(double[] data, double a, double b) {
		constructDensity(data, a, b);
		return histDensity.getHeights();
	}

	/**
	 * Same as umontreal.ssj.stat.density.DensityEstimator#evalDensity(double[],
	 * double[][], double, double) but with one evaluation point per bin. Since the
	 * histogram is constant within its bins, it is of no concern for this method,
	 * which evaluation point is actually selected in each bin.
	 * 
	 * @param data
	 *            the two-dimensional array carrying the observations of \f$m\f$
	 *            independent realizations of the underlying model.
	 * @param a
	 *            the left boundary of the interval.
	 * @param b
	 *            the right boundary of the interval.
	 * @return the histogram density estimator for each realization evaluated once
	 *         in each bin.
	 */
	public double[][] evalDensity(double[][] data, double a, double b) {
		int m = data.length;
		double[][] dens = new double[m][];
		for (int r = 0; r < m; r++)
			dens[r] = evalDensity(data[r], a, b);
		return dens;
	}

	/**
	 * Constructs a histogram estimator from \a tallyHist and evaluates it at \a x.
	 * 
	 * @param x
	 *            the evaluation point.
	 * @param tallyHist
	 *            the \ref umontreal.ssj.stat.TallyHistogram from which the
	 *            estimator is constructed.
	 * @return the histogram density estimator evaluated at \a x.
	 */
	public double evalDensity(double x, TallyHistogram tallyHist) {
		constructDensity(tallyHist);
		return histDensity.getHeights()[(int) ((x - getA()) / getH())];
	}

	/**
	 * Constructs a histogram density estimator from \a tallyHist and evaluates it
	 * at one point in each bin. Since the histogram is constant within its bins, it
	 * is of no concern for this method, which evaluation point is actually selected
	 * in each bin.
	 * 
	 * @param tallyHist
	 *            the \ref umontreal.ssj.stat.TallyHistogram from which the
	 *            estimator is constructed.
	 * @return the histogram density estimator evaluated once in each bin.
	 */
	public double[] evalDensity(TallyHistogram tallyHist) {
		constructDensity(tallyHist);
		return histDensity.getHeights();
	}

	/**
	 * Constructs a histogram estimator from \a tallyHist and evaluates it at each
	 * evaluation point in \a evalPoints.
	 * 
	 * @param evalPoints
	 *            the evaluation points.
	 * @param tallyHist
	 *            the \ref umontreal.ssj.stat.TallyHistogram from which the
	 *            estimator is constructed.
	 * @return the histogram density estimator evaluated at \a evalPoints.
	 */
	public double[] evalDensity(double[] evalPoints, TallyHistogram tallyHist) {
		int k = evalPoints.length;
		double[] density = new double[k];

		constructDensity(tallyHist);

		for (int j = 0; j < k; j++)
			density[j] = histDensity.getHeights()[(int) ((evalPoints[j] - getA()) / getH())];

		return density;
	}

	/**
	 * Same as #evalDensity(double[][], double, double) but here, each realization
	 * of the histogram is constructed from a \ref
	 * umontreal.ssj.stat.TallyHistogram.
	 * 
	 * @param tallyHistArray
	 *            the array of \ref umontreal.ssj.stat.TallyHistogram from which the
	 *            realizations of the estimator is constructed.
	 * @return the histogram density estimator for each realization evaluated once
	 *         in each bin.
	 */
	public double[][] evalDensity(TallyHistogram[] tallyHistArray) {
		int m = tallyHistArray.length;
		double[][] dens = new double[m][];
		for (int r = 0; r < m; r++)
			dens[r] = evalDensity(tallyHistArray[r]);
		return dens;
	}

	/**
	 * Constructs a histogram estimator from \a scaledHist and evaluates it at \a x.
	 * 
	 * @param x
	 *            the evaluation point.
	 * @param scaledHist
	 *            the \ref umontreal.ssj.stat.ScaledHistogram from which the
	 *            estimator is constructed.
	 * @return the histogram density estimator evaluated at \a x.
	 */

	public double evalDensity(double x, ScaledHistogram scaledHist) {
		constructDensity(scaledHist);

		return histDensity.getHeights()[(int) ((x - getA()) / getH())];
	}

	/**
	 * Constructs a histogram density estimator from \a scaledHist and evaluates it
	 * at one point in each bin. Since the histogram is constant within its bins, it
	 * is of no concern for this method, which evaluation point is actually selected
	 * in each bin.
	 * 
	 * @param scaledHist
	 *            the \ref umontreal.ssj.stat.ScaledHistogram from which the
	 *            estimator is constructed.
	 * @return the histogram density estimator evaluated once in each bin.
	 */
	public double[] evalDensity(ScaledHistogram scaledHist) {
		constructDensity(scaledHist);

		return histDensity.getHeights();
	}

	/**
	 * Constructs a histogram estimator from \a scaledHist and evaluates it at each
	 * evaluation point in \a evalPoints.
	 * 
	 * @param evalPoints
	 *            the evaluation points.
	 * @param scaledHist
	 *            the \ref umontreal.ssj.stat.ScaledHistogram from which the
	 *            estimator is constructed.
	 * @return the histogram density estimator evaluated at \a evalPoints.
	 */
	public double[] evalDensity(double[] evalPoints, ScaledHistogram scaledHist) {
		int k = evalPoints.length;
		double[] density = new double[k];

		constructDensity(scaledHist);

		for (int j = 0; j < k; j++)
			density[j] = histDensity.getHeights()[(int) ((evalPoints[j] - getA()) / getH())];

		return density;
	}

	/**
	 * Same as #evalDensity(double[][], double, double) but here, each realization
	 * of the histogram is constructed from a \ref
	 * umontreal.ssj.stat.ScaledHistogram.
	 * 
	 * @param scaledHistArray
	 *            the array of \ref umontreal.ssj.stat.ScaledHistogram from which
	 *            the realizations of the estimator is constructed.
	 * @return the histogram density estimator for each realization evaluated once
	 *         in each bin.
	 */
	public double[][] evalDensity(ScaledHistogram[] scaledHistArray) {
		int m = scaledHistArray.length;
		double[][] dens = new double[m][];
		for (int r = 0; r < m; r++)
			dens[r] = evalDensity(scaledHistArray[r]);
		return dens;
	}

	// TODO: do those methods make sense???
	// public double computeIV(double[][] densities, double[] variance) {
	// return computeIV(densities,getA(),getB(),variance);
	// }
	//
	// public void computeIV(ArrayList<double[][]> listDensity, ArrayList<Double>
	// listIV) {
	// computeIV(listDensity, getA(), getB(), listIV);
	// }
	//
	// public double[] computeMISE(ContinuousDistribution dist, double[] evalPoints,
	// double[][] density,double[] variance, double[] sqBias, double[] mse) {
	// return computeMISE(dist, evalPoints, density,
	// getA(),getB(),variance,sqBias,mse);
	// }
	//
	// public void computeMISE(ContinuousDistribution dist, double[] evalPoints,
	// ArrayList<double[][]> listDensity,
	// ArrayList<double[]> listMISE) {
	// computeMISE(dist, evalPoints, listDensity, getA(), getB(), listMISE);
	// }

}
