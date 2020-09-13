package umontreal.ssj.stat.density;

import umontreal.ssj.stat.ScaledHistogram;
import umontreal.ssj.stat.TallyHistogram;

/**
 * Histogram density estimator for a univariate density.
 * 
 * This class provides methods to construct, manipulate, and evaluate a
 * histogram from a set of \f$n\f$ real-valued observations \f$x_0,\dots,
 * x_{n-1}\f$. The constructor of a `DEHistogram` object takes the data, a
 * finite interval \f$[a,b]\f$ over which we want the histogram, and the number
 * of bins (intervals), and it constructs the histogram as a density estimator
 * over this interval. For that, it computes the number of observations in each
 * bin and rescales the heights of the bins so that the total area of the
 * histogram is equal to the proportion of the observations that fall into the
 * interval \f$[a,b]\f$. When all the observations are in \f$[a,b]\f$, then this
 * area should be 1. The density estimator is represented internally as a \ref
 * umontreal.ssj.stat.ScaledHistogram, which can be recovered by the
 * `getScaledHistogram()` method.
 *
 * 
 * More specifically, if we ask for \f$s\f$ bins over the interval \f$[a,b]\f$,
 * then the interval is partitioned over \f$s\f$ intervals of equal lengths
 * \f$h=(b-a)/s\f$, and the density estimator is defined by
 * 
 * \f[ \hat{f}_{n}(x) = \hat{f}_{n,h}(x) = \frac{n_j}{nh}, \quad\text{ for }
 * x\in[a+(j-1)h, a+jh), \qquad j=1,\dots,s. \f]
 * 
 * where \f$n_j\f$ denotes the number of observations that fall in interval
 * \f$j\f$. Observe that changing \f$[a,b]\f$, \f$h\f$, or \f$s\f$ changes the
 * structure of the histogram completely. Thus, after altering any of these
 * parameters the histogram needs to be constructed afresh.
 * 
 * 
 * The constructor can take as input an array that contains the raw data,
 * together with the parameters \f$a, b, s\f$. It can also take a 
 * \ref umontreal.ssj.stat.TallyHistogram or a 
 * \ref umontreal.ssj.stat.ScaledHistogram instead.
 * 
 * This class also offers several static methods so that the user can simply evaluate the density based on
 * a set of observations, a `TallyHistogram` or a `ScaledHistogram` 
 * without having to construct a histogram object.
 * 
 * Since histograms are constant within one bin, it is sometimes sufficient to
 * evaluate them only once per bin. To this end, this class provides methods to
 * evaluate the density estimator once in each bin, without having to explicitly
 * specify an array of evaluation points.
 * 
 * It needs to be added that, due to the fact that the partition into \f$s\f$
 * bins relies on half-open intervals, the boundary \f$b\f$ is not included in
 * any of these intervals. Since the probability of an observation being exactly
 * equal to \f$b\f$ is zero, we can effectively ignore this subtlety.
 * 
 * This class also offers static methods so that the user can simply evaluate
 * the density based on a set of observations, a \ref
 * umontreal.ssj.stat.TallyHistogram or a \ref
 * umontreal.ssj.stat.ScaledHistogram without having to construct a histogram.
 *
 */

public class DEHistogram extends DensityEstimator {

	/** the actual histogram */
	private ScaledHistogram histDensity;

	/**
	 * Constructs a histogram estimator over the interval \f$[a,b]\f$ with \a
	 * numBins number of bins.
	 * 
	 * @param a
	 *            the left boundary of the histogram.
	 * @param b
	 *            the right boundary of the histogram.
	 * @param numBins
	 *            the number of bins.
	 */
	public DEHistogram(double a, double b, int numBins) {
		histDensity = new ScaledHistogram(a, b, numBins);
	}

	/**
	 * Constructs a histogram over the interval \f$[a,b]\f$ with \a numBins number
	 * of bins from the observations passed in \a data.
	 * 
	 * @param data
	 *            the observations from the underlying model.
	 * 
	 * @param a
	 *            the left boundary of the histogram.
	 * @param b
	 *            the right boundary of the histogram.
	 * @param numBins
	 *            the number of bins.
	 */
	public DEHistogram(double[] data, double a, double b, int numBins) {
		// this.data = data;
		TallyHistogram tallyHist = new TallyHistogram(a, b, numBins);
		tallyHist.fillFromArray(data);
		histDensity = new ScaledHistogram(tallyHist, tallyHist.getProportionInBoundaries());
	}

	/**
	 * Constructs a histogram from a \ref umontreal.ssj.stat.TallyHistogram \a
	 * tallyHist.
	 * 
	 * @param tallyHist
	 *            a \ref umontreal.ssj.stat.TallyHistogram from which the estimator
	 *            is constructed.
	 */
	public DEHistogram(TallyHistogram tallyHist) {

		histDensity = new ScaledHistogram(tallyHist, tallyHist.getProportionInBoundaries());
	}

	/**
	 * Constructs a histogram from a \ref umontreal.ssj.stat.ScaledHistogram \a
	 * scaledHist.
	 * 
	 * @param scaledHist
	 *            a \ref umontreal.ssj.stat.ScaledHistogram from which the estimator
	 *            is constructed.
	 */
	public DEHistogram(ScaledHistogram scaledHist) {
		histDensity = scaledHist;
	}

	/**
	 * Takes the observations in \a data and constructs and redefines the histogram
	 * with these observations.
	 * 
	 * @param data
	 *            the observations to define the histogram.
	 * 
	 */
	@Override
	public void setData(double[] data) {
		TallyHistogram tallyHist = new TallyHistogram(histDensity.getA(), histDensity.getB(), histDensity.getNumBins());
		tallyHist.fillFromArray(data);
		histDensity = new ScaledHistogram(tallyHist, tallyHist.getProportionInBoundaries());
	}

	/**
	 * Returns the underlying `ScaledHistogram`.
	 * 
	 * @return underlying `ScaledHistogram` object.
	 */
	public ScaledHistogram getScaledHistogram() {
		return histDensity;
	}

	/**
	 * Gives the number of bins \f$s\f$.
	 * 
	 * @return the number of bins.
	 */
	public int getNumBins() {
		return histDensity.getNumBins();
	}

	/**
	 * Gives the bin width \f$h\f$.
	 * 
	 * @return the bin width.
	 */
	public double getH() {
		return (getB() - getA()) / getNumBins();
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
	 * {@inheritDoc}
	 */
	@Override
	public String toString() {
		return "Histogram estimator with " + getNumBins() + " bins.";
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public double evalDensity(double x) {
		return histDensity.getHeights()[(int) ((x - getA()) / getH())];
	}

	/**
	 * Evaluates the histogram estimator at the points in \a evalPoints.
	 * 
	 * @param evalPoints
	 *            the evaluation points.
	 * 
	 * @return the histogram density estimator evaluated at the points \a x.
	 */
	@Override
	public double[] evalDensity(double[] evalPoints) {
		int k = evalPoints.length;
		double[] density = new double[k];

		for (int j = 0; j < k; j++)
			density[j] = histDensity.getHeights()[(int) ((evalPoints[j] - getA()) / getH())];

		return density;

	}

	/**
	 * Evaluates the histogram density estimator at one point in each bin. Since the
	 * histogram is constant within its bins, it is of no concern for this method,
	 * which evaluation point is actually selected in each bin.
	 * 
	 * 
	 * @return the histogram density estimator evaluated at one point in each bin.
	 */
	public double[] evalDensity() {
		return histDensity.getHeights();
	}

	/**
	 * Sets the number of bins according to the provided binwidth @f$h@f$ via
	 * \f$\lfloor (b-a)/h \rfloor\f$. Note that the actual binwidth might differ
	 * from the provided @f$h@f$ due to the floor-function.
	 * 
	 * @param h The desired binwidth.
	 */
	public void setH(double h) {
		int numBins = (int) ((getB() - getA()) / h);
		TallyHistogram tallyHist = new TallyHistogram(histDensity.getA(), histDensity.getB(), numBins);
		tallyHist.fillFromArray(data);
		histDensity = new ScaledHistogram(tallyHist, tallyHist.getProportionInBoundaries());

	}

	// STATIC METHODS

	/**
	 * Evaluates the histogram with \a numBins bins over \f$[a,b]\f$ which is
	 * defined by the observations \a data at the evaluation point \a x.
	 * 
	 * @param x
	 *            the evaluation point.
	 * @param data
	 *            the observations.
	 * @param a
	 *            the left boundary of the histogram.
	 * @param b
	 *            the right boundary of the histogram.
	 * @param numBins
	 *            the number of bins.
	 * @return the histogram defined be the above parameters evaluated at \a x.
	 */
	public static double evalDensity(double x, double[] data, double a, double b, int numBins) {
		TallyHistogram tallyHist = new TallyHistogram(a, b, numBins);
		tallyHist.fillFromArray(data);
		return evalDensity(x, tallyHist);
	}

	/**
	 * Evaluates the histogram with \a numBins bins over \f$[a,b]\f$ which is
	 * defined by the observations \a data at each of the evaluation points \a
	 * evalPoints and returns the results in an array.
	 * 
	 * @param evalPoints
	 *            the evaluation points.
	 * @param data
	 *            the observations.
	 * @param a
	 *            the left boundary of the histogram.
	 * @param b
	 *            the right boundary of the histogram.
	 * @param numBins
	 *            the number of bins.
	 * @return the histogram defined by the above parameters evaluated at the points
	 *         in \a evalPoints
	 */
	public static double[] evalDensity(double[] evalPoints, double[] data, double a, double b, int numBins) {
		TallyHistogram tallyHist = new TallyHistogram(a, b, numBins);
		tallyHist.fillFromArray(data);
		return evalDensity(evalPoints, tallyHist);
	}

	/**
	 * Evaluates the histogram with \a numBins bins over \f$[a,b]\f$ which is
	 * defined by the observations \a data once in each bin and returns the results
	 * in an array.
	 * 
	 * @param data
	 *            the observations.
	 * @param a
	 *            the left boundary of the histogram.
	 * @param b
	 *            the right boundary of the histogram.
	 * @param numBins
	 *            the number of bins.
	 * @return the histogram defined by the above parameters evaluated once in each
	 *         bin.
	 */

	public static double[] evalDensity(double[] data, double a, double b, int numBins) {
		TallyHistogram tallyHist = new TallyHistogram(a, b, numBins);
		tallyHist.fillFromArray(data);
		return evalDensity(tallyHist);
	}

	/**
	 * Same as `evalDensity(double[], double[][], double, double, int)` but here,
	 * the density is evaluated once in each bin instead of at a given array of
	 * evaluation points.
	 * 
	 * @param data
	 *            the observations.
	 * @param a
	 *            the left boundary of the histogram.
	 * @param b
	 *            the right boundary of the histogram.
	 * @param numBins
	 *            the number of bins.
	 * @return the histogram defined by the above parameters evaluated once in each
	 *         bin.
	 * @return the histogram for each realization evaluated once in each bin.
	 */
	public static double[][] evalDensity(double[][] data, double a, double b, int numBins) {
		int m = data.length;
		double[][] dens = new double[m][numBins];
		TallyHistogram tallyHist = new TallyHistogram(a, b, numBins);
		for (int r = 0; r < m; r++) {
			tallyHist.fillFromArray(data[r]);
			dens[r] = evalDensity(tallyHist);
		}
		return dens;
	}

	/**
	 * Assume that we have \f$m\f$ independent realizations of the underlying model.
	 * For each such realization this method evaluates a histogram with \a numBins
	 * bins over \f$[a,b]\f$ once in each bin. The independent realizations are
	 * passed via the 2-dimensional \f$m\times n\f$ array \a data, where \f$n\f$
	 * denotes the number of observations per realization. Hence, its first index
	 * identifies the independent realization while its second index identifies a
	 * specific observation of this realization.
	 * 
	 * The result is returned as a \f$m\times k\f$ matrix, where \f$k \f$ is the
	 * number of evaluation points, i.e., the length of \a evalPoints. The first
	 * index, again, identifies the independent realization whereas the second index
	 * corresponds to the point of \a evalPoints at which the histogram was
	 * evaluated.
	 * 
	 * @param evalPoints
	 *            the evaluation points.
	 * @param data
	 *            the two-dimensional array of observations.
	 * @param a
	 *            the left boundary of the histogram.
	 * @param b
	 *            the right boundary of the histogram.
	 * @param numBins
	 *            the number of bins.
	 * @return the histogram for each realization evaluated at \a evalPoints.
	 * @return
	 */
	public static double[][] evalDensity(double evalPoints[], double[][] data, double a, double b, int numBins) {
		int m = data.length;
		double[][] dens = new double[m][numBins];
		TallyHistogram tallyHist = new TallyHistogram(a, b, numBins);
		for (int r = 0; r < m; r++) {
			tallyHist.fillFromArray(data[r]);
			dens[r] = evalDensity(evalPoints, tallyHist);
		}
		return dens;
	}

	/**
	 * Evaluates a histogram estimator defined by \a tallylist at \a x.
	 * 
	 * @param x
	 *            the evaluation point.
	 * @param tallyHist
	 *            the \ref umontreal.ssj.stat.TallyHistogram which defines the
	 *            histogram
	 * @return the histogram density estimator evaluated at \a x.
	 */
	public static double evalDensity(double x, TallyHistogram tallyHist) {
		ScaledHistogram hist = new ScaledHistogram(tallyHist, tallyHist.getProportionInBoundaries());
		double h = (hist.getB() - hist.getA()) / (double) hist.getNumBins();
		return hist.getHeights()[(int) ((x - hist.getA()) / h)];
	}

	/**
	 * Evaluates a histogram estimator defined by \a tallyHist at the evaluation
	 * points \a evalPoints.
	 * 
	 * @param evalPoints
	 *            the evaluation points.
	 * @param tallyHist
	 *            the \ref umontreal.ssj.stat.TallyHistogram which defines the
	 *            histogram.
	 * @return the histogram estimator defined by \a tallyHist evaluated at each
	 *         point in \a evalPoints.
	 */
	public static double[] evalDensity(double[] evalPoints, TallyHistogram tallyHist) {

		ScaledHistogram hist = new ScaledHistogram(tallyHist, tallyHist.getProportionInBoundaries());
		int k = evalPoints.length;
		double[] density = new double[k];
		double h = (hist.getB() - hist.getA()) / (double) hist.getNumBins();
		for (int j = 0; j < k; j++) {
			density[j] = hist.getHeights()[(int) ((evalPoints[j] - hist.getA()) / h)];
		}
		return density;
	}

	/**
	 * Evaluates a histogram estimator defined by \a tallyHist once in each bin.
	 * 
	 * @param tallyHist
	 *            the \ref umontreal.ssj.stat.TallyHistogram which defines the
	 *            histogram.
	 * @return the histogram density estimator defined by \a tallyHist evaluated
	 *         once in each bin.
	 */
	public static double[] evalDensity(TallyHistogram tallyHist) {
		ScaledHistogram hist = new ScaledHistogram(tallyHist, tallyHist.getProportionInBoundaries());
		return hist.getHeights();
	}

	/**
	 * This method considers a histogram for each of the \f$m\f$ \ref
	 * umontreal.ssj.stat.TallyHistogram from \a tallyHistArray and evaluates it at
	 * \f$k\f$ evaluation points \a evalPoints. The result is returned in a
	 * \f$m\times k\f$ array.
	 * 
	 * @param evalPoints
	 *            the evaluation points.
	 * @param tallyHistArray
	 *            the array of the \ref umontreal.ssj.stat.TallyHistogram which
	 *            define the histograms.
	 * @return the histogram for each element of \a tallyHistArray evaluated at each
	 *         point of \a evalPoints.
	 */
	public static double[][] evalDensity(double[] evalPoints, TallyHistogram[] tallyHistArray) {
		int m = tallyHistArray.length;
		int k = evalPoints.length;
		double[][] density = new double[m][k];
		for (int r = 0; r < m; r++)
			density[r] = evalDensity(evalPoints, tallyHistArray[r]);
		return density;
	}

	/**
	 * Same as `evalDensity(double[], TallyHistogram[])` but this method evaluates
	 * each histogram only once in each bin.
	 * 
	 * @param tallyHistArray
	 *            the array of the \ref umontreal.ssj.stat.TallyHistogram which
	 *            define the histograms.
	 * @return the histogram for each element of \a tallyHistArray evaluated once in
	 *         each bin.
	 * 
	 */
	public static double[][] evalDensity(TallyHistogram[] tallyHistArray) {
		int m = tallyHistArray.length;
		double[][] dens = new double[m][];
		for (int r = 0; r < m; r++)
			dens[r] = evalDensity(tallyHistArray[r]);
		return dens;
	}

	/**
	 * Evaluates the histogram defined by the \ref
	 * umontreal.ssj.stat.ScaledHistogram at \a x.
	 * 
	 * @param x
	 *            the evaluation point.
	 * @param scaledHist
	 *            the \ref umontreal.ssj.stat.ScaledHistogram which defines the
	 *            histogram.
	 * @return the histogram density estimator evaluated at \a x.
	 */

	public static double evalDensity(double x, ScaledHistogram scaledHist) {
		double h = (scaledHist.getB() - scaledHist.getA()) / (double) scaledHist.getNumBins();
		return scaledHist.getHeights()[(int) ((x - scaledHist.getA()) / h)];
	}

	/**
	 * Evaluates the histogram defined by the \ref
	 * umontreal.ssj.stat.ScaledHistogram at each evaluation point in \a evalPoints.
	 * 
	 * @param evalPoints
	 *            the evaluation points.
	 * @param scaledHist
	 *            the \ref umontreal.ssj.stat.ScaledHistogram which defines the
	 *            histogram.
	 * @return the histogram density estimator evaluated at each point in \a
	 *         evalPoints.
	 */
	public static double[] evalDensity(double[] evalPoints, ScaledHistogram scaledHist) {
		int k = evalPoints.length;
		double[] density = new double[k];
		double h = (scaledHist.getB() - scaledHist.getA()) / (double) scaledHist.getNumBins();
		for (int j = 0; j < k; j++)
			density[k] = scaledHist.getHeights()[(int) ((evalPoints[j] - scaledHist.getA()) / h)];
		return density;
	}

	/**
	 * Same as `evalDensity(double[], ScaledHistogram)` but evaluation is done once
	 * in each bin.
	 * 
	 * @param scaledHist
	 *            \ref umontreal.ssj.stat.ScaledHistogram which defines the
	 *            histogram.
	 * @return the histogram density estimator evaluated at each point in \a
	 *         evalPoints.
	 */
	public static double[] evalDensity(ScaledHistogram scaledHist) {
		return scaledHist.getHeights();
	}

	/**
	 * This method considers a histogram for each of the \f$m\f$ \ref
	 * umontreal.ssj.stat.ScaledHistogram from \a scaledHistArray and evaluates it
	 * at \f$k\f$ evaluation points \a evalPoints. The result is returned in a
	 * \f$m\times k\f$ array.
	 * 
	 * @param evalPoints
	 *            the evaluation points.
	 * @param scaledHistArray
	 *            the array of the \ref umontreal.ssj.stat.ScaledHistogram which
	 *            define the histograms.
	 * @return the histogram for each element of \a scaledHistArray evaluated at
	 *         each point of \a evalPoints.
	 */
	public static double[][] evalDensity(double[] evalPoints, ScaledHistogram[] scaledHistArray) {
		int m = scaledHistArray.length;
		double[][] density = new double[m][];

		for (int r = 0; r < m; r++)
			density[m] = evalDensity(evalPoints, scaledHistArray[r]);

		return density;
	}

	/**
	 * Same as `evalDensity(double[], ScaledHistogram[])` but the histograms are
	 * evaluated once in each bin.
	 * 
	 * @param scaledHistArray
	 *            the array of the \ref umontreal.ssj.stat.ScaledHistogram which
	 *            define the histograms.
	 * @return the histogram for each element of \a scaledHistArray evaluated once in
	 *         each bin.
	 */
	public static double[][] evalDensity(ScaledHistogram[] scaledHistArray) {
		int m = scaledHistArray.length;
		double[][] density = new double[m][];
		for (int r = 0; r < m; r++)
			density[r] = evalDensity(scaledHistArray[r]);

		return density;
	}

}
