
package umontreal.ssj.stat.density;

import java.util.ArrayList;

import umontreal.ssj.probdist.ContinuousDistribution;
import umontreal.ssj.stat.PgfDataTable;

/**
 * This abstract class implements a univariate density estimator (DE). To this
 * end, it provides basic tools to evaluate the DE at one point \f$x\f$ or at an
 * array of points \f$\{x_1, x_2, \dots, x_k\} \f$. More precisely, the single
 * point evaluation #evalDensity(double, double[], double, double) is abstract,
 * since it will definitely differ between realizations. For the evaluation on a
 * set of points one can use #evalDensity(double[], double[], double, double).
 * Note that, even though a default implementation is provided, very often
 * specific estimators will have more efficient evaluation algorithms. So,
 * overriding this method can be benificial in many cases. Furthermore, this
 * class includes a method to plot the estimated density.
 * 
 * There are also several more involved methods covered by this class, most of
 * which are concerned with the convergence behavior of DEs. Nevertheless, they
 * can be useful in many other cases beyond convergence behavior. As these
 * usually require more than one realization of a DE, or even a list of DE's,
 * they are implemented as static methods. For instance, #evalDensity(double[],
 * double[][], double, double) allows to evaluate several independent
 * replications of a DE at an array of evaluation points, and
 * #evalDensity(ArrayList, double[], double[][], double, double, ArrayList) does
 * the same but for more than one DE.
 * 
 * For measuring the performance of a DE, we need to confine ourselves to
 * estimation over a finite interval \f$[a,b]\f$. One standard way to assess the
 * quality of a DE is via the mean integrated square error (MISE). It can be
 * further decomposed into the integrated variance (IV), and the integrated
 * square bias (ISB)
 * 
 * \f[ \textrm{MISE} = \int_a^b\mathbb{E} [\hat{f}(x) - f(x)]^2\mathrm{d}x =
 * \int_a^b\textrm{Var}[\hat{f}(x)]\mathrm{d}x + \int_a^b \left(
 * \mathbb{E}[\hat{f}(x)] - f(x) \right)^2\mathrm{d}x, \f]
 * 
 * where \f$f\f$ denotes the true density and \f$\hat{f}\f$ the DE.
 * 
 * This class also provides methods to estimate the empirical IV, see
 * #computeIV(double[][], double, double, double[]) for one individual DE and
 * #computeIV(ArrayList, double, double, ArrayList) for several DEs. Note that
 * these mehtods merely compute an estimate, since computing an exact integral
 * is out of reach for this class.
 * 
 * Note that the MISE can only be computed in situations where either the ISB is
 * zero or the true density is known. In the first case, the IV is the same as
 * the MISE, of course. For the second case this class provides the possibility
 * to estimate the empirical MISE for a single DE via
 * #computeMISE(ContinuousDistribution, double[], double[][], double, double,
 * double[], double[], double[])} as well as for several DEs via
 * #computeMISE(ContinuousDistribution, double[], ArrayList, double, double,
 * ArrayList). Again, this merely gives an estimate of the empirical MISE.
 */

public abstract class DensityEstimator {

	/**
	 * Constructs a density estimator over the interval \f$[a,b]\f$ based on the
	 * observations \a data, if necessary, and evaluates it at \a x.
	 * 
	 * @param x
	 *            the evaluation point.
	 * @param data
	 *            the observations for constructing the density estimator.
	 * @param a
	 *            the left boundary of the interval.
	 * @param b
	 *            the right boundary of the interval
	 * @return the density estimator evaluated at \f$x\f$.
	 */

	public abstract double evalDensity(double x, double data[], double a, double b);

	/**
	 * Constructs a density estimator over the interval \f$[a,b]\f$ based on the
	 * observations \a data, if necessary, and evaluates it at the points in \a
	 * evalPoints. By default, this method calls #evalDensity(double, double[],
	 * double, double) for each entry of \a evalPoints. Many density estimators can
	 * handle evaluation at a set of points more efficiently than that. If so, it is
	 * suggested to override this method in the implementation of the corresponding
	 * estimator. Note that the construction of the density estimator -- if
	 * applicable -- would need to be handled by this method as well.
	 * 
	 * @param evalPoints
	 *            the evaluation points.
	 * @param data
	 *            the observations for constructing the density estimator.
	 * @param a
	 *            the left boundary of the interval.
	 * @param b
	 *            the right boundary of the interval
	 * @return the density estimator evaluated at the points \a x.
	 */
	public double[] evalDensity(double[] evalPoints, double[] data, double a, double b) {
		int k = evalPoints.length;
		double[] dens = new double[k];
		for (int j = 0; j < k; j++)
			dens[j] = evalDensity(evalPoints[j], data, a, b);
		return dens;
	}

	/**
	 * This method is particularly designed to evaluate the density estimator in
	 * such a way that the result can be easily used to estimate the empirical IV
	 * and other convergence-related quantities.
	 * 
	 * Assume that we have \f$m\f$ independent realizations of the underlying model.
	 * For each such realization this method constructs a density estimator over
	 * \f$[a,b]\f$ and evaluates it at the points from \a evalPoints. The
	 * independent realizations are passed via the 2-dimensional \f$m\times
	 * n\f$array \a data, where \f$n\f$ denotes the number of observations per
	 * realization. Hence, its first index identifies the independent realization
	 * while its second index identifies a specific observation of this realization.
	 * 
	 * The result is returned as a \f$m\times k\f$ matrix, where \f$k \f$ is the
	 * number of evaluation points, i.e., the length of \a evalPoints. The first
	 * index, again, identifies the independent realization whereas the second index
	 * corresponds to the point of \a evalPoints at which the density estimator was
	 * evaluated.
	 * 
	 *
	 * @param evalPoints
	 *            the evaluation points.
	 * @param data
	 *            the two-dimensional array carrying the observations of \f$m\f$
	 *            independent realizations of the underlying model.
	 * @param a
	 *            the left boundary of the interval.
	 * @param b
	 *            the right boundary of the interval.
	 * @return the density estimator for each realization evaluated at \a evalPoints.
	 */
	public double[][] evalDensity(double[] evalPoints, double[][] data, double a, double b) {
		int m = data.length;
		double[][] dens = new double[m][];
		for (int r = 0; r < m; r++)
			dens[r] = evalDensity(evalPoints, data[r], a, b);
		return dens;
	}

	/**
	 * This function is particularly designed for experiments with many different
	 * types of density estimators, as it evaluates all of these estimators at the
	 * points in \a evalPoints. To this end, the user passes a list of density
	 * estimators in \a listDE as well as \f$m\f$ independent realizations of the
	 * underlying model consisting of \f$n\f$ observations each in the \f$m\times
	 * n\f$ array \a data.
	 * 
	 * This method then calls #evalDensity(double[], double[][], double, double) for
	 * each density estimator in \a listDE, thus evaluating the respective density
	 * estimator at the \f$k\f$ points in \a evalPoints and adds the resulting
	 * \f$m\times k\f$ array to \a listDensity.
	 * 
	 * @param listDE
	 *            the list of density estimators.
	 * @param evalPoints
	 *            the evaluation points.
	 * @param data
	 *            the two-dimensional array carrying the observations of \f$m\f$
	 *            independent realizations of the underlying model.
	 * @param a
	 *            the left boundary of the interval.
	 * @param b
	 *            the right boundary of the interval.
	 * @param listDensity
	 *            a list to which the evaluations at \a evalPoints of each density
	 *            estimator in \a listDE are added.
	 * 
	 * @remark **Florian:** I kept the return type as "void" instead of
	 *         "ArrayList<double[][]>" and pass the corresponding list \a
	 *         listDensity to allow for more flexibility when working with it.
	 */
	public static void evalDensity(ArrayList<DensityEstimator> listDE, double[] evalPoints, double[][] data, double a,
			double b, ArrayList<double[][]> listDensity) {
		// ArrayList<double[][]> listDensity = new ArrayList<double[][]>();
		for (DensityEstimator de : listDE)
			listDensity.add(de.evalDensity(evalPoints, data, a, b));
		// return listDensity;
	}

	/**
	 * This method computes the empirical variance based on the values given in \a
	 * data. More precisely, \a density is a \f$m\times k\f$ matrix, whose entries
	 * correspond to \f$m\f$ independent realizations of the density estimator, each
	 * evaluated at \f$k\f$ evaluation points. Such a matrix can, for instance, be
	 * obtained by #evalDensity(double[], double[][], double, double).
	 * 
	 * The empirical variance is computed at each of those \f$k \f$ evaluation
	 * points and returned in an array of size \f$k\f$.
	 * 
	 * @param density
	 *            the estimated density of \f$m\f$ independent realizations of the
	 *            estimator, each evaluated at \f$k\f$ evaluation points.
	 * @return the empirical variance at those \f$k\f$ evaluation points.
	 */
	public static double[] computeVariance(double[][] density) {
		int m = density.length; // number of indep. replications
		int k = density[0].length; // number of evaluation points

		double x, y;

		double meanDens[] = new double[k]; // Average value over the rep replicates
		double varDens[] = new double[k]; // Variance at each evaluation point
		// Arrays.fill(meanDens, 0.0);
		// Arrays.fill(varDens, 0.0);
		for (int r = 0; r < m; r++) {
			// Update the empirical mean and variance at each evaluation point.
			for (int j = 0; j < k; j++) {
				x = density[r][j];
				y = x - meanDens[j];
				meanDens[j] += y / (double) (r + 1);
				varDens[j] += y * (x - meanDens[j]);
			}
		}

		for (int j = 0; j < k; j++) // normalize
			varDens[j] /= (double) (m - 1);

		return varDens;
	}

	/**
	 * This method estimates the empirical IV over the interval \f$[a,b]\f$. Based
	 * on the density estimates of \f$m\f$ independent replications of the density
	 * estimator evaluated at \f$k\f$ evaluation points, which are provided by \a
	 * density, it computes the empirical variance at each evaluation point and
	 * stores it in \a variance.
	 * 
	 * To estimate the empirical IV, we sum up the variance at the evaluation points
	 * \f$x_1,x_2,\dots,x_k\f$ and multiply by \f$(b-a)/k\f$, i.e.
	 * 
	 * \f[ \int_a^b \hat{f}(x)\mathrm{d}x \approx \frac{b-a}{k} \sum_{j =
	 * 1}^k\hat{f}(x_j), \f]
	 * 
	 * 
	 * where \f$\hat{f}\f$ denotes the density estimator. In other words, we
	 * approximate the empirical IV by an equally weighted quadrature rule using the
	 * aforementioned evaluation points as integration nodes.
	 * 
	 * 
	 * Note that this is only an approximation of the true empirical IV and that the
	 * approximation quality significantly depends on the choice of evaluation
	 * points.
	 * 
	 * The data for the variance are given in the two-dimensional \f$m\times k\f$
	 * array \a density, which is also described in #computeVariance(double[][]) and
	 * can be obtained by #evalDensity(double[], double[][], double, double) . The
	 * boundaries of the interval are given by \a a and \a b. Note that the array \a
	 * variance needs to be of length \f$k\f$.
	 * 
	 * 
	 * @param density
	 *            the \f$m\times k\f$ array that contains the data of evaluating
	 *            \f$m\f$ replicates of the density estimator at \f$k\f$ evaluation
	 *            points
	 * @param a
	 *            the left boundary of the interval.
	 * @param b
	 *            the right boundary of the interval.
	 * @param variance
	 *            the array of length \f$k\f$ in which the variance at each
	 *            evaluation point is stored.
	 * @return the estimated empirical IV over \f$[a,b]\f$.
	 */
	public static double computeIV(double[][] density, double a, double b, double[] variance) {
		variance = computeVariance(density);
		int k = density[0].length;
		double iv = 0.0;
		for (double var : variance)
			iv += var;
		return iv * (b - a) / (double) k;
	}

	/**
	 * This method estimates the empirical IV over the interval \f$[a,b]\f$ for a
	 * collection of different estimators. In \a densityList the user passes a list
	 * of \f$m\times k\f$ arrays which contain the density estimates of \f$m\f$
	 * independent replications of each density estimator evaluated at \f$k\f$
	 * evaluation points. Such a list can be obtained via #evalDensity(ArrayList,
	 * double[], double[][], double, double, ArrayList) , for instance.
	 * 
	 * The method then calls #computeIV(double[][], double, double, double[]) for
	 * each element of \a densityList and adds the thereby obtained estimated
	 * empirical IV to the list that is being returned.
	 * 
	 * 
	 * @param listDensity
	 *            list containing \f$m\times k\f$ arrays that contain the data of
	 *            evaluating \f$m\f$ replicates of each density estimator at \f$k\f$
	 *            evaluation points \a evalPoints.
	 * @param a
	 *            the left boundary of the interval.
	 * @param b
	 *            the right boundary of the interval.
	 * @param listIV
	 *            the list to which the estimated empirical IV of each density
	 *            estimator will be added.
	 * 
	 * @remark **Florian:** I kept the return type as "void" instead of
	 *         "ArrayList<double[][]>" and pass the corresponding list \a
	 *         listDensity to allow for more flexibility when working with it.
	 */
	public static void computeIV(ArrayList<double[][]> listDensity, double a, double b, ArrayList<Double> listIV) {

		// ArrayList<Double> returnList = new ArrayList<Double>();
		int k = (listDensity.get(0))[0].length;
		double[] variance = new double[k];
		for (double[][] density : listDensity) {
			// returnList.add(computeIV(density, a, b, variance));
			listIV.add(computeIV(density, a, b, variance));
		}
	}

	/**
	 * In situations where the true density is known this method can estimate the
	 * empirical MISE over the interval \f$[a,b]\f$. This can be particularly
	 * interesting and useful for testing density estimators. Since it is necessary
	 * to compute either the ISB or the IV to get the MISE and as there is not much
	 * computational overhead to estimate the other, an array containing the
	 * estimated empirical IV, the ISB, and MISE in exactly this order is returned.
	 * Based on the density estimates of \f$m\f$ independent replications of the
	 * density estimator evaluated at \f$k\f$ evaluation points \a evalPoints, which
	 * are provided by \a density, it computes the empirical variance, the
	 * square-bias, and the mean square error (MSE) at each evaluation point and
	 * stores the result in \a variance, \a sqBias, and \a mse, respectively. It is
	 * important that the evaluation points in \a evalPoints are the same as the
	 * ones used to construct \a density.
	 * 
	 * To estimate the empirical IV and MISE we sum up the variance and the MSE at
	 * the \f$k\f$ evaluation points and multiply by \f$(b-a)/k\f$, i.e. we
	 * approximate the empirical IV by an equally weighted quadrature rule with \a
	 * evalPoints as integration nodes. The ISB is then computed as the difference
	 * of the MISE and the IV. Note that this is only an approximation of the true
	 * empirical values and that the approximation quality significantly depends on
	 * the choice of \a evalPoints.
	 * 
	 * The data for the variance and mse are given in the two-dimensional \f$m\times
	 * k\f$ array \a density, which is also described in
	 * #computeVariance(double[][]) and can be obtained by #evalDensity(double[],
	 * double[][], double, double), and the true density is passed via a \ref
	 * umontreal.ssj.probdist.ContinuousDistribution. The evaluation points are
	 * contained in \a evalPoints and the boundaries of the interval over which we
	 * estimate are given by \a a and \a b. Note that the arrays \a variance, \a
	 * sqBias, and \a mse need to be of length \f$k\f$.
	 * 
	 * @param dist
	 *            the true density.
	 * @param evalPoints
	 *            the \f$k\f$ evaluation points.
	 * @param density
	 *            the \f$m\times k\f$ array that contains the data of evaluating
	 *            \f$m\f$ replicates of the density estimator at \f$k\f$ evaluation
	 *            points \a evalPoints.
	 * @param a
	 *            the left boundary of the interval.
	 * @param b
	 *            the right boundary of the interval.
	 * @param variance
	 *            the array of length \f$k\f$ in which the variance at each
	 *            evaluation point is stored.
	 * @param sqBias
	 *            the array of length \f$k\f$ in which the square-bias at each
	 *            evaluation point is stored.
	 * @param mse
	 *            the array of length \f$k\f$ in which the MSE at each evaluation
	 *            point is stored.
	 * @return an array containing the estimated empirical IV, ISB, and MISE in
	 *         exactly this order.
	 */
	public static double[] computeMISE(ContinuousDistribution dist, double[] evalPoints, double[][] density, double a,
			double b, double[] variance, double[] sqBias, double[] mse) {
		int m = density.length;
		int k = evalPoints.length;

		double x, y, z;
		double trueDensity;

		double meanDens[] = new double[k]; // Average value over m
											// replicates

		for (int r = 0; r < m; r++) {
			// Update the empirical mean, sum of squares, and mse of
			// observations at each evaluation point.
			for (int j = 0; j < k; j++) {
				x = density[r][j];
				y = x - meanDens[j];
				trueDensity = dist.density(evalPoints[j]);
				z = x - trueDensity;

				meanDens[j] += y / (double) (r + 1);
				variance[j] += y * (x - meanDens[j]);
				mse[j] += z * z;
			}

		}

		double iv = 0.0;
		double mise = 0.0;
		for (int j = 0; j < k; j++) {
			variance[j] /= (double) (m - 1.0);
			mse[j] /= (double) m;
			sqBias[j] = mse[j] - variance[j];

			iv += variance[j];
			mise += mse[j];
		}

		double fact = (b - a) / (double) k;
		iv *= fact;
		mise *= (b - a) / ((double) k);

		double[] res = { iv, mise - iv, mise };
		return res;

	}

	/**
	 * This method estimates the empirical MISE over the interval \f$[a,b]\f$ for a
	 * collection of different estimators. This can be done when the true density is
	 * actually known and is particularly interesting and/or useful for testing
	 * density estimators.
	 * 
	 * In \a densityList the user passes a list of \f$m\times k\f$ arrays which
	 * contain the density estimates of \f$m\f$ independent replications of each
	 * density estimator evaluated at \f$k\f$ evaluation points. Such a list can be
	 * obtained by calling #evalDensity(ArrayList, double[], double[][], double,
	 * double), for instance.
	 * 
	 * The method then calls #computeMISE(ContinuousDistribution, double[],
	 * double[][], double, double, double[], double[], double[]) for each element of
	 * \a listDensity. This results in an array containing the estimated empirical
	 * IV, ISB, and MISE in exactly this order, which is then added to the list \a
	 * listMISE.
	 * 
	 * @param dist
	 *            the true density.
	 * @param evalPoints
	 *            the \f$k\f$ evaluation points.
	 * @param listDensity
	 *            list containing \f$m\times k\f$ arrays that contain the data of
	 *            evaluating \f$m\f$ replicates of each density estimator at \f$k\f$
	 *            evaluation points \a evalPoints.
	 * @param a
	 *            the left boundary of the interval.
	 * @param b
	 *            the right boundary of the interval.
	 * @param listMISE
	 *            a list to which the arrays containing the estimated empirical IV,
	 *            ISB, and MISE of each density estimator are added.
	 * @remark **Florian:** I kept the return type as "void" instead of
	 *         "ArrayList<double[][]>" and pass the corresponding list \a
	 *         listDensity to allow for more flexibility when working with it.
	 */

	public static void computeMISE(ContinuousDistribution dist, double[] evalPoints, ArrayList<double[][]> listDensity,
			double a, double b, ArrayList<double[]> listMISE) {
		// ArrayList<double[]> returnList = new ArrayList<double[]>();
		int k = evalPoints.length;
		// double[] variance = new double[k];
		// double[] sqBias = new double[k];
		// double[] mse = new double[k];
		double[] tmp = new double[k];
		for (double[][] density : listDensity) {
			// returnList.add(computeMISE(dist, evalPoints, density, a, b, tmp, tmp, tmp));
			listMISE.add(computeMISE(dist, evalPoints, density, a, b, tmp, tmp, tmp));
		}
		// return returnList;
	}

	/**
	 * Gives a short description of the estimator.
	 * 
	 * @return a short description.
	 */
	public abstract String toString();

	/**
	 * Gives a plot of the estimated density. The \f$x\f$-values are passed in \a
	 * evalPoints and the \f$y\f$-values in \a density. The user may also set the
	 * title of the plot via \a plotTitle as well as the names of the axes via \a
	 * axisTitles. The latter contains the name of the \f$x\f$ axis as first element
	 * and the name of the \f$y\f$ axis as second.
	 * 
	 * The plot itself is returned as a string, which forms a stand-alone LaTex file
	 * (including necessary headers) implementing a tikZ picture.
	 * 
	 * This function merely tailors and simplifies the methods provided by \ref
	 * umontreal.ssj.stat.PgfDataTable for the purpose of plotting a univariate
	 * function. If the user seeks to produce more sophisticated plots, please refer
	 * to the aforementioned class.
	 * 
	 * @param evalPoints
	 *            the \f$x\f$-values.
	 * @param density
	 *            the \f$y\f$-values.
	 * @param plotTitle
	 *            the title of the plot.
	 * @param axisTitles
	 * @return
	 */
	public static String plotDensity(double[] evalPoints, double[] density, String plotTitle, String[] axisTitles) {
		double[][] plotData = { evalPoints, density };
		PgfDataTable table = new PgfDataTable(plotTitle, "", axisTitles, plotData);
		StringBuffer sb = new StringBuffer("");
		sb.append(PgfDataTable.pgfplotFileHeader());
		sb.append(table.drawPgfPlotSingleCurve(plotTitle, "axis", 0, 1, -1, "", ""));
		sb.append(PgfDataTable.pgfplotEndDocument());
		return sb.toString();
	}

	/**
	 * Computes the mean and the standard deviation of the observations of @f$m@f$
	 * simulations given in \a data.
	 * 
	 * @remark **Florian:** this should probably go somewhere else (DEDerivative as
	 *         Private?)
	 * 
	 * @param data
	 *            the observations.
	 * @return the mean and standard deviation.
	 */

	protected static double[] estimateMeanAndStdDeviation(double[][] data) {
		double[] result = new double[2];
		int m = data.length;
		int n = data[0].length;
		double stdDeviation = 0.0;
		double x, y;
		double meanSum = 0.0;
		for (int r = 0; r < m; r++) {

			double mean = 0.0;
			double var = 0.0;
			for (int i = 0; i < n; i++) {
				x = data[r][i];
				y = x - mean;
				mean += y / ((double) (i + 1.0));
				var += y * (x - mean);
			}
			stdDeviation += Math.sqrt(var / ((double) n - 1.0));
			meanSum += mean;
		}
		result[0] = meanSum / (double) m;
		result[1] = stdDeviation / (double) m;
		return result;
	}

	/**
	 * Computes the Coefficient of determination \f$R^2\f$ of the observed data \a
	 * data and the estimated data \a dataEstimated.
	 * 
	 * For observed data \f$y=(y_1,y_2,\dots,y_n)\f$ and estimated data \f$
	 * f=(f_1,f_2,\dots,f_n)\f$ this is defined as \f[ R^2 = 1 -
	 * \frac{\textrm{SS}_{\text{res}}}{\textrm{SS}_{\text{tot}}}, \f] where \f$
	 * \textrm{SS}_{\text{res}} \f$ denotes the sum of squares of the residuals \f[
	 * \textrm{SS}_{\text{res}} = \sum_{i=1}^n (f_i - y_i)^2 \f] and where
	 * \f$\textrm{SS}_{\text{tot}}\f$ is the total sum of squares \f[
	 * \textrm{SS}_{\text{tot}} = \sum_{i=1}^n (y_i - \bar{y})^2. \f] The closer
	 * this quantity is to one, the better the approximation of \f$y\f$ by \f$f\f$.
	 * 
	 * @remark **Florian:** this should probably go somewhere else (stat? we do not
	 *         need it for dens. est. per-se).
	 * 
	 * @param data
	 *            the observed data
	 * @param dataEstimated
	 *            the estimated data
	 * @return the Coefficient of determination \f$R^2\f$
	 */

	protected static double coefficientOfDetermination(double[] data, double[] dataEstimated) {
		int i;
		int max = data.length;
		double maxInv = 1.0 / (double) max;
		double dataMean = 0.0;
		double SSres = 0.0;
		double SStot = 0.0;
		for (i = 0; i < max; i++)
			dataMean += data[i];
		dataMean *= maxInv;
		for (i = 0; i < max; i++) {
			SSres += (data[i] - dataEstimated[i]) * (data[i] - dataEstimated[i]);
			SStot += (data[i] - dataMean) * (data[i] - dataMean);
		}
		return 1.0 - SSres / SStot;
	}
}

// ************ ******************************************************
// Manipulating the interval --> handled in each realization individually.
// For KDE it's not really important and for Histogram, e.g., you should
// not be able to change [a,b] unless you change everything.
// TODO: explain that [a,b] is more important for some estimators than for
// others;
// need it for IV, etc.
//
// /**left boundary of the interval over which we want to estimate */
// double a;
// /**right boundary of the interval over which we want to estimate */
// double b;
//
// /**
// * Gives the left boundary {@link #a} of the interval over which we estimate.
// *
// * @return the left boundary of the interval
// */
// public double geta() {
// return a;
// }
//
// /**
// * Gives the right boundary {@link #b} of the interval over which we estimate.
// *
// * @return the right boundary of the interval
// */
// public double getb() {
// return b;
// }
//
// /**
// * Sets the interval @f$[a,b]@f$ over which we estimate.
// *
// * @param a
// * left boundary of the interval.
// * @param b
// * right boundary of the interval.
// */
// public void setRange(double a, double b) {
// this.a = a;
// this.b = b;
// }
//

// **********************************************************************
// Should not be done this way. Shall be handled by the constructor!
// Reason: it does not do anything for the KDE, DEDirect, etc.
// /**
// * Constructs the estimator from the data points in vector \a data.
// *
// * @param data
// * the data points.
// */
// public abstract void constructDensity(double[] data);
//
// /**
// * Returns the value of the estimator evaluated at point \a x. It assumes that
// * the density has been constructed before.
// *
// * @param x
// * the point at which the density is to be evaluated.
// * @return the value of the estimated density at x.
// */

// ***********************************************************************
// This does not work without "constructDensity(...)" that way.
// /**
// * This method constructs \f$m\f$ independent realizations of a density
// * estimator and evaluates them at \a evalPoints. More precisely, the matrix
// \a
// * data contains the results of \f$m\f$ independent simulations of the
// * underlying model. From each such simulation result, this method constructs
// a
// * density estimator, evaluates it at \a evalPoints, and writes the resulting
// * vector to the corresponding row of \a density.
// *
// * Note that this method works essentially different than #evalDensity(double)
// * and {@link #evalDensity(double[], double[])} in that it constructs \f$m\f$
// m density
// * estimators and evaluates them, instead of considering only one estimator,
// * which needs to be constructed beforehand.
// * @remark **Florian:**
// *
// * @param evalPoints
// * the points at which the estimators shall be evaluated.
// * @param density
// * the values of the densities at the evaluation points.
// */
// public void evalDensity( double[] evalPoints, double[][] density) {
// int m = density.length;
// int numEvalPoints = evalPoints.length;
// for (int rep = 0; rep < m; rep++) {
// density[rep] = new double[numEvalPoints];
// evalDensity(evalPoints, density[rep]);
// }
// }
//
// /**
// * Same as #evalDensity(double[][], double[], double[][]) but using \a k
// * equidistant evaluation points over @f$[a,b]@f$ generated by
// * #equidistantPoints(int).
// *
// * @param k
// * the number of equidistant points at which the estimators shall be
// * evaluated.
// * @param density
// * the values of the densities at the evaluation points.
// */
// public void evalDensity(int k, double[][] density) {
// evalDensity(getEquidistantPoints(k), density);
// }
