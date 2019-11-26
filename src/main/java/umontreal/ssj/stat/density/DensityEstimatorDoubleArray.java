package umontreal.ssj.stat.density;

import java.util.ArrayList;

import umontreal.ssj.probdist.ContinuousDistribution;
import umontreal.ssj.stat.PgfDataTable;

/**
 * Same as @ref DensityEstimator but here the observations of the underlying
 * model are \f$t\f$-dimensional. This can be useful for conditional density
 * estimators or likelihood-ratio estimators, for instance. Here, mostly tools
 * to evaluate the density are implemented, where evaluation at one point is to
 * be defined individually for each realization and is thus kept abstract.
 * Additional functionalities for experiments such as plotting the density,
 * computing the integrated variance or the mean integrated square error, etc.,
 * can still be used from the @ref DensityEstimator class.
 * 
 * @author florian
 *
 */
public abstract class DensityEstimatorDoubleArray {

	/**
	 * The data associated with this DensityEstimatorDoubleArray object, if any.
	 */
	protected double[][] data;

	/**
	 * Sets the observations for the density estimator to \a data. Note that, in
	 * some cases, this might require to completely reconstruct the density
	 * estimator.
	 * 
	 * @param data the desired observations.
	 */
	public abstract void setData(double[][] data);

	/**
	 * Gives the observations for this density estimator, if any.
	 * 
	 * @return the observations for this density estimator.
	 */
	public double[][] getData() {
		return data;
	}

	/**
	 * Evaluates the density estimator at \a x.
	 * 
	 * @param x the evaluation point.
	 *
	 * @return the density estimator evaluated at \f$x\f$.
	 */

	public abstract double evalDensity(double x);

	/**
	 * Evaluates the density estimator at the points in \a evalPoints. By default,
	 * this method calls `evalDensity(double)` for each entry of \a evalPoints. Many
	 * density estimators can handle evaluation at a set of points more efficiently
	 * than that. If so, it is suggested to override this method in the
	 * implementation of the corresponding estimator.
	 *
	 * @return the density estimator evaluated at the points \a evalPoints.
	 */
	public double[] evalDensity(double[] evalPoints) {
		int k = evalPoints.length;
		double[] dens = new double[k];
		for (int j = 0; j < k; j++)
			dens[j] = evalDensity(evalPoints[j]);
		return dens;
	}

	/**
	 * Sets the observations for the density estimator to \a data and evaluates the
	 * density at each point in \a evalPoints.
	 * 
	 * @param evalPoints the evaluation points.
	 * @param data       the observations.
	 * @return the density estimator defined by \a data evaluated at each point in
	 *         \a evalPoints.
	 */
	public double[] evalDensity(double[] evalPoints, double[][] data) {
		setData(data);
		return evalDensity(evalPoints);
	}

	/**
	 * This method is particularly designed to evaluate the density estimator in
	 * such a way that the result can be easily used to estimate the empirical IV
	 * and other convergence-related quantities.
	 * 
	 * Assume that we have \f$m\f$ independent realizations of the underlying model.
	 * For each such realization this method constructs a density and evaluates it
	 * at the points from \a evalPoints. The independent realizations are passed via
	 * the 3-dimensional \f$m\times n \times t\f$ array \a data, where \f$n\f$
	 * denotes the number of observations per realization and \f$t\f$ the model
	 * dimension. Hence, its first index identifies the independent realization
	 * while its second index identifies a specific observation of this realization.
	 * 
	 * The result is returned as a \f$m\times k\f$ matrix, where \f$k \f$ is the
	 * number of evaluation points, i.e., the length of \a evalPoints. The first
	 * index, again, identifies the independent realization whereas the second index
	 * corresponds to the point of \a evalPoints at which the density estimator was
	 * evaluated.
	 * 
	 *
	 * @param evalPoints the evaluation points.
	 * @param data       the three-dimensional array carrying the observations of
	 *                   \f$m\f$ independent realizations of the underlying model.
	 *
	 * @return the density estimator for each realization evaluated at \a
	 *         evalPoints.
	 */
	public double[][] evalDensity(double[] evalPoints, double[][][] data) {
		int m = data.length;
		double[][] density = new double[m][];
		for (int r = 0; r < m; r++)
			density[r] = evalDensity(evalPoints, data[r]);

		return density;
	}

	/**
	 * This function is particularly designed for experiments with many different
	 * types of density estimators, as it evaluates all of these estimators at the
	 * points in \a evalPoints. To this end, the user passes a list of density
	 * estimators in \a listDE as well as \f$m\f$ independent realizations of the
	 * underlying model consisting of \f$n\f$ observations each in the \f$m\times n
	 * \times t\f$ array \a data.
	 * 
	 * This method then calls `evalDensity(double[], double[][][])` for each density
	 * estimator in \a listDE, thus evaluating the respective density estimator at
	 * the \f$k\f$ points in \a evalPoints and adds the resulting \f$m\times k\f$
	 * array to \a listDensity.
	 * 
	 * @param listDE      the list of density estimators.
	 * @param evalPoints  the evaluation points.
	 * @param data        the three-dimensional array carrying the observations of
	 *                    \f$m\f$ independent realizations of the underlying model.
	 * @param listDensity a list to which the evaluations at \a evalPoints of each
	 *                    density estimator in \a listDE are added.
	 * 
	 */
	public static void evalDensity(ArrayList<DensityEstimatorDoubleArray> listDE, double[] evalPoints,
			double[][][] data, ArrayList<double[][]> listDensity) {
		for (DensityEstimatorDoubleArray de : listDE)
			listDensity.add(de.evalDensity(evalPoints, data));
	}

	/**
	 * Gives a short description of the estimator.
	 * 
	 * @return a short description.
	 */
	public abstract String toString();

}
