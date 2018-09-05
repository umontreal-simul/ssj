
package umontreal.ssj.stat.density;

import java.util.ArrayList;

import umontreal.ssj.probdist.ContinuousDistribution;

/**
 * This abstract class implements a univariate density estimator (DE). To this
 * end, it provides basic tools to evaluate the DE at one point \f$x\f$ or at an
 * array of points \f$\{x_1, x_2, \dots, x_k\} \f$. More precisely, the single
 * point evaluation #evalDensity(double) is abstract, since it will most likely
 * differ between realizations. For the evaluation on a set of points,
 * #evalDensity(double[], double[]), a default implementation is provided, but
 * very often specific estimators will have more efficient evaluation
 * algorithms. Furthermore, this class includes a method to plot the estimated
 * density.
 * 
 * There are also several more involved methods covered by this class, most of
 * which are concerned with the convergence behavior of DEs. Nevertheless, they
 * can be useful in many other cases beyond convergence behavior. As these
 * usually require more than one realization of a DE, they are implemented as
 * static methods. For instance, #evalDensity(ArrayList, double[]) allows to
 * evaluate several DEs passed in a list at the same evaluation points.
 * 
 * For measuring the performance of a DE, we need to confine ourselves to
 * estimation over a finite interval \f$[a,b]\f$. One standard way to assess the
 * quality of a DE is via the mean integrated square error (MISE). It can be
 * further decomposed into the integrated variance (IV), and the integrated
 * square bias (ISB) \f[ \textrm{MISE} = \int_a^b\mathbb{E} [\hat{f} -
 * f(x)]^2\mathrm{d}x = \int_a^b\textrm{Var}[\hat{f}(x)], \f] where \f$f\f$
 * denotes the true density and \f$\hat{f}\f$ the DE. The result of this can
 * subsequently be used to compute the integrated
 * 
 * 
 * over a finite interval \f$[a,b]\f$. The density is constructed from a set of
 * \f$n\f$ observations. Some of these observations can lie outside the interval
 * \f$[a,b]\f$.
 * 
 * The density can be evaluated at a single point \f$x\f$ or on a grid of
 * evaluation points.
 * 
 * Furthermore, this class implements basic routines for the computation of the
 * integrated variance (IV) as well as for calculating basic error measures,
 * such as the integrated square bias (ISB) and the mean integrated square error
 * (MISE) in case the true density is known.
 */
public abstract class DensityEstimator {
	
	/**
	 * Constructs a density estimator over the interval \f$[a,b]\f$ based on the
	 * observations \a data, if necessary, and evaluates it at \a x.
	 * @param x the evaluation point.
	 * @param data the observations for constructing the density estimator.
	 * @param a the left boundary of the interval.
	 * @param b the right boundary of the interval
	 * @return the density estimator evaluated at \f$x\f$.
	 */
	public abstract double evalDensity(double x, double data[], double a, double b);

	/**
	 * Constructs a density estimator over the interval \f$[a,b]\f$ based on the
	 * observations \a data, if necessary, and evaluates it at the points in \a x.
	 * @param x the evaluation points.
	 * 
	 * By default, this method calls #evalDensity(double, double[], double, double)
	 * for each entry of \a x. Many density estimators can handle evaluation
	 * at a vector of points more efficiently than that. If so, it is suggested to
	 * override this method in the implementation of the corresponding estimator.
	 * Note that the construction of the density estimator -- if applicable --
	 * needs to be handled by this method as well.
	 * 
	 * @param data the observations for constructing the density estimator.
	 * @param a the left boundary of the interval.
	 * @param b the right boundary of the interval
	 * @return the density estimator evaluated at the points \a x.
	 */
	public  double [] evalDensity(double[] x, double[] data, double a, double b) {
		int k = x.length;
		double[] dens = new double[k];
		for(int j = 0; j < k; j++)
			dens[j] = evalDensity(x[j],data,a,b);
		return dens;
	}

	/**
	 * This method is particularly designed to evaluate the density estimator in such a way that the result can be easily used to
	 * estimate the empirical IV and other convergence-related quantities.
	 * 
	 * Assume that we have \f$m\f$ independent realizations of the underlying model. For each such realization this method
	 * constructs a density estimator over \f$[a,b]\f$ and evaluates it at the points from \a x. The independent realizations
	 * are passed via the 2-dimensional \f$m\times n\f$array \a data, where \f$n\f$ denotes the number of observations per realization. 
	 * Hence, its first index identifies the independent realization while its second
	 * index identifies a specific observation of this realization.
	 * 
	 * The result is returned as a \f$m\times k\f$ matrix, where \f$k \f$ is the number of evaluation points, i.e., the length of \a x.
	 * The first index, again, identifies the independent realization whereas the second index indicates at which point of \a x the
	 * density estimator was evaluated.
	 * 
	 *
	 * @param x the evaluation points.
	 * @param data the two dimensional array carrying the observations of \f$m\f$ independent realizations of the underlying model.
	 * @param a the left boundary of the interval.
	 * @param b the right boundary of the interval.
	 * @return the density estimator for each realization evaluated at \a x.
	 */
	public double [][] evalDensity(double[] x, double[][] data, double a, double b){
		int m = data.length;
		double [][] dens = new double[m][];
		for(int r = 0; r < m; r++) 
			dens[r] = evalDensity(x,data[r],a,b);
		return dens;
	}
	
	/**
	 * This method computes the empirical variance based on the values given in
	 * \a data. More precisely, \a density is a \f$m\times k\f$ matrix, whose entries
	 * correspond to \f$m\f$ independent realizations of the density estimator, each evaluated
	 * at \f$k\f$ evaluation points. Such a matrix can, for instance, be obtained by
	 * #evalDensity(double[], double[][], double, double).
	 * 
	 * The empirical variance is computed at each of those \f$k \f$ evaluation points
	 * and returned in an array of size \f$k\f$.
	 * 
	 * @param density the estimated density of \f$m\f$ independent realizations of the 
	 * estimator, each evaluated at \f$k\f$ evaluation points. 
	 * @return the empirical variance at those \f$k\f$ evaluation points.
	 */
	public static double [] computeVariance(double[][] density) {
		int m = density.length; //number of indep. replications
		int numEvalPoints = density[0].length; //number of evaluation points

		double x, y;
	
		double meanDens[] = new double[numEvalPoints]; // Average value over the rep replicates
		double varDens[] = new double[numEvalPoints]; // Variance at each evaluation point
		// Arrays.fill(meanDens, 0.0);
		// Arrays.fill(varDens, 0.0);
		for (int r = 0; r < m; r++) {
			// Update the empirical mean and variance at each evaluation point.
			for (int j = 0; j < numEvalPoints; j++) {
				x = density[r][j];
				y = x - meanDens[j];
				meanDens[j] += y / (double) (r + 1);
				varDens[j] += y * (x - meanDens[j]);
			}
		}
		
		for(int j = 0; j < numEvalPoints; j++) //normalize
			varDens[j] /= (double)(m - 1);
		
		return varDens;
	}
	
	/**
	 * This method estimates the empirical IV over the interval \f$[a,b]\f$. Based on the density estimates of \f$m\f$ independent
	 * replications of the density estimator evaluated at \f$k\f$ evaluation points, which is provided by \a density, it computes
	 * the empirical variance at each evaluation point and stores it in \a variance.
	 * 
	 * To estimate the empirical IV, we sum up the variance at the \f$k\f$ evaluation points and multiply by \f$(b-a)/k\f$, i.e. we
	 * approximate the empirical IV by an equally weighted quadrature rule with aforementioned evaluation points as integration nodes.
	 * Note that this is only an approximation of the true empirical IV and that the approximation quality significantly depends on
	 * the choice of evaluation points.
	 * 
	 * The data for the variance are given in the two-dimensional \f$m\times k\f$ array \a density, which is also described in 
	 *  #computeVariance(double[][]) and can be obtained by  #evalDensity(double[], double[][], double, double). The boundaries of
	 *  the interval are given by \a a and \a b. Note that the array \a variance needs to be of length \f$k\f$.
	 * 
	 * 
	 * @param density the \f$m\times k\f$ array that contains the data of evaluating \f$m\f$ replicates of the density estimator
	 * at \f$k\f$ evaluation points
	 * @param a the left boundary of the interval.
	 * @param b the right boundary of the interval.
	 * @param variance the array of length \f$k\f$ in which the variance at each evaluation point is stored.
	 * @return the estimated empirical IV over \f$[a,b]\f$.
	 */
	public static double computeIV(double[][] density, double a, double b, double[] variance) {
		variance = computeVariance(density);
		int k = density[0].length;
		double iv = 0.0;
		for(double var : variance)
			iv += var;
		return iv * (b - a)/(double) k;
	}
	
	public static double[] computeMISE(ContinuousDistribution dist, double[][] density, double a, double b, doubl)
	
	/**
	 * Gives a short description of the estimator.
	 * 
	 * @return a short description.
	 */
	public abstract String toString();
	

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

	public String plotDensity(double[] evalPoints, double[] density) {
		return "Still todo";
	}

	/**
	 * TODO: check what is faster: void & passing or create list each call.
	 * TODO: add data[][]
	 * extension of #evalDensity(double[], double[]) to lists of DE's
	 * 
	 * @param listDE
	 * @param evalPoints
	 * @return
	 */
	public static ArrayList<double[]> evalDensity (ArrayList<DensityEstimator> listDE, double[] evalPoints) {
		ArrayList<double[]> listEvalDens = new ArrayList<double[]>();
		for (DensityEstimator de : listDE) {
			double[] density = new double[evalPoints.length];
			de.evalDensity(evalPoints, density);
			listEvalDens.add(density);
		}
		return listEvalDens;
	}

	
	/**
	 * TODO: return variance at each point
	 * TODO: make clear that it's only an estimate
	 * This method estimates the empirical IV based on \f$m\f$ realizations of the
	 * density estimator which have previously been evaluated at \f$k\f$
	 * predetermined points. The matrix \a estDensities has dimensions \f$m\times k
	 * \f$, i.e. each row contains the evaluations of one density estimator.
	 * 
	 * @param estDensities
	 *            the \f$m\times k\f$matrix containing the results of estimating
	 *            \f$m\f$ densities at \f$k\f$ evaluation points.
	 * @return the estimated empirical IV
	 */
	public static double computeIV(ArrayList<double[]> listEvalDensities, double a, double b) {
		int m = listEvalDensities.size();
		int k = listEvalDensities.get(0).length;

		double x, y;
		// int numEvalPoints = evalPoints.length();
		double meanDens[] = new double[k]; // Average value over the rep replicates
		double varDens[] = new double[k]; // Variance at each evaluation point
		// Arrays.fill(meanDens, 0.0);
		// Arrays.fill(varDens, 0.0);
		for (int rep = 0; rep < m; rep++) {
			// Estimate the density for this rep and evaluate it at the evaluation points
			// Update the empirical mean and sum of squares of centered observations at each
			// evaluation point.
			for (int j = 0; j < k; j++) {
				x = (listEvalDensities.get(rep))[j];
				y = x - meanDens[j];
				meanDens[j] += y / (double) (rep + 1);
				varDens[j] += y * (x - meanDens[j]);
			}
		}
		double sumVar = 0.0;
		for (int i = 0; i < k; ++i)
			sumVar += varDens[i];
		return sumVar * (b - a) / (double) (k * (m - 1)); // Empirical integrated variance.
	}

	/**
	 * For a \ref umontreal.ssj.mcqmctools.MonteCarloModelDensityKnown this method
	 * estimates the empirical ISB based on \f$m\f$ realizations of the density
	 * estimator which have previously been evaluated at the \f$k\f$ points stored
	 * in \a evalPoints. The matrix \a estDensities has dimensions \f$m\times k \f$,
	 * i.e. each row contains the evaluations of one density estimator.
	 * 
	 * This is particularly useful for testing the density estimator with a toy
	 * model.
	 * 
	 * @param model
	 *            the underlying model.
	 * @param estDensities
	 *            the \f$m\times k\f$matrix containing the results of evaluating
	 *            \f$m\f$ densities at \f$k\f$ evaluation points.
	 * @param evalPoints
	 *            the \f$k\f$ evaluation points at which the \f$m\f$ densities were
	 *            evaluated.
	 * @return the empirical ISB
	 */
	public static double computeISB(ContinuousDistribution dist, ArrayList<double[]> listEvalDensities,
			double[] evalPoints, double a, double b) {

		int k = evalPoints.length;
		int m = listEvalDensities.size();

		double meanDens[] = new double[k];
		double sqBiasDens[] = new double[k]; // squared bias

		double x, y, z;

		for (int rep = 0; rep < m; rep++) {
			// Update the empirical mean and sum of squares of centered
			// observations at each evaluation point.
			for (int j = 0; j < k; j++) {
				x = (listEvalDensities.get(rep))[j];
				y = x - meanDens[j];
				meanDens[j] += y / (double) (rep + 1);
			}

		}

		for (int j = 0; j < k; j++) {
			z = (meanDens[j] - dist.density(evalPoints[j]));
			sqBiasDens[j] = z * z;
		}

		double isb = 0.0;
		for (int i = 0; i < k; ++i)
			isb += sqBiasDens[i];

		return isb * (b - a) / ((double) k);

	}

	public static double computeISB(ContinuousDistribution dist, ArrayList<double[]> listEvalDensities,
			double[] evalPoints) {
		return computeISB(dist, listEvalDensities, evalPoints, evalPoints[0], evalPoints[evalPoints.length - 1]);
	}

	/**
	 * For a \ref MonteCarloModelDensityKnown this method estimates the empirical
	 * MISE based on \f$m\f$ realizations of the density estimator which have
	 * previously been evaluated at the \f$k\f$ points stored in \a evalPoints. The
	 * matrix \a estDensities has dimensions \f$m\times k \f$, i.e. each row
	 * contains the evaluations of one density estimator.
	 * 
	 * This is particularly useful for testing the density estimator with a toy
	 * model.
	 * 
	 * @param model
	 *            the underlying model.
	 * @param estDensities
	 *            the \f$m\times k\f$matrix containing the results of evaluating
	 *            \f$m\f$ densities at \f$k\f$ evaluation points.
	 * @param evalPoints
	 *            the \f$k\f$ evaluation points at which the \f$m\f$ densities were
	 *            evaluated.
	 * @return the empirical MISE
	 */

	public static double computeMISE(ContinuousDistribution dist, ArrayList<double[]> listEvalDensities,
			double[] evalPoints, double a, double b) {

		int m = listEvalDensities.size();
		int k = evalPoints.length;

		double mseDens[] = new double[k]; // mse at evalPoints

		double y;
		for (int rep = 0; rep < m; rep++) {

			for (int j = 0; j < k; j++) {
				y = listEvalDensities.get(rep)[j] - dist.density(evalPoints[j]); // model.density(x);
				mseDens[j] += y * y;
			}
		}

		double sumMISE = 0.0;
		for (int i = 0; i < k; ++i)
			sumMISE += mseDens[i];
		double fact = (b - a) / ((double) (k * m));
		return sumMISE * fact;
	}

	public static double computeMISE(ContinuousDistribution dist, ArrayList<double[]> listEvalDensities,
			double[] evalPoints) {
		return computeMISE(dist, listEvalDensities, evalPoints, evalPoints[0], evalPoints[evalPoints.length - 1]);
	}

	/**
	 * TODO: merge with computeMISE.
	 * For a \ref MonteCarloModelDensityKnown this method estimates the empirical
	 * IV, ISB, and MISE based on \f$m\f$ realizations of the density estimator,
	 * which have previously been evaluated at the \f$k\f$ points stored in \a
	 * evalPoints, and returns them in an array in this order. The matrix \a
	 * estDensities has dimensions \f$m\times k \f$, i.e. each row contains the
	 * evaluations of one density estimator.
	 * 
	 * Estimating the empirical IV, ISB, and MISE in parallel saves a few redundant
	 * computations that would arise if they were estimated individually.
	 * 
	 * This method is particularly useful for testing the density estimator with a
	 * toy model.
	 * 
	 * @param model
	 *            the underlying model.
	 * @param estDensities
	 *            the \f$m\times k\f$matrix containing the results of evaluating
	 *            \f$m\f$ densities at \f$k\f$ evaluation points.
	 * @param evalPoints
	 *            the \f$k\f$ evaluation points at which the \f$m\f$ densities were
	 *            evaluated.
	 * @return the estimated empirical IV, ISB, and MISE.
	 */
	public static double[] computeIVandISBandMISE(ContinuousDistribution dist, ArrayList<double[]> listEvalDensities,
			double[] evalPoints, double a, double b) {

		double x, y, z;
		double dens;
		int k = evalPoints.length;
		int m = listEvalDensities.size();

		double meanDens[] = new double[k]; // Average value over m
											// replicates
		double varDens[] = new double[k]; // Variance at each
											// evaluation point
		double mseDens[] = new double[k]; // MSE at each evaluation
											// point

		for (int rep = 0; rep < m; rep++) {
			// Update the empirical mean, sum of squares, and mse of
			// observations at each evaluation point.
			for (int j = 0; j < k; j++) {
				x = (listEvalDensities.get(rep))[j];
				y = x - meanDens[j];
				dens = dist.density(evalPoints[j]);
				z = x - dens;

				meanDens[j] += y / (double) (rep + 1);
				varDens[j] += y * (x - meanDens[j]);
				mseDens[j] += z * z;
			}

		}

		double iv = 0.0;
		double mise = 0.0;
		// double isb = 0.0;
		for (int i = 0; i < k; ++i) {
			iv += varDens[i];
			mise += mseDens[i];
			// isb += biasDens[i];
		}

		double fact = (b - a) / ((double) (k * (m - 1.0)));
		iv *= fact;
		mise *= (b - a) / ((double) (k * m));

		double[] res = { iv, mise - iv, iv };
		return res;

	}

	public static double[] computeIVandISBandMISE(ContinuousDistribution dist, ArrayList<double[]> listEvalDensities,
			double[] evalPoints) {
		return computeIVandISBandMISE(dist, listEvalDensities, evalPoints, evalPoints[0],
				evalPoints[evalPoints.length - 1]);
	}

	// *********************************************************************************************
	// This is not needed for now. If we decide do keep it, make it private?
	// Would also need to rethink [a,b] if we kept it.
	// /**
	// * Generates \a numPoints equidistant points over \f$[a,b]\f$ by fixing the
	// * distance between two points as \f$\delta = (b-a)/k\f$ and set the first
	// point
	// * as \f$a + \delta/2\f$.
	// *
	// * @param k
	// * the number of points to be returned.
	// * @return an array of equidistant points over \f$[a,b]\f$.
	// */
	// private double[] getEquidistantPoints(int k) {
	// double evalPoints[] = new double[k];
	// double delta = (b - a) / (k);
	// for (int j = 0; j < k; j++)
	// evalPoints[j] = a + delta * (0.5 + j);
	//
	// return evalPoints;
	// }

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
