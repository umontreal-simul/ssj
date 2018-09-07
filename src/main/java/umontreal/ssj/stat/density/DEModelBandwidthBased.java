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

import umontreal.ssj.mcqmctools.MonteCarloModelDensityKnown;
import umontreal.ssj.mcqmctools.MonteCarloModelDouble;

/** Recall
 * that asymptotically, in both cases, \f$\textrm{ISB} \approx B h^{\alpha}\f$.
 * Now, for MC and hence for RQMC too, it is known (see \cite tSCO15a) that \f$
 * B = R(f')/12 \f$ for histograms and \f$B = \mu_2(k)R(f'')/4 \f$ for KDEs,
 * where 
 * 
 * \f[ R(g) = \int_a^b g^2(x)\mathrm{d}\x\f] 
 * 
 * denotes the roughness
 * functional (see
 * umontreal.ssj.stat.DensityEstimator#roughnessFunctional(double[], double,
 * double) ) and where \f$\mu_k \f$ denotes the \f$k\f$th moment
 * \f$\int_{-\infty}^{\infty}x^k g(x)\mathrm{d}x \f$. Thus, in order to estimate
 * \f$B\f$ one needs to estimate derivatives of the sought density.
 */
public class DEModelBandwidthBased implements DensityEstimationModel{

	@Override
	public double estimateIV(DensityEstimator de) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double estimateISB(DensityEstimator de) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double estimateMISE(DensityEstimator de) {
		// TODO Auto-generated method stub
		return 0;
	}
	
	/**
	 * Depending on whether \a maxDerivative is even or odd, this method computes
	 * \f$ h_{\text{AMISE}}^{(r_0)}\f$ for \f$r_0=2\f$ or \f$r_0=1\f$, respectively.
	 * This is done by recursively applying the formula 
	 * ( {@link REF_stat_density_DensityDerivativeEstimator_hopt hopt} )
	 *  from above. The recursion
	 * starts at @f$r=\textit{maxDerivative} - 2@f$ with the initial value \a init
	 * for @f$R(f^{(r+2)})@f$. The second moment of \f$K\f$ is given by \a mu2 and
	 * the second moments of the derivatives \f$k^{(r)},k^{(r-2)},\dots,k^{(r_0)}
	 * \f$ are passed via \a mu2Derivatives. It is assumed that the elements of \a
	 * mu2Derivatives are sorted by order of derivative in a decreasing manner, the
	 * first one being \f$\mu_2(k^{(r)})\f$.
	 * 
	 * @param dde            the DDE to use in the recursion.
	 * @param data           observations from one simulation.
	 * @param mu2            the second moment of \f$k\f$.
	 * @param mu2Derivatives the second moments of the derivatives of \f$k\f$.
	 * @param evalPoints     the integration nodes.
	 * @param maxDerivative  the highest order of the density derivative considered.
	 * @param init           initial value of the roughness functional in the
	 *                       denominator.
	 * @return either \f$ h_{\text{AMISE}}^{(2)}\f$ or \f$
	 *         h_{\text{AMISE}}^{(1)}\f$, depending on whether \a maxDerivative is
	 *         even or odd.
	 */
	public static double hOptAmise(DensityDerivativeEstimator dde, double[] data, double mu2, double[] mu2Derivatives,
			double[] evalPoints, double a, double b, int maxDerivative, double init) {

		int n = data.length;
		int order = maxDerivative - 2;
		// Index of the right derivative in mu2Derivatives
		int derIndex = 0;

		int numEvalPoints = evalPoints.length;

		double[] estDensDerivative = new double[numEvalPoints];
		// Arrays.fill(estDensDerivative, 0.0);

		double h = hAmiseR(order, mu2, mu2Derivatives[derIndex], init, n);
		// System.out.println("hAmiseGaussian = " + h);
		order -= 2;
		derIndex++;
		double rF = init; //the roughness functional
		// only enters loop when the last computed h is not already h^(2) or
		// h^(1)
		while (order > 0) {

			dde.setOrder(order + 2);
			dde.setH(h);
			
			estDensDerivative = dde.evalDensity(evalPoints, data,a,b);
			rF = roughnessFunctional(estDensDerivative,a,b);

			h = hAmiseR(order, mu2, mu2Derivatives[derIndex], rF, n);
			order -= 2;
			derIndex++;

		}

		return h;
	}

	

	/**
	 * Computes the asymptotically optimal multiplicative constant in the ISB for a
	 * histogram as described above. The bandwidth in these density derivate
	 * estimates is computed recursively via
	 * {@link #hOptAmise(DensityDerivativeEstimator, double[], double, double[], double[], int, double)}
	 * starting with initial value \a init for the roughness functional of the
	 * derivative of the density of order \a maxDerivative. Each occuring
	 * integration, e.g. for evaluating a roughness-functional, is carried out by a
	 * quadrature rule using the integration nodes \a evalPoints.
	 * 
	 * Analogously to
	 * {@link #hOptAmise(DensityDerivativeEstimator, double[], double, double[], double[], int, double)}
	 * \a mu2 is the second moment of the kernel and \a mu2Derivatives contains the
	 * second moments of the derivatives of \f$k\f$ decreasing w.r.t. the order of
	 * the derivatives, starting with order \a maxDerivative - 2.
	 * 
	 * Note that for a histogram estimator \a maxDerivative has to be odd!
	 * 
	 * @param dde            the DDE we use for determining the optimal bandwidth.
	 * @param data           the observations gained from @f$m@f$ independent
	 *                       simulations.
	 * @param mu2            the second moment of \f$k\f$.
	 * @param mu2Derivatives the second moments of the derivatives of \f$k\f$.
	 * @param evalPoints     the integration nodes.
	 * @param maxDerivative  the highest order of the density derivative considered.
	 * @param init           initial value of the roughness functional in the
	 *                       denominator.
	 * @param de             the histogram estimator for the sought density.
	 * @return the asymptotically optimal multiplicative constant in the ISB.
	 * 
	 * @remark **Florian:** Maybe we should throw an exception, when \a maxDerivative is not odd
	 */
	// TODO: add exception if maxDerivative is odd!
	public static double computeB(DensityDerivativeEstimator dde, double[][] data, double mu2, double[] mu2Derivatives,
			double[] evalPoints, int maxDerivative, double init, DEHistogram de) {
		int m = data.length;
		double B = 0.0;
		double a = de.geta();
		double b = de.getb();

		int numEvalPoints = evalPoints.length;
		double[] estDensDerivative = new double[numEvalPoints];
//		Arrays.fill(estDensDerivative, 0.0);

		for (int r = 0; r < m; r++) {
			double h = hOptAmise(dde, data[r], mu2, mu2Derivatives, evalPoints, maxDerivative, init);
			// System.out.println("h = " + h);
			dde.setOrder(1);
			dde.setH(h);
			dde.constructDensity(data[r]);
			dde.evalDensity(evalPoints, estDensDerivative);
			double roughnessFunctional = 0.0;
			for (int i = 0; i < numEvalPoints; i++)
				roughnessFunctional += estDensDerivative[i] * estDensDerivative[i];
			roughnessFunctional *= (b - a) / (double) numEvalPoints;

			B += roughnessFunctional;
		}

		return B / (12.0 * (double) m);

	}

	/**
	 * Same as
	 * {@link #computeB(DensityDerivativeEstimator, double[][], double, double[], double[], int, double, DEHistogram)}
	 * but with \a numEvalPoints equidistant integration nodes.
	 * 
	 * @param dde            the DDE we use for determining the optimal bandwidth.
	 * @param data           the observations gained from @f$m@f$ independent
	 *                       simulations.
	 * @param mu2            the second moment of \f$k\f$.
	 * @param mu2Derivatives the second moments of the derivatives of \f$k\f$.
	 * @param numEvalPoints  the number of equidistant integration nodes.
	 * @param maxDerivative  the highest order of the density derivative considered.
	 * @param init           initial value of the roughness functional in the
	 *                       denominator.
	 * @param de             the histogram estimator for the sought density.
	 * @return the asymptotically optimal multiplicative constant in the ISB.
	 */

	public static double computeB(DensityDerivativeEstimator dde, double[][] data, double mu2, double[] mu2Derivatives,
			int numEvalPoints, int maxDerivative, double init, DEHistogram de) {
		return computeB(dde, data, mu2, mu2Derivatives, dde.getEquidistantPoints(numEvalPoints), maxDerivative, init,
				de);
	}

	/**
	 * Computes the asymptotically optimal multiplicative constant in the ISB for a
	 * KDE as described above. The bandwidth in these density derivate estimates is
	 * computed recursively via
	 * {@link #hOptAmise(DensityDerivativeEstimator, double[], double, double[], double[], int, double)}
	 * starting with initial value \a init for the roughness functional of the
	 * derivative of the density of order \a maxDerivative. Each occuring
	 * integration, e.g. for evaluating a roughness-functional, is carried out by a
	 * quadrature rule using the integration nodes \a evalPoints.
	 * 
	 * Analogously to
	 * {@link #hOptAmise(DensityDerivativeEstimator, double[], double, double[], double[], int, double)}
	 * \a mu2 is the second moment of the kernel and \a mu2Derivatives contains the
	 * second moments of the derivatives of \f$k\f$ decreasing w.r.t. the order of
	 * the derivatives, starting with order \a maxDerivative - 2.
	 * 
	 * Note that for a KDE \a maxDerivative has to be even!
	 * 
	 * @param dde            the DDE we use for determining the optimal bandwidth.
	 * @param data           the observations gained from @f$m@f$ independent
	 *                       simulations.
	 * @param                * @param mu2 the second moment of \f$k\f$.
	 * @param mu2Derivatives the second moments of the derivatives of \f$k\f$.
	 * @param evalPoints     the integration nodes.
	 * @param maxDerivative  the highest order of the density derivative considered.
	 * @param init           initial value of the roughness functional in the
	 *                       denominator.
	 * @param de             the KDE for the sought density.
	 * @return the asymptotically optimal multiplicative constant in the ISB.
	 */
	// TODO: add exception that maxDerivative has to be even!
	public static double computeB(DensityDerivativeEstimator dde, double[][] data, double mu2, double[] mu2Derivatives,
			double[] evalPoints, int maxDerivative, double init, DEKernelDensity de) {
		int m = data.length;
		double B = 0.0;
		double a = de.geta();
		double b = de.getb();

		int numEvalPoints = evalPoints.length;
		double[] estDensDerivative = new double[numEvalPoints];
//		Arrays.fill(estDensDerivative, 0.0);

		for (int r = 0; r < m; r++) {
			double h = hOptAmise(dde, data[r], mu2, mu2Derivatives, evalPoints, maxDerivative, init);
			// System.out.println("h = " + h);
			dde.setOrder(2);
			dde.setH(h);
			dde.constructDensity(data[r]);
			dde.evalDensity(evalPoints, estDensDerivative);
			double roughnessFunctional = 0.0;
			for (int i = 0; i < numEvalPoints; i++)
				roughnessFunctional += estDensDerivative[i] * estDensDerivative[i];
			roughnessFunctional *= (b - a) / (double) numEvalPoints;

			B += roughnessFunctional;
		}

		return 0.25 * mu2 * B / (double) m;

	}

	/**
	 * Same as
	 * {@link #computeB(DensityDerivativeEstimator, double[][], double, double[], double[], int, double, DEKernelDensity)}
	 * but with \a numEvalPoints equidistant integration nodes.
	 * 
	 * @param dde            the DDE we use for determining the optimal bandwidth.
	 * @param data           the observations gained from @f$m@f$ independent
	 *                       simulations.
	 * @param mu2            the second moment of \f$k\f$.
	 * @param mu2Derivatives the second moments of the derivatives of \f$k\f$.
	 * @param numEvalPoints  the number of equidistant integration nodes.
	 * @param maxDerivative  the highest order of the density derivative considered.
	 * @param init           initial value of the roughness functional in the
	 *                       denominator.
	 * @param de             the KDE for the sought density.
	 * @return the asymptotically optimal multiplicative constant in the ISB.
	 */

	public static double computeB(DensityDerivativeEstimator dde, double[][] data, double mu2, double[] mu2Derivatives,
			int numEvalPoints, int maxDerivative, double init, DEKernelDensity de) {
		return computeB(dde, data, mu2, mu2Derivatives, dde.getEquidistantPoints(numEvalPoints), maxDerivative, init,
				de);
	}
	
	/**
	 * Computes the asymptotically optimal multiplicative constant in the ISB for a
	 * histogram. This works exactly as
	 * {@link #computeB(DensityDerivativeEstimator, double[][], double, double[], double[], int, double, DEHistogram)}
	 * but with the kernel of the DDE being a Gaussian.
	 * 
	 * In this case, the second moment of the kernel is 1 and those of the
	 * derivatives can be obtained via {@link #densFunctionalGaussian(int, double)}.
	 * 
	 * Note that for a histogram estimator \a maxDerivative has to be odd!
	 * 
	 * @param dde
	 *            the DDE we use for determining the optimal bandwidth.
	 * @param data
	 *            the observations gained from @f$m@f$ independent simulations.
	 * @param evalPoints
	 *            the integration nodes.
	 * @param maxDerivative
	 *            the highest order of the density derivative considered.
	 * @param init
	 *            initial value of the roughness functional in the denominator.
	 * @param de
	 *            the histogram estimator for the sought density.
	 * @return the asymptotically optimal multiplicative constant in the ISB.
	 */
	public static double computeB(DDEGaussian dde, double[][] data, double[] evalPoints, int maxDerivative, double init,
			DEHistogram de) {
		int t = (maxDerivative - 1) / 2; // maxDerivative = 2t + 1
		double[] mu2Derivatives = new double[t - 1];
		for (int r = t - 1; r >= 1; r--)
			// the first order is 2t-1, this is init.
			// the first in mu2Derivatives is 2t-3, the second 2t-5,..., the last 2*1-1.
			mu2Derivatives[t - 1 - r] = densityFunctionalGaussian(2 * r - 1, 1.0);
		return computeB(dde, data, 1.0, mu2Derivatives, evalPoints, maxDerivative, init, de);

	}

	/**
	 * Same as
	 * {@link #computeB(DDEGaussian, double[][], double[], int, double, DEHistogram)}
	 * but with \a numEvalPoints equidistant integration nodes.
	 * 
	 * @param dde
	 *            the DDE we use for determining the optimal bandwidth.
	 * @param data
	 *            the observations gained from @f$m@f$ independent simulations.
	 * @param numEvalPoints
	 *            the number of equidistant integration nodes.
	 * @param maxDerivative
	 *            the highest order of the density derivative considered.
	 * @param init
	 *            initial value of the roughness functional in the denominator.
	 * @param de
	 *            the histogram estimator for the sought density.
	 * @return the asymptotically optimal multiplicative constant in the ISB.
	 */
	public static double computeB(DDEGaussian dde, double[][] data, int numEvalPoints, int maxDerivative, double init,
			DEHistogram de) {
		return computeB(dde, data, dde.getEquidistantPoints(numEvalPoints), maxDerivative, init, de);
	}

	/**
	 * Computes the asymptotically optimal multiplicative constant in the ISB for a
	 * KDE. This works exactly as
	 * {@link #computeB(DensityDerivativeEstimator, double[][], double, double[], double[], int, double, DEKernelDensity)}
	 * but with the kernel of the DDE being a Gaussian.
	 * 
	 * In this case, the second moment of the kernel is 1 and those of the
	 * derivatives can be obtained via {@link #densFunctionalGaussian(int, double)}.
	 * 
	 * Note that for a histogram estimator \a maxDerivative has to be odd!
	 * 
	 * @param dde
	 *            the DDE we use for determining the optimal bandwidth.
	 * @param data
	 *            the observations gained from @f$m@f$ independent simulations.
	 * @param evalPoints
	 *            the integration nodes.
	 * @param maxDerivative
	 *            the highest order of the density derivative considered.
	 * @param init
	 *            initial value of the roughness functional in the denominator.
	 * @param de
	 *            the KDE for the sought density.
	 * @return the asymptotically optimal multiplicative constant in the ISB.
	 */
	public static double computeB(DDEGaussian dde, double[][] data, double[] evalPoints, int maxDerivative, double init,
			DEKernelDensity de) {
		int t = maxDerivative / 2; // maxDerivative = 2t
		double[] mu2Derivatives = new double[t - 2];
		for (int r = t - 1; r >= 2; r--)
			// the first order is 2t-2, this is init.
			// the first in mu2Derivatives is 2t-4, the second 2t-6,..., the last 2*2-2.
			mu2Derivatives[t - 1 - r] = densityFunctionalGaussian(2 * r - 2, 1.0);
		return computeB(dde, data, 1.0, mu2Derivatives, evalPoints, maxDerivative, init, de);

	}
	
	/**
	 * Same as
	 * {@link #computeB(DDEGaussian, double[][], double[], int, double, DEKernelDensity)}
	 * but with \a numEvalPoints equidistant integration nodes.
	 * 
	 * @param dde
	 *            the DDE we use for determining the optimal bandwidth.
	 * @param data
	 *            the observations gained from @f$m@f$ independent simulations.
	 * @param numEvalPoints
	 *            the number of equidistant integration nodes.
	 * @param maxDerivative
	 *            the highest order of the density derivative considered.
	 * @param init
	 *            initial value of the roughness functional in the denominator.
	 * @param de
	 *            the KDE for the sought density.
	 * @return the asymptotically optimal multiplicative constant in the ISB.
	 */
	public static double computeB(DDEGaussian dde, double[][] data, int numEvalPoints, int maxDerivative, double init,
			DEKernelDensity de) {
		return computeB(dde, data, dde.getEquidistantPoints(numEvalPoints), maxDerivative, init, de);
	}

	protected double alpha;
	/**<model parameters for ISB */
	protected double beta;
	/**<model parameters for IV */
	protected double delta;
	/**<model parameters for IV */
	protected double C;
	/**<model parameters for IV */
	protected double B;
	/**<model parameters for ISB */

	protected double a;
	/**<left boundary of the interval over which we estimate. */
	protected double b;

	/**<right boundary of the interval over which we estimate. */

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
	 * value of \a tH.
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
	 * Sets the current value of {@link #alpha} to \a alpha.
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
	 * Sets the current value of {@link #B} to \a B.
	 * 
	 * @param B
	 *            the desired value for {@link #B}
	 */
	public void setB(double B) {
		this.B = B;
	}

	/**
	 * Sets the current value of {@link #beta} to \a beta$.
	 * 
	 * @param beta
	 *            the desired value for {@link #beta}.
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
	 * Sets the current value of {@link #C} to \a C.
	 * 
	 * @param C
	 *            the desired value for {@link #C}.
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
	 * Sets the current value of {@link #delta} to \a delta.
	 * 
	 * @param delta
	 *            the desired value for {@link #delta}.
	 */
	public void setDelta(double delta) {
		this.delta = delta;
	}

	/**
	 * Gives the estimated IV based on the (local) model assumption
	 * \f$Cn^{-\beta}h^{-\delta} \f$. Note that this requires the parameters \a C,
	 * \a beta, and \a delta to be set for this estimator.
	 * 
	 * @param n
	 *            the number of observations.
	 * @return the estimated IV.
	 */
	public double computeEstimatedIV(int n) {
		return C * Math.pow(n, -beta) * Math.pow(h, -delta);
	}

	/**
	 * Gives the estimated ISB for when the exact density of the underlying model is
	 * not known, based on its asymptotic value \f$Bh^{\alpha}\f$. Note that this
	 * requires the parameter \a B to be set.
	 * 
	 * @return the estimated ISB.
	 */
	public double computeEstimatedISB() {
		return B * Math.pow(h, alpha);
	}

	/**
	 * Same as {@link #computeEstimatedISB()} with dummy arguments to overload
	 * {@link #computeISB(MonteCarloModelDensityKnown, double[][], double[])}
	 * for situations, where the true density is not known.
	 * 
	 * @return the estimated ISB.
	 */
	public double computeISB(MonteCarloModelDouble model, double[][] estDensities, double[] evalPoints) {
		return computeEstimatedISB();
	}


	/**
	 * Computes the estimated MISE, i.e. the sum of
	 * {@link #computeEstimatedIV(int)} and {@link #computeDensityISB()}.
	 * 
	 * @param n
	 *            the number of observations.
	 * @return the estimated MISE
	 */
	public double computeEstimatedMISE(int n) {
		return computeEstimatedISB() + computeEstimatedIV(n);
	}

	/**
	 * Computes the semi-empirical MISE for situations, where the true density is
	 * not known. I.e., it takes the sum of the empirical IV
	 * {@link #computeIV(double[][])} and the estimated ISB
	 * {@link #computeEstimatedISB()}. 
	 * 
	 * The estimate of the empirical IV is based on \f$m\f$ realizations of the density estimator,
	 * which have previously been evaluated at the \f$k\f$ points stored in \a
	 * evalPoints. The matrix \a estDensities has
	 * dimensions \f$m\times k \f$, i.e. each row contains the evaluations of one
	 * density estimator.
	 * 
	 * @param the \f$m\times k\f$matrix containing the results of evaluating
	 *            \f$m\f$ densities at \f$k\f$ evaluation points.
	 *            the observations to construct the density.
	 * @return the semi-empirical MISE
	 */
	public double computeMISE(double[][] estDensities) {
		double iv = computeIV(estDensities);
		return iv + computeEstimatedISB();
	}


	/**
	 * Same as {@link #computeMISE(double[][], double[])} but with a dummy
	 * argument \a model to overload
	 * {@link #computeMISE(MonteCarloModelDensityKnown, double[][], double[])}
	 * for cases where the true density is not known.
	 * 
	 * @param model
	 * @param estDensities
	 *           the \f$m\times k\f$matrix containing the results of evaluating
	 *            \f$m\f$ densities at \f$k\f$ evaluation points.
	 *            the observations to construct the density.
	 * @return an estimate for the MISE
	 */
	public double computeMISE(MonteCarloModelDouble model, double[][] estDensities, double[] evalPoints) {
		return computeMISE(estDensities);
	}

	/**
	 * This method estimates the 
	 * IV, ISB, and MISE based on \f$m\f$ realizations of the density estimator,
	 * which have previously been evaluated at the \f$k\f$ points stored in \a
	 * evalPoints, and returns them in an array in this order. The matrix \a estDensities has
	 * dimensions \f$m\times k \f$, i.e. each row contains the evaluations of one
	 * density estimator.
	 * 
	 * The estimate for the IV is obtained from an estimate of the empirical IV {@link #computeIV(double[][])},
	 * the ISB is estimated by {@link #computeEstimatedISB()}, and the MISE-estimate is obtained by their sum.
	 * 
	 * 
	 * @param estDensities
	 *            the \f$m\times k\f$matrix containing the results of evaluating
	 *            \f$m\f$ densities at \f$k\f$ evaluation points.
	 * @return estimates for the IV, the ISB, and the MISE.
	 */
	public double[] computeIVandISBandMISE(double[][] estDensities) {
		double iv = computeIV(estDensities);
		double isb = computeEstimatedISB();
		double[] res = {iv,isb,iv+isb};
		return res;
	}
	
	/**
	 * Same as {@link #computeIVandISBandMISE(double[][])} but with dummy arguments to overload
	 * {@link #computeIVandISBandMISE(MonteCarloModelDensityKnown, double[][], double[])} when the
	 * true density is not known.
	 * 
	 * @param model
	 * @param estDensities
	 * @param evalPoints
	 * @return
	 */
	public double[] computeIVandISBandMISE(MonteCarloModelDouble model, double[][] estDensities,
			double[] evalPoints) {
		return computeIVandISBandMISE(estDensities);
	}

}
