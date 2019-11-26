package umontreal.ssj.stat.density.florian.examples;

import java.util.Arrays;

import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.stat.density.ConditionalDensityEstimator;

/**
 * This class implements a conditional density estimator for the displacement of
 * a cantilever beam
 * 
 * @f$f(E,X,Y)@f$ as explained in \link https://www.sfu.ca/~ssurjano/canti.html
 *                \endlink. More precisely,
 * 
 *                @f[ f(E,X,Y) = \frac{4L^3}{Ewt}\sqrt{ \frac{Y^2}{t^4} +
 *                \frac{X^2}{w^4} }, @f] 
 * 
 *                where @f$L,w@f$ and @f$t@f$ are constants. The random
 *                variables @f$E,X,Y@f$ are assumed to be normally distributed
 *                and independent.
 * 
 *                In fact, this class derives three estimators (each variable is
 *                hidden once), and consider their weighted sum as the CDE.
 *                Weights are supposed to be positive and sum up to one. If a negative
 *                weight is passed, it is automatically set to zero.
 * @author florian
 *
 */
public class CDECantilever extends ConditionalDensityEstimator {

	/**
	 * Weights for the variables E, X, and Y (in this order).
	 */
	private double[] weights;

	private double alpha, beta, gamma; // alpha = 4L^3/(wt); beta = w^4; gamma = t^4;

	private double muE, sigmaE, muX, sigmaX, muY, sigmaY;

	/**
	 * Constructor.
	 * 
	 * @param L length of the beam
	 * @param t thickness of the beam
	 * @param w width of the beam
	 * @param muE mean of @f$E@f$.
	 * @param sigmaE standard deviation of @f$E@f$.
	 * @param muX mean of @f$X@f$.
	 * @param sigmaX standard deviation of @f$X@f$.
	 * @param muY mean of @f$Y@f$.
	 * @param sigmaY standard deviation of @f$Y@f$.
	 * @param weights the weights of the convex combination.
	 */


	public CDECantilever(double L, double t, double w, double muE, double sigmaE, double muX, double sigmaX, double muY,
			double sigmaY, double[] weights) {
		this.weights = new double[3];
		Arrays.fill(weights, 1.0 / 3.0);
		this.alpha = 4.0 * L * L * L / (w * t);
		this.beta = w * w * w * w;
		this.gamma = t * t * t * t;
		this.muE = muE;
		this.sigmaE = sigmaE;
		this.muX = muX;
		this.sigmaX = sigmaX;
		this.muY = muY;
		this.sigmaY = sigmaY;
		setWeights(weights);
	}

	/**
	 * Setter for the weights. Normalizes the weights to 1 in \f$\ell^1\f$.
	 * 
	 * @param weights
	 */
	public void setWeights(double[] weights) {
		this.weights = new double[3];
		for (int i = 0; i < 3; i++)
			this.weights[i] = weights[i];
//		normalizeWeights();
	}

	

	private double deltaSqX(double x, double[] data) {
		return (beta * (x * x * data[0] * data[0] / (alpha * alpha) - data[2] * data[2] / gamma));
	}

	private double deltaSqY(double x, double[] data) {
		return (gamma * (x * x * data[0] * data[0] / (alpha * alpha) - data[1] * data[1] / beta));
	}

	@Override
	public double evalEstimator(double x, double[] data) {
		double val = 0.0;
		double arg;
		double temp;

		if (weights[0] >= 0) {
			arg = alpha * Math.sqrt(data[1] * data[1] / beta + data[2] * data[2] / gamma) / x;
			val += weights[0] * NormalDist.density(muE, sigmaE, arg) * arg / x;
		}

		temp = deltaSqX(x, data);

		if ((weights[1] >= 0) && (temp >= 0)) {
			arg = Math.sqrt(temp);
			temp = weights[1] * x * data[0] * data[0] * beta / (alpha * alpha * arg);
			val += temp * (NormalDist.density(muX, sigmaX, arg) + NormalDist.density(muX, sigmaX, -arg));
		}


		temp = deltaSqY(x, data);
		if (weights[2] >= 0) {
			arg = Math.sqrt(temp);
			temp = weights[2] * x * data[0] * data[0] * gamma / (alpha * alpha * arg);
			val += temp * (NormalDist.density(muY, sigmaY, arg) + NormalDist.density(muY, sigmaY, -arg));
		}

		return val;
	}

	@Override
	public String toString() {
		return "CDECantilever";
	}

}
