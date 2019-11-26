package umontreal.ssj.stat.density.florian.examples;

import umontreal.ssj.stat.density.ConditionalDensityEstimator;

/**
 * This class implements a likelihood-ratio density estimator for the displacement of
 * a cantilever beam
 * 
 * @f$f(E,X,Y)@f$ as explained in \link https://www.sfu.ca/~ssurjano/canti.html
 *                \endlink. More precisely,
 * 
 *                @f[ f(E,X,Y) = \frac{4L^3}{Ewt}\sqrt{ \frac{Y^2}{t^4} +
 *                \frac{X^2}{w^4} }, @f] / D_0
 * 
 *                where @f$L,w@f$ and @f$t@f$ are constants. The random
 *                variables @f$E,X,Y@f$ are assumed to be normally distributed
 *                and independent.
 * 
 * 
 *                This estimator allows for one free parameter @f$p@f$p.
 * @author florian
 *
 */

public class LRCantilever extends ConditionalDensityEstimator {
	private double p; // Scaling factor

	private double alpha, beta, gamma; // alpha = 4L^3/(wt); beta = w^4; gamma = t^4;

	private double muE, sigmaE, muX, sigmaX, muY, sigmaY;

	/**
	 * Constructor with @f$p=1/3@f$.
	 * @param L length of the beam
	 * @param t thickness of the beam
	 * @param w width of the beam
	 * @param muE mean of @f$E@f$.
	 * @param sigmaE standard deviation of @f$E@f$.
	 * @param muX mean of @f$X@f$.
	 * @param sigmaX standard deviation of @f$X@f$.
	 * @param muY mean of @f$Y@f$.
	 * @param sigmaY standard deviation of @f$Y@f$.
	 */
	public LRCantilever(double L, double t, double w, double muE, double sigmaE, double muX, double sigmaX, double muY,
			double sigmaY) {
		p = 1.0 / 3.0;
		this.alpha = 4.0 * L * L * L / (w * t);
		this.beta = w * w * w * w;
		this.gamma = t * t * t * t;
		this.muE = muE;
		this.sigmaE = sigmaE;
		this.muX = muX;
		this.sigmaX = sigmaX;
		this.muY = muY;
		this.sigmaY = sigmaY;
	}

	/**
	 * Constructor
	 * @param L length of the beam
	 * @param t thickness of the beam
	 * @param w width of the beam
	 * @param muE mean of @f$E@f$.
	 * @param sigmaE standard deviation of @f$E@f$.
	 * @param muX mean of @f$X@f$.
	 * @param sigmaX standard deviation of @f$X@f$.
	 * @param muY mean of @f$Y@f$.
	 * @param sigmaY standard deviation of @f$Y@f$.
	 * @param p free parameter.
	 */
	public LRCantilever(double L, double t, double w, double muE, double sigmaE, double muX, double sigmaX, double muY,
			double sigmaY, double p) {
		this(L, t, w, muE, sigmaE, muX, sigmaX, muY, sigmaY);
		this.p = p;
	}

	@Override
	public double evalEstimator(double x, double[] data) {
		double E = data[0];
		double X = data[1];
		double Y = data[2];
		if ((alpha / E * Math.sqrt(X * X / beta + Y * Y / gamma)) > x)
			return 0;
		else
			return (  (- E * (p-1) * (E - muE)/(sigmaE*sigmaE)- X * p * (X - muX)/(sigmaX*sigmaX) - Y * p * (Y - muY)/(sigmaY*sigmaY) + 3.0 *p -1.0) / x  );
	}

	public double getP() {
		return p;
	}

	public void setP(double p) {
		this.p = p;
	}

}
