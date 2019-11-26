package umontreal.ssj.mcqmctools.florian.examples;

import cern.colt.matrix.DoubleMatrix2D;
import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.randvarmulti.MultinormalPCAGen;
import umontreal.ssj.rng.RandomStream;

/**
 * This class implements the so-called short column function @f$f(Y,M,P)@f$ as
 * explained in \link https://www.sfu.ca/~ssurjano/shortcol.html \endlink.
 * More precisely,
 * 
 * @f[ f(Y,M,P) = 1 - \frac{4M}{bh^2Y}-\frac{P^2}{b^2h^2Y^2}, @f]
 * 
 * where @f$b@f$ and @f$h@f$ are the deterministic width and depth of the cross
 * section (in mm), respectively. The random variable @f$Y@f$ is assumed to be
 * lognormally distributed while @f$M@f$ and @f$P@f$ are normally distributed
 * with covariance Matrix @f$\Sigma@f$. @f$Y@f$ is independent of @f$M@f$
 * and @f$P@f$.
 * 
 * @author florian
 *
 */

public class ShortColumnFunction implements MonteCarloModelDouble {

	double h;
	double b;
	double muY, muM, muP;
	double sigmaY;
	DoubleMatrix2D trafoMat;
	double performance;
	double[][] sigma;

	/**
	 * Constructor
	 * 
	 * @param h      height of cross section in mm.
	 * @param b      depth of cross section in mm.
	 * @param muY    mean of @f$Y@f$.
	 * @param muM    mean of @f$M@f$.
	 * @param muP    mean of @f$P@f$.
	 * @param sigmaY standard deviation of @f$Y@f$.
	 * @param sigma  covariance matrix @f$\Sigma@f$ of @f$M@f$ and @f$P@f$.
	 */
	public ShortColumnFunction(double h, double b, double muY, double muM, double muP, double sigmaY,
			DoubleMatrix2D sigma) {
		this.h = h;
		this.b = b;
		this.muY = muY;
		this.muM = muM;
		this.muP = muP;
		this.sigmaY = sigmaY;
		this.trafoMat = MultinormalPCAGen.decompPCA(sigma);
//		this.trafoMat = new CholeskyDecomposition(sigma));
//		this.trafoMat = DoubleMatrix2D.identity(2);

	}

	/**
	 * Same as above, but now @f$\Sigma@f$ is given as a two-dimensionla array
	 * instead of a DoubleMatrix2D.
	 * 
	 * @param h
	 * @param b
	 * @param muY
	 * @param muM
	 * @param muP
	 * @param sigmaY
	 * @param sigma
	 */
	public ShortColumnFunction(double h, double b, double muY, double muM, double muP, double sigmaY,
			double[][] sigma) {
		this.h = h;
		this.b = b;
		this.muY = transformMu(muY, sigmaY);
		this.muM = muM;
		this.muP = muP;
		this.sigmaY = transformSigma(muY, sigmaY);
		this.trafoMat = MultinormalPCAGen.decompPCA(sigma);
//		this.trafoMat = DoubleFactory2D.dense.make(sigma);
//		System.out.println(trafoMat);
		this.sigma = sigma;
	}

	private double transformMu(double mu, double sigma) {
		return (Math.log(mu) - 0.5 * Math.log(sigma * sigma / (mu * mu) + 1.0));
	}

	private double transformSigma(double mu, double sigma) {
		return (Math.sqrt(Math.log(1.0 + sigma * sigma / (mu * mu))));
	}

	@Override
	public void simulate(RandomStream stream) {
		double Y = Math.exp(NormalDist.inverseF(muY, sigmaY, stream.nextDouble()));
		double[] u = new double[2];
		double[] z = new double[2];
		u[0] = NormalDist.inverseF01(stream.nextDouble());
		u[1] = NormalDist.inverseF01(stream.nextDouble());
		for (int j = 0; j < 2; j++) {
			z[j] = 0;
			for (int c = 0; c < 2; c++)
//				z[j] += trafoMat.getQuick(j, c) * u[c];
				z[j] += sigma[j][c] * u[c];
		}
		z[0] += muM;
		z[1] += muP;
		double fac = 1.0 / (b * h * h * Y);
//		System.out.println("(" + z[0] + ", " + z[1] + ")");
		performance = 1.0 - 4.0 * z[0] * fac - z[1] * z[1] * fac / (b * Y);
	}

	@Override
	public double getPerformance() {
		return performance;
	}

}
