package umontreal.ssj.mcqmctools.florian.examples;

import cern.colt.matrix.DoubleMatrix2D;
import umontreal.ssj.mcqmctools.MonteCarloModelDoubleArray;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.randvarmulti.MultinormalPCAGen;
import umontreal.ssj.rng.RandomStream;

/**
 * This class can be used to generate the random variables @f$Y,M,P@f$ that are
 * used in the so-called short column function @f$f(Y,M,P)@f$ as explained in
 * \link https://www.sfu.ca/~ssurjano/shortcol.html \endlink. More precisely,
 * 
 * @f[ f(Y,M,P) = 1 - \frac{4M}{bh^2Y}-\frac{P^2}{b^2h^2Y^2}, @f]
 * 
 * where @f$b@f$ and @f$h@f$ are the deterministic width and depth of the cross
 * section (in mm), respectively.
 * 
 * This class can be used to, e.g., estimate the density of @f$f(Y,M,P)@f$ with
 * a conditional density estimator.
 * 
 * The random variable @f$Y@f$ is assumed to be lognormally distributed
 * while @f$M@f$ and @f$P@f$ are normally distributed with covariance
 * Matrix @f$\Sigma@f$. @f$Y@f$ is independent of @f$M@f$ and @f$P@f$.
 * 
 * @author florian
 *
 */

public class ShortColumnFunctionVars implements MonteCarloModelDoubleArray {

	double muY, muM, muP;
	double sigmaY;
	DoubleMatrix2D trafoMat;
	double[] performance = new double[3];
	double[][] sigma;

	/**
	 * Constructor
	 * @param muY    mean of @f$Y@f$.
	 * @param muM    mean of @f$M@f$.
	 * @param muP    mean of @f$P@f$.
	 * @param sigmaY standard deviation of @f$Y@f$.
	 * @param sigma  covariance matrix @f$\Sigma@f$ of @f$M@f$ and @f$P@f$.
	 */
	public ShortColumnFunctionVars(double h, double b, double muY, double muM, double muP, double sigmaY,
			DoubleMatrix2D sigma) {

		this.muY = transformMu(muY, sigmaY);
		this.muM = muM;
		this.muP = muP;
		this.sigmaY = transformSigma(muY, sigmaY);
		this.trafoMat = MultinormalPCAGen.decompPCA(sigma);
	}

	/**
	 * Same as above, but now @f$\Sigma@f$ is given as a two-dimensionla array
	 * instead of a DoubleMatrix2D.
	 * 
	 * @param muY
	 * @param muM
	 * @param muP
	 * @param sigmaY
	 * @param sigma
	 */
	public ShortColumnFunctionVars(double muY, double muM, double muP, double sigmaY, double[][] sigma) {
		this.muY = transformMu(muY, sigmaY);
		this.muM = muM;
		this.muP = muP;
		this.sigmaY = transformSigma(muY, sigmaY);
		this.trafoMat = MultinormalPCAGen.decompPCA(sigma);
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
		performance[0] = Math.exp(NormalDist.inverseF(muY, sigmaY, stream.nextDouble()));
		double[] u = new double[2];
		u[0] = NormalDist.inverseF01(stream.nextDouble());
		u[1] = NormalDist.inverseF01(stream.nextDouble());
		for (int j = 0; j < 2; j++) {
			performance[j + 1] = 0.0;
			for (int c = 0; c < 2; c++)
//				performance[j+1] += trafoMat.getQuick(j, c) * u[c];
				performance[j + 1] += sigma[j][c] * u[c];

		}
		performance[1] += muM;
		performance[2] += muP;

	}

	@Override
	public double[] getPerformance() {
		return performance;
	}

	@Override
	public int getPerformanceDim() {
		return performance.length;
	}

	@Override
	public String toString() {
		return "ShortColumn";
	}

}
