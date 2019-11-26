package umontreal.ssj.stat.density.florian.examples;

import java.io.FileWriter;
import java.io.IOException;

import umontreal.ssj.hups.LMScrambleShift;
import umontreal.ssj.hups.PointSet;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.RQMCPointSet;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.mcqmctools.MonteCarloModelDoubleArray;
import umontreal.ssj.mcqmctools.RQMCExperiment;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.density.ConditionalDensityEstimator;
import umontreal.ssj.stat.list.ListOfTallies;

/**
 * This class implements a likelihood-ratio density estimator (LR) for the Short
 * column function
 * 
 * @f$f(Y,M,P)@f$ as explained in \link
 *                https://www.sfu.ca/~ssurjano/shortcol.html \endlink. More
 *                precisely,
 * 
 *                @f[ f(Y,M,P) = 1 - \frac{4M}{bh^2Y}-\frac{P^2}{b^2h^2Y^2}, @f]
 * 
 *                where @f$b@f$ and @f$h@f$ are the deterministic width and
 *                depth of the cross section (in mm), respectively. The random
 *                variable @f$Y@f$ is assumed to be lognormally distributed
 *                while @f$M@f$ and @f$P@f$ are normally distributed with
 *                covariance Matrix @f$\Sigma@f$. @f$Y@f$ is independent
 *                of @f$M@f$ and @f$P@f$.
 * 
 *                The final LR allows for a free parameter @f$p@f$ and
 *                data to construct this estimators can be obtained by the class
 *                #ShortColumnFunctionVars.
 * @author florian
 *
 */

public class LRShortColumnFunction extends ConditionalDensityEstimator {

	double b, h;
	double muY, sigmaY;
	double p;

	/**
	 * Constructor that sets @f$p=1/2@f$.
	 * @param b width of the cross section in mm.
	 * @param h height of the cross section in mm.
	 * @param muY mean of @f$Y@f$.
	 * @param sigmaY standard deviation of @f$Y@f$.
	 */
	public LRShortColumnFunction(double b, double h, double muY, double sigmaY) {
		this.b = b;
		this.h = h;
		this.muY = transformMu(muY, sigmaY);
		this.sigmaY = transformSigma(muY, sigmaY);
		this.p = 0.5;
	}

	/**
	 * Constructor
	 * @param b width of the cross section in mm.
	 * @param h height of the cross section in mm.
	 * @param muY mean of @f$Y@f$.
	 * @param sigmaY standard deviation of @f$Y@f$.
	 * @param p free parameter.
	 */
	public LRShortColumnFunction(double b, double h, double muY, double sigmaY, double p) {
		this.b = b;
		this.h = h;
		this.muY = transformMu(muY, sigmaY);
		this.sigmaY = transformSigma(muY, sigmaY);
		this.p = p;
	}

	private double transformMu(double mu, double sigma) {
		return (Math.log(mu) - 0.5 * Math.log(sigma * sigma / (mu * mu) + 1.0));
	}

	private double transformSigma(double mu, double sigma) {
		return (Math.sqrt(Math.log(1.0 + sigma * sigma / (mu * mu))));
	}

	private double g(double[] data) {
		return (1.0 - 4.0 * data[1] / (b * h * h * data[0]) - data[2] * data[2] / (b * b * h * h * data[0] * data[0]));
	}

	@Override
	public double evalEstimator(double x, double[] data) {
		if (g(data) >= x)
			return 0;
		if (data[2] < 0)
			System.out.println("p<0");
		double fac1 = data[1] - 0.25 * b * h * h * data[0];
		double fac2 = (2.0 * (500.0 + data[2]) - data[1]) / 120000.0;
		double fac3 = Math.abs(data[2]) * (2000.0 + data[1] - 8.0 * data[2]) / 60000.0;
		return ((fac1 * fac2) + 0.5 * (fac3) + 0.5 + 1.0) / x;

	}

	

}
