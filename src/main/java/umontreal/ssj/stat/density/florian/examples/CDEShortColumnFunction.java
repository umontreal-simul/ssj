package umontreal.ssj.stat.density.florian.examples;

import java.io.FileWriter;
import java.io.IOException;

import umontreal.ssj.hups.LMScrambleShift;
import umontreal.ssj.hups.PointSet;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.RQMCPointSet;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.mcqmctools.MonteCarloModelDoubleArray;
import umontreal.ssj.mcqmctools.RQMCExperiment;

import umontreal.ssj.probdist.LognormalDist;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.density.ConditionalDensityEstimator;
import umontreal.ssj.stat.list.ListOfTallies;

/**
 * This class implements a conditional density estimator for the Short column
 * function
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
 *                The hidden variable is @f$Y@f$ and the realizations
 *                of @f$(Y,M,P)@f$ can be obtained by the class
 *                #ShortColumnFunctionVars.
 * @author florian
 *
 */
public class CDEShortColumnFunction extends ConditionalDensityEstimator {

	double b, h;
	double muY, sigmaY;

	/**
	 * Constructor
	 * @param b width of the cross section in mm.
	 * @param h height of the cross section in mm.
	 * @param muY mean of @f$Y@f$.
	 * @param sigmaY standard deviation of @f$Y@f$.
	 */
	public CDEShortColumnFunction(double b, double h, double muY, double sigmaY) {
		this.b = b;
		this.h = h;
		this.muY = transformMu(muY, sigmaY);
		this.sigmaY = transformSigma(muY, sigmaY);
	}

	private double transformMu(double mu, double sigma) {
		return (Math.log(mu) - 0.5 * Math.log(sigma * sigma / (mu * mu) + 1.0));
	}

	private double transformSigma(double mu, double sigma) {
		return (Math.sqrt(Math.log(1.0 + sigma * sigma / (mu * mu))));
	}

	@Override
	public double evalEstimator(double x, double[] data) {
//		double 1mX = 1.0 - x;
		// Sqrt[ 4 M^2 +( P^2 h^2 (1 -z)) ]
		double temp1 = Math.sqrt(4.0 * data[1] * data[1] + (data[2] * data[2] * h * h * (1.0 - x)));
		double arg = (2.0 * data[1] + temp1) / (b * h * h * (1.0 - x));

		double val = LognormalDist.density(muY, sigmaY, arg);
		val *= (arg / (1.0 - x) - 0.5 * data[2] * data[2] / (b * (1.0 - x) * temp1));
		return val;
	}

	
}
