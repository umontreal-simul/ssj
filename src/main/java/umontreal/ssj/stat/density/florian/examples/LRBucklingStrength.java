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
 * This class implements a likelihood-ratio density estimator for the buckling
 * strength @f$X@f$ of a steel plate as given by Schields and Zhang '16. The
 * model is given by
 * 
 * \f[ X=\left( \frac{2.1}{\lambda}-\frac{0.9}{\lambda^2} \right) \left(
 * 1-\frac{0.75\delta_0}{\lambda}\right) \left(1-\frac{2\eta t}{b} \right), \f]
 * 
 * 
 * where \lambda = (b/t) \sqrt{\sigma_0/eta} and where the
 * parameters @f$b,t,\sigma_0,E,\delta,\eta@f$ are random variables. More
 * precisely, $t$ and $b_0$ follow a lognormal distribution and all the others
 * are normally distributed.
 * 
 * The realizations of the random variables at which this estimator can be constructed can 
 * be obtained with the class #BucklingStrengthVars.
 * 
 * This estimator allows for free parameters which can be passed as weights.
 * @remark **Florian:** We need to re-check what the free parameters are. Do they need to sum up to 1?
 */

public class LRBucklingStrength extends ConditionalDensityEstimator {

	double mud, mue; 
	double sigmad, sigmae;
	double [] weights;
	
	/**
	 * Constructor
	 * @param mud mean of @f$\delta_0@f$.
	 * @param sigmad standard deviation of @f$\delta_0@f$
	 * @param mue mean of @f$\eta@f$.
	 * @param sigmae standard deviation of @f$\eta@f$.
	 * @param weights free parameters.
	 */
	public LRBucklingStrength(double mud, double sigmad, double mue, double sigmae,double[] weights) {
		this.mue = mue;
		this.sigmae = sigmae;
		this.mud = mud;
		this.sigmad = sigmad;
		this.weights = weights;
	}
	
	private double strength(double[] data) {
		return ((2.1 / data[1] - 0.9 / (data[1] * data[1])) * (1.0 - 0.75 * data[2] / data[1])
				* (1.0 - data[0] * data[3]));
	}

	@Override
	public double evalEstimator(double x, double[] data) {
		if (strength(data) > x)
			return 0;
		else {
		
			return ( weights[0] * (-4.0  * data[1]/3.0 + data[2]) * (mud  - data[2])/(sigmad * sigmad) + weights[1] * (data[3] - 1.0/data[0]) * (mue - data[3])/(sigmae * sigmae)+1 ) / x;
		}
	}
	
	public String toString() {
		return "LRBucklingStrength";
	}

	public double[] getWeights() {
		return weights;
	}

	public void setWeights(double [] weights) {
		this.weights = weights;
	}

	
}
