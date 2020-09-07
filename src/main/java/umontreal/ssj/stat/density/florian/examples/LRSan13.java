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
 * This class implements a likelihood-ratio density estimator for a specific
 * stochastic activity network with 9 nodes and 13 links, taken from Elmaghraby
 * (1977) and used again in L'Ecuyer and Lemieux (2000), "Variance Reduction via
 * Lattice Rules". The goal is to estimate the distribution of the length of the
 * longest path.
 * 
 * The necessary values to construct this estimator can be obtained with #
 * San13Vars.
 * 
 * @author florian
 *
 */

public class LRSan13 extends ConditionalDensityEstimator {
	double paths[];

	public void LLSan() {
		paths = new double[6];
	}

	private double computePaths(double[] data) {
		double[] paths = new double[6];
		double maxPath;
		paths[0] = data[0] + data[2] + data[5] + data[10];
		paths[1] = data[0] + data[3] + data[6] + data[11] + data[12];
		paths[2] = data[0] + data[3] + data[7] + data[8] + data[12];
		paths[3] = data[0] + data[3] + data[7] + data[9] + data[10];
		paths[4] = data[0] + data[4] + data[10];
		paths[5] = data[1] + data[5] + data[10];
		maxPath = paths[0];
		for (int p = 1; p < 6; p++)
			if (paths[p] > maxPath)
				maxPath = paths[p];
		return maxPath;
	}

	@Override
	public double evalEstimator(double x, double[] data) {
		if (computePaths(data) > x)
			return 0.0;
		return ((13.0 - data[0]) / (3.25 * 3.25) * data[0] + (5.5 - data[1]) / (1.375 * 1.375) * data[1] - data[2] / 7.0
				+ (5.2 - data[3]) / (1.3 * 1.3) * data[3] - data[4] / 16.5 - data[5] / 14.7 - data[6] / 10.3
				- data[7] / 6.0 - data[8] * 0.25 - data[9] * 0.05 + (3.2 - data[10]) / (0.8 * 0.8) * data[10]
				+ (3.2 - data[11]) / (0.8 * 0.8) * data[11] - data[12] / 16.5 + 13.0) / x;
	}

	@Override
	public String toString() {
		return "LRSan13";
	}

	

}
