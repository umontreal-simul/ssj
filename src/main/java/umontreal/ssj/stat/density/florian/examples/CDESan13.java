package umontreal.ssj.stat.density.florian.examples;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

import umontreal.ssj.hups.LMScrambleShift;
import umontreal.ssj.hups.PointSet;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.RQMCPointSet;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.mcqmctools.MonteCarloModelDoubleArray;
import umontreal.ssj.mcqmctools.RQMCExperiment;
import umontreal.ssj.probdist.ContinuousDistribution;
import umontreal.ssj.probdist.DistributionFactory;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.density.ConditionalDensityEstimator;
import umontreal.ssj.stat.list.ListOfTallies;

/**
 * This class implements a CDE for a specific stochastic activity network with 9
 * nodes and 13 links, taken from Elmaghraby (1977) and used again in L'Ecuyer
 * and Lemieux (2000), "Variance Reduction via Lattice Rules". The goal is to
 * estimate the distribution of the length of the longest path.
 * 
 * The necessary values to construct this CDE can be obtained with #
 * San13VarsCDE.
 * 
 * @author florian
 *
 */
public class CDESan13 extends ConditionalDensityEstimator {

	public ContinuousDistribution[] dist = new ContinuousDistribution[5];

	/**
	 * A constructor which reads the link length distributions in a file. Be sure to
	 * provide all 13 distributions as for the original San13 example, from which
	 * the sample data is generated.
	 * 
	 * @param fileName file containing the distributions for each link.
	 * @throws IOException
	 */
	public CDESan13(String fileName) throws IOException {
		readDistributions(fileName);

	}

	private void readDistributions(String fileName) throws IOException {
		// Reads data and construct arrays.
		BufferedReader input = new BufferedReader(new FileReader(fileName));
		Scanner scan = new Scanner(input);
		int i = 0;
		for (int k = 0; k < 13; k++) {

			if (k == 4 || k == 5 || k == 6 || k == 8 || k == 9) {
				dist[i] = DistributionFactory.getContinuousDistribution(scan.nextLine());
//			System.out.println("distr. " + i + " is " + dist[i].toString());
				i++;
			} else
				scan.nextLine();
		}
		scan.close();
	}

	@Override
	public double evalEstimator(double x, double[] data) {
		int t = data.length;
		double dens = 0.0;
		for (int j = 0; j < t; j++) {
			double prod = 1.0;
			for (int k = 0; k < t; k++) {
				if (k != j)
//					System.out.println("Zeigs:" + data[k]);
					prod *= dist[k].cdf(x - data[k]);
			}
			dens += dist[j].density(x - data[j]) * prod;
		}

		return dens;
	}

	public String toString() {
		return "CDESan13";
	}

}
