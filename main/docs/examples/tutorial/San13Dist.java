package tutorial;

import java.io.*;
import java.util.Scanner;
import umontreal.ssj.charts.HistogramChart;
import umontreal.ssj.charts.EmpiricalChart;
import umontreal.ssj.probdist.ContinuousDistribution;
import umontreal.ssj.probdist.DistributionFactory;
import umontreal.ssj.rng.*;
import umontreal.ssj.mcqmctools.*;
import umontreal.ssj.stat.TallyStore;

/**
 * This class simulates a specific stochastic activity network with 9 nodes and
 * 13 links, taken from Elmaghraby (1977) and used again in L'Ecuyer and Lemieux
 * (2000), "Variance Reduction via Lattice Rules". The goal is to estimate the
 * distribution of the length T of the shortest path and make a histogram.
 * This program is very specific to this example and uses a very naive way to
 * compute the shortest path, by enumerating all six paths!
 */

public class San13Dist implements MonteCarloModelDouble {

    int dim = 13;  // model dimension (number of uniforms needed for a simulation). 
	double[] V = new double[dim];
	ContinuousDistribution[] dist = new ContinuousDistribution[dim];
	// We consider the 6 paths that can lead to the sink.
	double[] paths = new double[6];  // Path lengths.
	double maxPath; // Length of the current longest path.

	
	// The constructor reads link length distributions in a file.
	public San13Dist (String fileName) throws IOException {
		readDistributions(fileName);
	}

	public void readDistributions(String fileName) throws IOException {
		// Reads data and construct arrays.
		BufferedReader input = new BufferedReader(new FileReader(fileName));
		Scanner scan = new Scanner(input);
		for (int k = 0; k < 13; k++) {
			dist[k] = DistributionFactory.getContinuousDistribution(scan
					.nextLine());
		}
		scan.close();
	}

	public int getDimension () {
	   return dim;
    }

	public double deterministicLengths () {
		for (int k = 0; k < 13; k++) {
			V[k] = dist[k].getMean();
			if (V[k] < 0.0)
				V[k] = 0.0;
		}
		return computePathsAndT();
	}

	public void simulate (RandomStream stream) {
		for (int k = 0; k < 13; k++) {
			V[k] = dist[k].inverseF(stream.nextDouble());
			if (V[k] < 0.0)
				V[k] = 0.0;
		}
		computePathsAndT();
	}
	
    // Compute the lengths of all paths and returns the longest length T
	public double computePathsAndT () {
		// Path lengths
		paths[0] = V[1] + V[5] + V[10];
		paths[1] = V[0] + V[2] + V[5] + V[10];
		paths[2] = V[0] + V[4] + V[10];
		paths[3] = V[0] + V[3] + V[7] + V[9] + V[10];
		paths[4] = V[0] + V[3] + V[7] + V[8] + V[12];
		paths[5] = V[0] + V[3] + V[6] + V[11] + V[12];
		maxPath = paths[0];
		for (int p = 1; p < 6; p++)
			if (paths[p] > maxPath)
				maxPath = paths[p];
		return maxPath;
	}
	
	// Returns the length T of longest path.
	public double getPerformance() {
		return maxPath;
	}

	public String toString() {
		String s = "SAN network with 9 nodes and 13 links, from Elmaghraby (1977)\n"
				   + "Estimate distribution of length T of longest path \n";
		return s;
	}

	public static void main(String[] args) throws IOException {
		int n = 100000;
        int groupSize = 100;   
		San13Dist san = new San13Dist("san13a.dat");
		TallyStore statT = new TallyStore("TallyStore for SAN13 example");
		System.out.println (MonteCarloExperiment.simulateRunsDefaultReportStudent 
			(san, n, new LFSR113(), statT, 0.95, 4));
		statT.quickSort();

		Writer file;
		TallyStore statTaggregated = statT;
        if (groupSize > 1) 
        	statTaggregated = statT.aggregate (groupSize);
        double [] aggreg = statTaggregated.getArray();
		EmpiricalChart cdf = new EmpiricalChart("Empirical cdf of $T$", 
				 "Values of $T$", "cdf", aggreg, statTaggregated.numberObs());
		double[] bounds2 = { 0, 200, 0, 1.0 };
		cdf.setManualRange(bounds2);
		cdf.view(800, 500);
		String cdfLatex = cdf.toLatex(12.0, 8.0);
		file = new FileWriter("san13cdf.tex");
		file.write(cdfLatex);
		file.close();
		
		HistogramChart hist = new HistogramChart("Distribution of $T$",
				"Values of $T$", "Frequency", statT.getArray(), n);
		double[] bounds = { 0, 200, 0, 12000 };
		hist.setManualRange(bounds);
		(hist.getSeriesCollection()).setBins(0, 40, 0, 200);
		hist.view(800, 500);
		String histLatex = hist.toLatex(12.0, 8.0);
		file = new FileWriter("san13chart.tex");
		file.write(histLatex);
		file.close();


		// Print p-th quantile
		double p = 0.99;
		int index = (int)Math.round (p * n);
		double xip = statT.getArray()[index];
		System.out.printf("%5.3g -th quantile: %9.6g \n", p, xip);
	}
}
