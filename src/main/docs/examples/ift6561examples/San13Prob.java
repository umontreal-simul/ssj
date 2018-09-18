package ift6561examples;
import java.io.*;
import java.util.Scanner;
import umontreal.ssj.rng.*;
import umontreal.ssj.mcqmctools.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.stat.*;

/**
 * This class simulates a specific stochastic activity network with 9 nodes and
 * 13 links, taken from Elmaghraby (1977) and used again in L'Ecuyer and Lemieux
 * (2000), "Variance Reduction via Lattice Rules". The goal is to estimate the
 * probability that the length of the longest path exceed a given constant x.
 * This program is very specific to this example and uses a very naive way to
 * compute the shortest path, by enumerating all six paths!
 */

public class San13Prob implements MonteCarloModelDouble {

    int dim = 13;  // dimension.
	double[] V = new double[dim];
	ContinuousDistribution[] dist = new ContinuousDistribution[13];
	// We consider the 6 paths that can lead to the sink.
	double[] paths = new double[6];
	double maxPath; // Length of the current longest path.
	double x;  // Threshold for longest path: we estimate P[T > x].

	
	// The constructor reads link length distributions in a file.
	public San13Prob(double x, String fileName) throws IOException {
		this.x = x;
		readDistributions(fileName);
	}

	public void readDistributions(String fileName) throws IOException {
		// Reads data and construct arrays.
		BufferedReader input = new BufferedReader(new FileReader(fileName));
		Scanner scan = new Scanner(input);
		for (int k = 0; k < 13; k++) {
			dist[k] = DistributionFactory.getContinuousDistribution(scan
					.nextLine());
			// gen[k] = new RandomVariateGen (stream, dist);
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

    // Returns the indicator that T > x.
	public double getIndicatorTx() {
		if (maxPath > x)
			return 1.0;
		else
			return 0.0;
	}

	public void setx(double x) {
		this.x = x;
	}

	public String toString() {
		String s = "SAN network with 9 nodes and 13 links, from Elmaghraby (1977)\n"
				+ "Estimate distribution on length T of longest path \n"
				+ "and P[T > x} for x = " + x + "\n";
		return s;
	}

	public static void main(String[] args) throws IOException {
		int n = 10000000;
		double x = 90.0;
		San13Prob san = new San13Prob(x, "san13a.dat");
		TallyStore statT = new TallyStore("SAN13 example");
		
		// Deterministic approximation:
		System.out.println("Deterministic approximation: T = "
				+ san.deterministicLengths () + "\n\n");
		
		// Run a simulation experiment with n runs.
		MonteCarloExperiment.simulateRunsDefaultReportStudent (san, n, new LFSR113(),
				statT, 0.95, 4);
 	}
}
