package ift6561examples;

import java.io.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.mcqmctools.*;
/**
 * 
 */

/**
 * @author Pierre L'Ecuyer
 * 
 */
public class San13CMC extends San13Prob {

	double estimate; // Cond. prod. of exceeding x.

	public San13CMC(double x, String fileName) throws IOException {
		super(x, fileName);
	}

	public void simulate(RandomStream stream) {
		for (int k = 0; k < 13; k++) {
			V[k] = dist[k].inverseF(stream.nextDouble());
			if (V[k] < 0.0)
				V[k] = 0.0;
		}
		double prod = 1.0;
		// Incomplete path lengths
		paths[0] = V[1] + V[10];
		paths[1] = V[0] + V[2] + V[10];
		if (paths[0] > paths[1])
			paths[1] = paths[0];
		prod *= dist[5].cdf(x - paths[1]);
		paths[2] = V[0] + V[10];
		prod *= dist[4].cdf(x - paths[2]);
		paths[3] = V[0] + V[3] + V[7] + V[10];
		prod *= dist[9].cdf(x - paths[3]);
		paths[4] = V[0] + V[3] + V[7] + V[12];
		prod *= dist[8].cdf(x - paths[4]);
		paths[5] = V[0] + V[3] + V[11] + V[12];
		prod *= dist[6].cdf(x - paths[5]);
		estimate = 1.0 - prod;
	}

	public double getValue() {
		return estimate;
	}

	public static void main(String[] args) throws IOException {
		int n = 1000000;
		San13CMC san = new San13CMC(90.0, "san13a.dat");
		MonteCarloExperiment.simulateRunsDefaultReportStudent (san, n, new LFSR113(),
				new Tally ("SAN13 example with CMC"), 0.95, 4);
	}
}
