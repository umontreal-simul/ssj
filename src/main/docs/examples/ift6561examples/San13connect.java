package ift6561examples;
import java.io.*;
// import java.util.Scanner;
import umontreal.ssj.rng.*;
// import umontreal.ssj.probdist.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.mcqmctools.*;

/**
 * This class simulates a specific stochastic activity network with 9 nodes and
 * 13 links, taken from Elmaghraby (1977) and used again in L'Ecuyer and Lemieux
 * (2000), "Variance Reduction via Lattice Rules". The goal here is to estimate the
 * probability that nodes 1 and 9 are disconncted.
 * This program is very specific to this example and uses a very naive way to
 * compute the shortest path, by enumerating all six paths!
 */

public class San13connect implements MonteCarloModelDouble {

	double[] r = new double[13];
	boolean[] Y = new boolean[13];
	// We consider the 6 paths that can lead to the sink.
	boolean connected;

	// The constructor.
	public San13connect(double rj) {
		for (int k = 0; k < 13; k++) r[k] = rj;
		// r[4] = 0.99; r[5] = 0.99; r[6] = 0.999;
	}

	public void simulate(RandomStream stream) {
		for (int k = 0; k < 13; k++) {
			Y[k] = (stream.nextDouble() < r[k]);
		}
		connected = false;
		// Path lengths
		if ((Y[1] & Y[5] & Y[10]) |
			(Y[0] & Y[2] & Y[5] & Y[10]) |
			(Y[0] & Y[4] & Y[10]) |
			(Y[0] & Y[3] & Y[7] & Y[9] & Y[10]) |
			(Y[0] & Y[3] & Y[7] & Y[8] & Y[12]) |
			(Y[0] & Y[3] & Y[6] & Y[11] & Y[12])) connected = true;
	}

	public double getPerformance() {
		if (connected)
			return 0.0;
		else
			return 1.0;
	}

	public String toString() {
		String s = "SAN network with 9 nodes and 13 links, from Elmaghraby (1977)\n"
				+ "Estimate prob that nodes 1 and 9 are disconnected";
		return s;
	}

	public static void main(String[] args) throws IOException {
		int n = 1000000 * 1000;
		double rj = 0.999;
		San13connect san = new San13connect(rj);
		Tally statC = new Tally("SAN13 reliability example");
		MonteCarloExperiment.simulateRunsDefaultReportStudent (san, n, new LFSR113(),
				statC, 0.95, 4);
		System.out.println(statC.report(0.95, 8));
	}
}
