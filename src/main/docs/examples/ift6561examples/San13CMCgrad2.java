package ift6561examples;

import java.io.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.stat.*;
// import umontreal.ssj.probdist.*;

/**
 * @author Pierre L'Ecuyer
 * 
 */
public class San13CMCgrad2 extends San13Prob {

    double estimate; // Cond. prod. of exceeding x.
	double epsilon;  // Step size.
	int numx;        // Number of steps (values of x).
	int[] indexV = new int[8];  // The 8 generated variates.
	int[] indexCut = new int[5];  // The 5 arcs of the cut.
	Tally[] statsProb;
	Tally[] statsDeriv;  // Stats on derivative w.r.t. x.

    public San13CMCgrad2 (double x0, double epsilon, int numx, String fileName) throws IOException {
		super(x0, fileName);
		// The CVs: 0, 1, 2, 3, 7, 10, 11, 12
		indexV[0] = 0; indexV[1] = 1;  indexV[2] = 2;  indexV[3] = 3;
		indexV[4] = 7; indexV[5] = 10; indexV[6] = 11; indexV[7] = 12;
		indexCut[0] = 5; indexCut[1] = 4;  indexCut[2] = 9;  indexCut[3] = 8;  indexCut[4] = 6;
		this.epsilon = epsilon;
		this.numx = numx;
		double x;
		statsProb = new Tally[numx];
		statsDeriv = new Tally[numx];  // Stats on derivative.
		for (int j = 0; j < numx; j++)
		{
			x = x0 + j * epsilon;
			statsProb[j] = new Tally("Probability at x = " + Double.toString(x));
			statsDeriv[j] = new Tally("Derivative at x = " + Double.toString(x));
		}
	}

	public void simulateRuns (int n, RandomStream stream) {
		for (int j = 0; j < numx; j++)
		{
			statsProb[j].init();
			statsDeriv[j].init();
		}
		for (int i = 0; i < n; i++)
			simulate(stream);
	}

	public void simulate(RandomStream stream)
	{
		int k;
		double x1;
		double sum;
		// Generate: 0, 1, 2, 3, 7, 10, 11, 12
		for (int j = 0; j < 8; j++)
		{
			k = indexV[j];
			V[k] = dist[k].inverseF(stream.nextDouble());
			if (V[k] < 0.0)
				V[k] = 0.0;
		}
		// Incomplete path lengths
		paths[0] = V[1] + V[10];
		paths[1] = V[0] + V[2] + V[10];
		if (paths[0] > paths[1])
			paths[1] = paths[0];
		paths[2] = V[0] + V[10];
		paths[3] = V[0] + V[3] + V[7] + V[10];
		paths[4] = V[0] + V[3] + V[7] + V[12];
		paths[5] = V[0] + V[3] + V[11] + V[12];
		for (int j = 0; j < numx; j++) {
			x1 = x + (double)j * epsilon;
			double prod = 1.0;
			for (int i = 0; i < 5; i++) 
			   prod *= dist[indexCut[i]].cdf(x1 - paths[i+1]);
		    statsProb[j].add(1.0 - prod);
			sum = 0.0;
   		    if (prod > 0.0) 
				for (int i = 0; i < 5; i++)
				   sum += dist[indexCut[i]].density(x1 - paths[i+1]) 
				   	      * prod / dist[indexCut[i]].cdf(x1 - paths[i+1]);
			statsDeriv[j].add(sum);
		}
	}

	public String toString() {
		String s = "SAN network with 9 nodes and 13 links, from Elmaghraby (1977)\n"
				+ "Estimate prob longest path > x and gradient w.r.t. x.";
		return s;
	}

	public static void main(String[] args) throws IOException {
		int n = 100000;
		San13CMCgrad2 san = new San13CMCgrad2 (90.0, 0.2, 11, "san13a.dat");
		san.simulateRuns(n, new LFSR113());
		for (int j = 0; j < san.numx; j++)
			System.out.println(san.statsProb[j].reportAndCIStudent(0.95, 6));
		System.out.println("--------------------------------------------\n");
		for (int j = 0; j < san.numx; j++)
			System.out.println(san.statsDeriv[j].reportAndCIStudent(0.95, 6));
	}
}
