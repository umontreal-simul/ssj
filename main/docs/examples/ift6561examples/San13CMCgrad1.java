package ift6561examples;

import java.io.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.stat.*;
// import umontreal.ssj.probdist.*;

/**
 * @author Pierre L'Ecuyer
 * 
 */
public class San13CMCgrad1 extends San13Prob {

    double estimate; // Cond. prod. of exceeding x.
	int[] indexV = new int[8];  // The 8 control variates.
	double[] meansCV = new double[8];  // CV expectations.
	Tally statsDeriv3 = new Tally("Deriv. theta3");  // Stats on derivative.
	Tally statsDeriv5 = new Tally("Deriv. theta5");  // Stats on derivative.
	Tally statsProb = new TallyStore("Original CMC estimator");

    public San13CMCgrad1 (double x, String fileName) throws IOException {
		super(x, fileName);
		// The CVs: 0, 1, 2, 3, 7, 10, 11, 12
		indexV[0] = 0; indexV[1] = 1;  indexV[2] = 2;  indexV[3] = 3;
		indexV[4] = 7; indexV[5] = 10; indexV[6] = 11; indexV[7] = 12;
	}

	public void simulateRuns (int n, RandomStream stream) { 
      	statsDeriv3.init();
		statsDeriv5.init();
		statsProb.init();
		for (int i = 0; i < n; i++)
			simulate(stream);
	}

	public void simulate(RandomStream stream)
	{
		int k;
		// Generate: 0, 1, 2, 3, 7, 10, 11, 12
		for (int j = 0; j < 8; j++)
		{
			k = indexV[j];
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
		statsProb.add(estimate);
		if (prod == 0.0) {
			statsDeriv3.add(0.0);
			statsDeriv5.add(0.0);
		}
		else {
		   if (paths[0] < paths[1])
			   statsDeriv3.add(V[2]/dist[2].getMean() * dist[5].density(x - paths[1]) 
			  	              * prod / dist[5].cdf(x - paths[1]));
		   else
			   statsDeriv3.add(0.0);
		   statsDeriv5.add (((x - paths[2]) * dist[4].density(x - paths[2]) / dist[4].getMean())
				   * (prod / dist[4].cdf(x - paths[2])));
		}
	}


	public String toString() {
		String s = "SAN network with 9 nodes and 13 links, from Elmaghraby (1977)\n"
				+ "Estimate prob longest path > x = " + x  
				+ ",\n and gradient w.r.t. theta_3 and theta_5.\n";
		return s;
	}

	public static void main(String[] args) throws IOException {
		int n = 100000;
		San13CMCgrad1 san = new San13CMCgrad1 (90.0, "san13a.dat");
		san.simulateRuns(n, new LFSR113());
		System.out.println (san.statsProb.reportAndCIStudent(0.95, 6));
		System.out.println (san.statsDeriv3.reportAndCIStudent(0.95, 6));
		System.out.println (san.statsDeriv5.reportAndCIStudent(0.95, 6));
	}
}
