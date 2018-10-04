package ift6561examples;

import java.io.*;
import umontreal.ssj.hups.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.mcqmctools.*;

// This is question 3 of homework 2 of ift6561, A-2018.

public class ProductExpCosRQMC implements MonteCarloModelDouble {

	int s;
	double a, b;
	double prod;

	// Constructor.
	public ProductExpCosRQMC(int s, double a, double b) {
		this.s = s; 	this.a = a;		this.b = b;
	}

	// Generates and returns X, without IS.
	public void simulate (RandomStream stream) {
		prod = 1.0;
		double u;
		for (int j = 0; j < s; j++) {
			u = stream.nextDouble();
			prod *= Math.exp(a * u) * Math.cos(b * u);
		}
	}

	// Generates and returns X, without IS.
	public double getPerformance () {
		return prod;
	}

	// Descriptor of this model.
	public String toString () {
		return "Test function for MC and RQMC: product of exponentials and cosine functions.";
	}

	public static void main(String[] args) throws IOException {
		int s = 3;
		int n = 100000;
		int m = 20;                     // Number of RQMC randomizations.
		RandomStream stream = new LFSR113();
		DigitalNet p = new SobolSequence(16, 31, s); // n = 2^{16} points in s dim.
		// PointSetRandomization rand = new LMScrambleShift(stream);
		PointSetRandomization rand = new RandomShift(stream);

		System.out.println (RQMCExperiment.makeComparisonExperimentMCvsRQMC 
			    (new ProductExpCosRQMC(s, 2.0, 0.5), stream, p, rand, n, m));
		System.out.println (RQMCExperiment.makeComparisonExperimentMCvsRQMC 
				(new ProductExpCosRQMC(s, 2.0, 50.0), stream, p, rand, n, m));
		System.out.println (RQMCExperiment.makeComparisonExperimentMCvsRQMC 
				(new ProductExpCosRQMC(s, 10.0, 0.5), stream, p, rand, n, m));
	}
}
