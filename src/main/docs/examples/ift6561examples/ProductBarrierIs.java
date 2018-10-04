package ift6561examples;

import java.io.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.mcqmctools.MonteCarloExperiment;
import umontreal.ssj.stat.*;

// This is question 4 of homework 2 of ift6561, A-2018.

public class ProductBarrierIs extends ProductBarrier {

	double maxU1;    // U1 is generated in [0, maxU1)

	// Constructor.
	public ProductBarrierIs(double C, double K, double a, double b, double mu1, double sigma1,
	        double mu2, double sigma2) {
		super(C, K, a, b, mu1, sigma1, mu2, sigma2);
		maxU1 = dist1.cdf(a / C); 
	}

	// Generates and returns X, without IS.
	public void simulate(RandomStream stream) {
		double CW1 = C * dist1.inverseF(maxU1 * stream.nextDouble());
		double minU2 = dist2.cdf(b / CW1);  // U2 must be larger than this.
		double u2 = minU2 + (1.0 - minU2) * stream.nextDouble();
		double X = CW1 * dist2.inverseF(u2);
		payoff = (X - K) * maxU1 * (1.0 - minU2);  // Unbiased IS estimator.
	}

	// Descriptor of model
	public String toString() {
		return "Simplified financial option with barriers, with IS";
	}
	
	public static void main(String[] args) throws IOException {
		double C = 100, K = 102, a = 100, b = 102;
		double mu1 = 0.01, sigma1 = 0.05;
		double mu2 = 0.01, sigma2 = 0.05;
		int n = 100000;
		RandomStream stream = new LFSR113();
		Tally statX = new TallyStore("Option payoffs");  // To store the n observations of X.

		ProductBarrier pb = new ProductBarrier(C, K, a, b, mu1, sigma1, mu2, sigma2);
		System.out.println (MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX));
		pb = new ProductBarrierIs(C, K, a, b, mu1, sigma1, mu2, sigma2);
		System.out.println (MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX));

		b = 112;
		pb = new ProductBarrier(C, K, a, b, mu1, sigma1, mu2, sigma2);
		System.out.println (MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX));
		pb = new ProductBarrierIs(C, K, a, b, mu1, sigma1, mu2, sigma2);
		System.out.println (MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX));
    }
}
