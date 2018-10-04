package ift6561examples;

import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.mcqmctools.*;

// This is question 4 of homework 2 of ift6561, A-2018.

public class ProductBarrier implements MonteCarloModelDouble {

	double C, K, a, b;
	double mu1, sigma1;
	double mu2, sigma2;
	ContinuousDistribution dist1;  // For W_1
	ContinuousDistribution dist2;  // For W_2
	double payoff;   // Value of X to return.

	// Constructor.
	public ProductBarrier(double C, double K, double a, double b, double mu1, double sigma1,
	        double mu2, double sigma2) {
		this.C = C;
		this.K = K;
		this.a = a;
		this.b = b;
		dist1 = new LognormalDist(mu1, sigma1);
		dist2 = new LognormalDist(mu2, sigma2);
	}

	// Generates payoff X, without IS.
	public void simulate(RandomStream stream) {
		payoff = 0.0;
		double CW1 = C * dist1.inverseF(stream.nextDouble());
		if (CW1 > a)
			return;
		double X = CW1 * dist2.inverseF(stream.nextDouble());
		if (X > b)
			payoff = X - K;
	}

	// Returns payoff X
	public double getPerformance() {
		return payoff;
	}

	// Descriptor of model
	public String toString() {
		return "Simplified financial option with barriers, Devoir 2, A-2018, no IS";
	}
	
}
