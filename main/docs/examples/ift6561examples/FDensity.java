package ift6561examples;
import umontreal.ssj.probdist.NormalDist;

import java.io.IOException;

import umontreal.ssj.probdist.*;
import umontreal.ssj.probdist.FisherFDist;
import umontreal.ssj.rng.RandomStream;

class FDensity {
	int n1, n2; // Degrees of freedom.
	int dim;
	double x; 
    double max;
    
    // Generates a F(d1, d2) random variable using d1 + d2 uniforms,
    // in order to estimate the density from 0 to the umax quantile.
	public FDensity (int n1, int n2, double umax) {
		this.n1 = n1;
		this.n2 = n2;
		dim = n1 + n2;
	    // double z = BetaSymmetricalDist.inverseF (1.0, 8, umax);
	    // max = ((n2 * z) / (n1 * (1 - z)));
        max = FisherFDist.inverseF (n1, n2, umax);
	    // max = 10.0;
		System.out.println("F(max) = " + cdf(max) + "\n");
		System.out.println("max = " + max + "\n\n");
	}

	// Returns the dimension (number s of required uniforms) for this model.
	public int getDimension() {
		return dim;
	};

	// Simulates the model for one run
	public double simulate(RandomStream stream)  {
		double sum1 = 0.0;
		double sum2 = 0.0;
		double z;
		for (int j = 0; j < n1; ++j) {
			z = NormalDist.inverseF01(stream.nextDouble());
			sum1 += z * z;
		}
		for (int j = 0; j < n2; ++j) {
			z = NormalDist.inverseF01(stream.nextDouble());
			sum2 += z * z;
		}
		return x = (n2 * sum1) / (n1* sum2);
		}

	// Recovers the realization of the performance measure $X$.
	public double getValue() {
		return x;
	}

	// Recovers the cdf of the last realization of the performance measure $X$.
	public double getValueU01() {
	   return FisherFDist.cdf (n1, n2, x);
	}

	// Recovers the density of X evaluated at x.
	public double density(double x) {
	   return FisherFDist.density (n1, n2, x);
	};

	// Recovers the cdf of $X$ evaluated at x.
	public double cdf(double x) {
		// return BetaDist.cdf (n1/2.0, n2/2.0, (n1*x)/(n1*x + n2));
		// return BetaSymmetricalDist.cdf (1.0, 6, (n1*x)/(n1*x + n2));
	    return FisherFDist.cdf (n1, n2, x);
	};

	// Returns a short description of the model and its parameters.
	public String toString() {
		String s = "FisherF dist.\n ";
		return s;
	}

	public double getMin() {
		return 0.0;
	}

	public double getMax() {
		return max;
	};

	public static void main(String[] args) throws IOException {
		new FDensity (2, 2, 0.9999);
	}
      
}
