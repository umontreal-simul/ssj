package ift6561examples;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.mcqmctools.*;

// This was intended to be a faster implementation for the special case
// of a GMB process, but it turns out that the more general
// AsianOption with sp = GeometricBrownianMotion  is at least as fast!
//
public class OldAsianGBM implements MonteCarloModelDouble {
	double strike; // Strike price.
	int d; // Number of observation times.
	double discount; // Discount factor exp(-r * zeta[t]).
	double[] muDelta; // Differences * (r - sigma^2/2).
	double[] sigmaSqrtDelta; // Square roots of differences * sigma.
	double[] logS; // Log of the GBM process: logS[t] = log (S[t]).

	// double[] normals; // To store the standard normals.

	// obsTimes[0..d+1] must contain obsTimes[0]=0.0, plus the s observation
	// times.
	public OldAsianGBM(double r, double sigma, double strike, double s0, int d,
			double[] obsTimes) {
		this.strike = strike;
		this.d = d;
		// normals = new double[d];
		discount = Math.exp(-r * obsTimes[d]);
		double mu = r - 0.5 * sigma * sigma;
		muDelta = new double[d];
		sigmaSqrtDelta = new double[d];
		logS = new double[d + 1];
		double delta;
		for (int j = 0; j < d; j++) {
			delta = obsTimes[j + 1] - obsTimes[j];
			muDelta[j] = mu * delta;
			sigmaSqrtDelta[j] = sigma * Math.sqrt(delta);
		}
		logS[0] = Math.log(s0);
	}

	// Generates the process S.
	public void simulate(RandomStream stream) {
		// gen.nextArrayOfDouble(normals, 0, d);
		for (int j = 0; j < d; j++)
			logS[j + 1] = logS[j] + muDelta[j] + sigmaSqrtDelta[j]
					* NormalDist.inverseF01(stream.nextDouble());
	}

	// Computes and returns the discounted option payoff.
	public double getPerformance() {
		double average = 0.0; // Average of the GBM process.
		for (int j = 1; j <= d; j++)
			average += Math.exp(logS[j]);
		average /= d;
		if (average > strike)
			return discount * (average - strike);
		else
			return 0.0;
	}

	public static void main(String[] args) {
		int d = 12;
		double[] obsTimes = new double[d + 1];
		obsTimes[0] = 0.0;
		for (int j = 1; j <= d; j++)
			obsTimes[j] = (double) j / (double) d;
		OldAsianGBM asian = new OldAsianGBM(0.05, 0.5, 100.0, 100.0, d, obsTimes);
		Tally statValue = new Tally("Stats on value of Asian option");

		int n = 1000000;
		MonteCarloExperiment.simulateRunsDefaultReportStudent (asian, n,
				new MRG32k3a(), statValue, 0.95, 4);
	}
}
