package ift6561examples;

import umontreal.ssj.stochprocess.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.stat.*;
import umontreal.ssj.mcqmctools.*;
import umontreal.ssj.util.Chrono;

public class TestAsianOptionGBMCV {

	// Testing experiments for an Asian Option under a GBM process, with and without a control variate.
	public static void main(String[] args) {
		int numObsTimes = 12;
		double T1 = 1.0 / 12.0;
		double T = 1.0;
		double strike = 100.0;
		double s0 = 100.0;
		double r = 0.05;
		double sigma = 0.5;
		AsianOptionGBMCV asian = new AsianOptionGBMCV (r, numObsTimes, T1, T, strike);
		RandomStream noise = new LFSR113();
		NormalGen gen = new NormalGen(noise);
		GeometricBrownianMotion gbmSeq = new GeometricBrownianMotion(s0, r,
				sigma, new BrownianMotion(0, 0, 1, gen));
		asian.setProcess(gbmSeq);
		TallyStore statValueMC = new TallyStore ("Stats on payoff with crude MC");
		// TallyStore statValueCV = new TallyStore ("Stats on values of CV");

		int n = 1000000;

		Chrono timer = new Chrono();
		System.out.println (MonteCarloExperiment.simulateRunsDefaultReport (asian, n, noise, statValueMC));

		// Extract positive payoffs and put them in collector statValuePosMC.
		TallyStore statValuePosMC = statValueMC.extractSubrange (0.00000000001, 1.0E200);
		statValuePosMC.setName ("Stats on positive payoffs with crude MC");
		System.out.println(statValuePosMC.report(0.95, 4));

		double[] mean = new double[2];
		double[] variance = new double[2];
		System.out.println (MonteCarloExperiment.simulateRunsDefaultReportCV 
				(asian, n, noise, mean, variance, 0.95, 4, timer));
	}

}
