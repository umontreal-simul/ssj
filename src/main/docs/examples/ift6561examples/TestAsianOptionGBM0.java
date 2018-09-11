package ift6561examples;
import umontreal.ssj.stochprocess.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.stat.*;
import umontreal.ssj.util.Chrono;

public class TestAsianOptionGBM0 {

	// For testing with a GBM process, MC only.
	public static void main(String[] args) {
		int numObsTimes = 12;
		double T1 = 1.0 / 12.0;
		double T = 1.0;
		double strike = 100.0;
		double s0 = 100.0;
		double r = 0.05;
		double sigma = 0.5;
		AsianOption asian = new AsianOption(r, numObsTimes, T1, T, strike);
		RandomStream noise = new LFSR113();
		NormalGen gen = new NormalGen(noise);
		GeometricBrownianMotion gbmSeq = new GeometricBrownianMotion(s0, r,
				sigma, new BrownianMotion(0, 0, 1, gen));
		asian.setProcess(gbmSeq);
		Tally statValueMC = new Tally("Stats on payoff with crude MC");
		int n = 1000000;

		Chrono timer = new Chrono();
		asian.simulateRuns(n, new MRG32k3a(), statValueMC);
		statValueMC.setConfidenceIntervalStudent();
		System.out.println(asian.toString());
		System.out.println(statValueMC.report(0.95, 4));
		System.out.println("Total CPU time:      " + timer.format() + "\n");
	}

}
