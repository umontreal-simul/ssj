package ift6561examples;
import umontreal.ssj.stochprocess.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.hups.*;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.stat.*;
import umontreal.ssj.util.*;
import umontreal.ssj.mcqmctools.*;

// QMC and RQMC examples in book introduction, in 2 and 12 dimensions.
public class TestAsianOptionGBM2 {

	AsianOption asian;
	RandomStream noise = new MRG32k3a();
	PointSet pStrat;
	PointSet pSobol;
	PointSet pLCG;
	SubsetOfPointSet pLCG0;
	PointSetRandomization dShift = new LMScrambleShift(noise);
	PointSetRandomization rShift = new RandomShift(noise);
	// Next one is for the deterministic case.
	PointSetRandomization noShift = new EmptyRandomization();
	double secondsMC;
	double varianceMC;

	public TestAsianOptionGBM2(AsianOption asian) {
		this.asian = asian;
	}

	// Sobol point set will have 2^k points.
	// Korobov (LCG) point set will have n1 points.
	// Multiplier for LCG point set will be a1.
	public void createPointSets(int k, int n1, int a1) {
		int dim = asian.getNumObsTimes();
		pLCG = new LCGPointSet(n1, a1);
		pLCG0 = new SubsetOfPointSet(pLCG);
		// The LCG point set without (0,0).
		pLCG0.selectPointsRange(1, n1);
		// We create a Sobol' sequence in d-1 dimensions, the we add
		// a first coordinate equal to i/n for the point i.
		// In two dimensions, this gives the Hammersley point set.
		DigitalSequenceBase2 p0 = new SobolSequence(k, 31, dim - 1);
		pSobol = p0.toNetShiftCj();
		// if (dim <= 2) pStrat = new StratifiedUnitCube (k/2, dim);
	}

	// Make an RQMC experiment and compare with MC.
	public void experOneProcessRQMC(int m, RandomStream noise) {
		System.out.println("\n******************************************");
		System.out.println("Asian option with sequential sampling \n");
		Tally statRQMC = new Tally("Stats on payoff with RQMC");
        if (asian.getNumObsTimes() <= 2) {
 		    RQMCExperiment.simulReplicatesRQMCDefaultReportCompare(asian, pStrat, rShift,
			   	   m, statRQMC, varianceMC, secondsMC);
        }
        RQMCExperiment.simulReplicatesRQMCDefaultReportCompare(asian, pLCG, rShift,
				m, statRQMC, varianceMC, secondsMC);
		// RQMCExperiment.simulateRQMCDefaultReportCompare(asian, pLCGBaker,
		// rShift, m, noise, statRQMC, varianceMC, secondsMC);
		RQMCExperiment.simulReplicatesRQMCDefaultReportCompare(asian, pSobol,
				dShift, m, statRQMC, varianceMC, secondsMC);
		Chrono timer = new Chrono();
		RQMCExperiment.simulReplicatesRQMC(asian, pLCG0, noShift, m, statRQMC);
		System.out.println(pLCG0.toString());
		System.out.println("Total CPU time:      " + timer.format() + "\n");
		System.out.printf("Average (deterministic): %10.6g%n%n", statRQMC
				.average());
	}

	// Main program: QMC and RQMC experiment with Asian option.
	public static void main(String[] args) {
		int numObsTimes = 12;
		double T1 = 1.0 / (double) numObsTimes;
		double T = 1.0;
		double strike = 100.0;
		double s0 = 100.0;
		double r = 0.05;
		double sigma = 0.5;
		RandomStream noise = new LFSR113();
		NormalGen gen = new NormalGen(noise);
		AsianOption asian = new AsianOption(r, numObsTimes, T1, T, strike);
		asian.setProcess(new GeometricBrownianMotion(s0, r, sigma,
				new BrownianMotion(0, 0, 1, gen)));
		TestAsianOptionGBM2 test = new TestAsianOptionGBM2(asian);

		Tally statValueMC = new Tally("Stats on payoff with crude MC");
		int n = 10000000; // 10 million runs for Monte Carlo.
		Chrono timer = new Chrono();
		MonteCarloExperiment.simulateRunsDefaultReportStudent (asian, n, noise,
				statValueMC, 0.95, 4, timer);
		// We memorize CPU time and variance to compare with RQMC.
		test.secondsMC = timer.getSeconds() / (double) n;
		test.varianceMC = statValueMC.variance();

		// RQMC experiments.
		test.createPointSets(7, 101, 12);
		test.experOneProcessRQMC(1000, noise);

		test.createPointSets(10, 1024, 115);
		test.experOneProcessRQMC(1000, noise);

		test.createPointSets(16, 65536, 12421);
		test.experOneProcessRQMC(1000, noise);
		
		// test.createPointSets(20, 1048576, 286857);
		// test.experOneProcessRQMC(1000, noise);
	}

}
