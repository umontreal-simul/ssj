package tutorial;

import umontreal.ssj.rng.*;

import java.util.ArrayList;

import umontreal.ssj.hups.*;
import umontreal.ssj.mcqmctools.*;
import umontreal.ssj.stat.Tally;

public class NewTestAsianRQMC extends Asian implements MonteCarloModelDouble {

	public NewTestAsianRQMC(double r, double sigma, double strike, double s0, int s,
	        double[] zeta) {
		super(r, sigma, strike, s0, s, zeta);
	}

	@Override
	public void simulate(RandomStream stream) {
		generatePath(stream);
	}

	@Override
	public double getPerformance() {
		return getPayoff();
	}

	@Override
	public int getDimension() {
		return s;
	}

	public static void main(String[] args) {
		int s = 12;
		int dim = s;
		double[] zeta = new double[s + 1];
		for (int j = 0; j <= s; j++)
			zeta[j] = (double) j / (double) s;
		NewTestAsianRQMC model = new NewTestAsianRQMC(0.05, 0.5, 100.0, 100.0, s, zeta);
		Tally statValue = new Tally("value of Asian option");
		// Tally statQMC = new Tally("QMC averages for Asian option");

		// RQMCExperimentDouble exper = new RQMCExperimentDouble();

		int n = 10000;
		System.out.println("Ordinary MC:\n");
		RQMCExperimentDouble.simulateRunsDefaultReportStudent(model, n, new MRG32k3a(), statValue,
		        0.95, s);

		int base = 2;
		int numSets = 8;
		int numSkipReg = 1;
		int i;
		int m = 20;                     // Number of RQMC randomizations.
		PointSet p;  // [] pointSet = new PointSet[numSets];
		RQMCPointSet[] theRQMCSet = new RQMCPointSet[numSets];
		RandomStream noise = new MRG32k3a();
		PointSetRandomization rand;
		int[] N = { 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288,
		        1048576, 2097152 }; // 13
		int[] a = { 149, 275, 857, 1731, 2431, 6915, 12545, 19463, 50673, 96407, 204843, 443165,
		        768165 };
		int mink = 9;

		RQMCExperimentSeries experSeries = new RQMCExperimentSeries(theRQMCSet, base);
		ArrayList<RQMCPointSet[]> listRQMC = new ArrayList<RQMCPointSet[]>();

		// Korobov lattice + baker
		rand = new RandomShift(noise);
		for (i = 0; i < numSets; ++i) {
			p = new BakerTransformedPointSet(new KorobovLattice(N[i], a[i], dim));
			theRQMCSet[i] = new RQMCPointSet(p, rand);
		}
		listRQMC.add(theRQMCSet);

		experSeries.testVarianceRate (model, m); 
		System.out.println (experSeries.reportVarianceRate (numSkipReg, true)); 
		System.out.println ((experSeries.toPgfDataTable ("Korobov+baker")).drawPgfPlotSingleCurve 
		      ("Korobov+baker", 0, 2, true, "no marks"));

		// Sobol + LMS + shift
		rand = new LMScrambleShift(noise);
		for (i = 0; i < numSets; ++i) {
			p = (new SobolSequence(i + mink, 31, dim)).toNetShiftCj();
			theRQMCSet[i] = new RQMCPointSet(p, rand);
		}
		listRQMC.add(theRQMCSet);
		
		ArrayList<PgfDataTable> listCurves = null;
		System.out.println (experSeries.testVarianceRateManyPointTypes (model, listRQMC, m, numSkipReg,
				true, true, true, listCurves));
		System.out.println (PgfDataTable.drawPgfPlotManyCurves ("two curves", 0, 2, listCurves, true, "no marks"));
		
	}

}
