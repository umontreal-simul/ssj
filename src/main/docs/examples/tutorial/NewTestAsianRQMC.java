package tutorial;

import java.io.*;
import java.util.ArrayList;
import umontreal.ssj.rng.*;
import umontreal.ssj.hups.*;
import umontreal.ssj.mcqmctools.*;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.PgfDataTable;

public class NewTestAsianRQMC extends AsianGBM implements MonteCarloModelDouble {

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

	public int getDimension() {
		return d;
	}

	public String toString() {
		return "Asian option under GBM, for testing";
	}
	
	public static void main(String[] args) throws IOException {
	// public static void main() {
		int s = 12;   // Dimension (number of time steps).
		int dim = s;
		double[] zeta = new double[s + 1];  // The observation times.
		for (int j = 0; j <= s; j++)
			zeta[j] = (double) j / (double) s;
		NewTestAsianRQMC model = new NewTestAsianRQMC(0.05, 0.5, 100.0, 100.0, s, zeta);
		Tally statValue = new Tally("value of Asian option");
		// Tally statQMC = new Tally("QMC averages for Asian option");

		// RQMCExperimentDouble exper = new RQMCExperimentDouble();
		int n = 10000;
		System.out.println("Ordinary MC:\n");
		System.out.println (MonteCarloExperiment.simulateRunsDefaultReportStudent
				(model, n, new MRG32k3a(), statValue, 0.95, s));

		int base = 2;        // Basis for the loglog plots.
		int numSets = 10;    // Number of sets in the series.
		int numSkipReg = 1;  // Number of sets skipped for the regression.
		int i;
		int m = 50;          // Number of RQMC randomizations.
		int[] N = { 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288,
		        1048576, 2097152 }; // 13
		int[] a = { 149, 275, 857, 1731, 2431, 6915, 12545, 19463, 50673, 96407, 204843, 443165,
		        768165 };
		
        // Create a list of series of RQMC point sets.
		ArrayList<RQMCPointSet[]> listOfSeries = new ArrayList<RQMCPointSet[]>();
		PointSet p; 
		PointSetRandomization rand;
		RandomStream noise = new MRG32k3a();

		// Korobov lattice + baker
		RQMCPointSet[] theRQMCSetLattice = new RQMCPointSet[numSets];
		rand = new RandomShift(noise);
		for (i = 0; i < numSets; ++i) {
			p = new BakerTransformedPointSet(new KorobovLattice(N[i], a[i], dim));
			theRQMCSetLattice[i] = new RQMCPointSet(p, rand);
		}
		theRQMCSetLattice[0].setLabel("Lattice+baker+shift");
		listOfSeries.add(theRQMCSetLattice);

		// Sobol + LMS + shift
		RQMCPointSet[] theRQMCSetSobol = new RQMCPointSet[numSets];
		int mink = 9;     // Smallest power of 2 considered.
		rand = new LMScrambleShift(noise);
		for (i = 0; i < numSets; ++i) {
			p = (new SobolSequence(i + mink, 31, dim)); 
			theRQMCSetSobol[i] = new RQMCPointSet(p, rand);
		}
		theRQMCSetSobol[0].setLabel("Sobol+LMS+Shift");
		listOfSeries.add(theRQMCSetSobol);
		
		RQMCExperimentSeries experSeries = new RQMCExperimentSeries (theRQMCSetLattice, base);
		experSeries.testVarianceRate (model, m); 
		System.out.println (experSeries.reportVarianceRate (numSkipReg, true)); 
		System.out.println ((experSeries.toPgfDataTable ("Korobov+baker")).drawPgfPlotSingleCurve 
		      ("Korobov+baker", "axis", 3, 4, 2, "", "marks=*"));

		// Perform an experiment with a list of series of RQMC point sets.
		ArrayList<PgfDataTable> listCurves = new ArrayList<PgfDataTable>();
		System.out.println (experSeries.testVarianceRateManyPointTypes (model, listOfSeries, m, numSkipReg,
				true, true, true, listCurves));
		System.out.println ("\n Now printing data table for the two curves  *****  \n\n\n");
		// Prints the data of each curve as a table.		
		for (PgfDataTable curve : listCurves)
			System.out.println (curve.formatTable());
		// Produces LaTeX code to draw these curves with pgfplot. 
		String plot = PgfDataTable.drawPgfPlotManyCurves 
				("Korobov and Sobol", "loglogaxis", 0, 2, listCurves, 2, "", " ");
	    System.out.println (plot);

	    // Produces a complete LaTeX file with the plots.
	    String pfile = (PgfDataTable.pgfplotFileHeader () + plot + 
   	    		PgfDataTable.pgfplotEndDocument ());
		Writer file = new FileWriter("testPlotAsianRQMC.tex");
		file.write (pfile);
		file.close();
	}

}
