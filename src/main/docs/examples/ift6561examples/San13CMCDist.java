package ift6561examples;
import java.io.*;

import umontreal.ssj.charts.HistogramChart;
import umontreal.ssj.rng.*;
import umontreal.ssj.stat.TallyStore;
import umontreal.ssj.mcqmctools.*;


/**
 * Here we compute the empirical distribution of the CMC estimator.   
 *      MAKES NO SENSE !!!!!!!!!!!
 */

public class San13CMCDist extends San13CMC {

	// The constructor reads link length distributions in a file.
	public San13CMCDist(String fileName) throws IOException {
		super(0.0, fileName);
	}

	public double getValue() {
		return maxPath;
	}

	public String toString() {
		return "This one makes no sense !!!!\n\n"
				+ "SAN network with 9 nodes and 13 links, from Elmaghraby (1977)\n"
				+ "Estimate distribution of length of longest path.\n";
	}

	public static void main(String[] args) throws IOException {
		int n = 100000;
		San13CMCDist san = new San13CMCDist("san13a.dat");
		TallyStore statT = new TallyStore("TallyStore for SAN13CMC example");
		MonteCarloExperiment.simulateRunsDefaultReportStudent (san, n, new LFSR113(),
				statT, 0.95, 4);
		statT.quickSort();
		HistogramChart hist = new HistogramChart("Distribution of $T$",
				"Values of $T$", "Frequency", statT.getArray(), n);
		double[] bounds = { 0, 200, 0, 12000 };
		hist.setManualRange(bounds);
		(hist.getSeriesCollection()).setBins(0, 40, 0, 200);
		hist.view(800, 500);
		String histLatex = hist.toLatex(12.0, 8.0);
		Writer file = new FileWriter("san13CMCchart.tex");
		file.write(histLatex);
		file.close();

		// Print p-th quantile
		double p = 0.99;
		int index = (int)Math.round (p * n);
		double xip = statT.getArray()[index];
		System.out.printf("%5.3g -th quantile: %9.6g \n", p, xip);
	}
}
