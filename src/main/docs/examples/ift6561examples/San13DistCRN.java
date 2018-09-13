package ift6561examples;
import java.io.*;

import umontreal.ssj.charts.HistogramChart;
import umontreal.ssj.rng.*;
import umontreal.ssj.stat.TallyStore;
import umontreal.ssj.util.Chrono;


/**
 * Here we compute the empirical distribution of the shortest path lengths and
 * we construct a histogram.
 */

public class San13DistCRN extends San13Dist {

	// The constructor reads link length distributions in a file.
	public San13DistCRN(String fileName) throws IOException {
		super(fileName);
	}

	public void updatePaths() {
		V[2] *= 10.0/7.0;
		V[4] *= 18.5/16.5;
		// Path lengths
		paths[1] = V[0] + V[2] + V[5] + V[10];
		paths[2] = V[0] + V[4] + V[10];
		if (paths[1] > maxPath)
				maxPath = paths[1];
		if (paths[2] > maxPath)
			maxPath = paths[2];
		// return maxPath;
	}
	
	public static void main(String[] args) throws IOException {
		int n = 100000;
		San13DistCRN san = new San13DistCRN("san13a.dat");
		TallyStore statIRN = new TallyStore("Tally for IRN");
		TallyStore statCRN = new TallyStore("Tally for CRN");
		TallyStore statCRNPos = new TallyStore("Tally for CRN, positive values");
		RandomStream stream = new LFSR113();
		double x1, delta;
		Chrono timer = new Chrono();
	    for (int i=0; i<n; i++) {
	         stream.resetNextSubstream();
	    	 san.simulate(stream);
	    	 x1 = san.getPerformance();
	    	 san.updatePaths();
		     delta = san.getPerformance() - x1;
	    	 statCRN.add(delta);
			 if (delta > 0.0000000001)
				statCRNPos.add(delta);
	    	 
	    	 san.simulate(stream);  // with new random numbers.
	    	 san.updatePaths();     // and new parameters.
	    	 statIRN.add(san.getPerformance()-x1);
	    }
		System.out.println("Total CPU time:      " + timer.format() + "\n");
		statCRN.setConfidenceIntervalStudent();
		System.out.println(statCRN.report(0.95, 4));
		System.out
				.printf("Variance per run: %9.5g%n", statCRN.variance() * n);
		System.out.println(statCRNPos.report(0.95, 4));
		statIRN.setConfidenceIntervalStudent();
		System.out.println(statIRN.report(0.95, 4));
		System.out
				.printf("Variance per run: %9.5g%n", statIRN.variance() * n);
		System.out
		.printf("Variance ratio: %9.5g%n", statIRN.variance() / statCRN.variance());

		HistogramChart hist = new HistogramChart("Dist. of positive Delta, CRN",
				"Values of Delta", "Frequency", statCRNPos.getArray(), statCRNPos.numberObs());
		double[] bounds = { 0, 20, 0, 2200 };
		hist.setManualRange(bounds);
		(hist.getSeriesCollection()).setBins(0, 40, 0, 20);
		hist.view(800, 500);
		String histLatex = hist.toLatex(12.0, 8.0);
		Writer file = new FileWriter("san13CRNchart.tex");
		file.write(histLatex);
		file.close();
		
		hist = new HistogramChart("Distribution of Delta IRN",
				"Values of Delta", "Frequency", statIRN.getArray(), n);
		double[] bounds2 = { -150, 150, 0, 16000 };
		hist.setManualRange(bounds2);
		(hist.getSeriesCollection()).setBins(0, 30, -150, 150);
		hist.view(800, 500);
		histLatex = hist.toLatex(12.0, 8.0);
		file = new FileWriter("san13IRNchart.tex");
		file.write(histLatex);
		file.close();
	}
}
