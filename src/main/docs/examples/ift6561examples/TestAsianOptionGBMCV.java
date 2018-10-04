package ift6561examples;

import java.io.*;
import umontreal.ssj.stochprocess.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.stat.*;
import umontreal.ssj.stat.list.lincv.ListOfTalliesWithCV;
import umontreal.ssj.charts.HistogramChart;
import umontreal.ssj.mcqmctools.*;
import umontreal.ssj.util.Chrono;

public class TestAsianOptionGBMCV {

	// Testing experiments for an Asian Option under a GBM process, with and without a control variate.
	// The program first makes an MC experiment to estimate the distribution of the payoff without the CV.
	// It computes a confidence interval on the mean, then makes a histogram of the positive payoffs.
	// It also produces a standalone Latex file with the histogram.
	// After that, it makes another experiment with the CV, and produces a report that compares the variances 
	// and gives a 95% confidence interval.
	public static void main(String[] args) throws IOException {
		int numObsTimes = 12;
		double T1 = 1.0 / 12.0;
		double T = 1.0;
		double strike = 100.0;
		double s0 = 100.0;
		double r = 0.05;
		double sigma = 0.5;
		int n = 1000000;   // Number of simulation runs.

		AsianOptionGBMCV asian = new AsianOptionGBMCV (r, numObsTimes, T1, T, strike);
		RandomStream noise = new LFSR113();
		NormalGen gen = new NormalGen(noise);
		GeometricBrownianMotion gbmSeq = new GeometricBrownianMotion(s0, r,
				sigma, new BrownianMotion(0, 0, 1, gen));
		asian.setProcess(gbmSeq);
		TallyStore statValueMC = new TallyStore ("Stats on payoff with crude MC");

		Chrono timer = new Chrono();
		System.out.println (MonteCarloExperiment.simulateRunsDefaultReport 
				(asian, n, noise, statValueMC));

		// Extract positive payoffs, put them in collector statValuePosMC, and print report.
		TallyStore statValuePosMC = statValueMC.extractSubrange (0.00000000001, 1.0E200);
		statValuePosMC.setName ("Stats on positive payoffs with crude MC");
		System.out.println(statValuePosMC.report(0.95, 4));
		double fractionZero = 1.0 - (double) statValuePosMC.numberObs()
				/ (double) n;
		System.out.printf("Proportion of zero payoffs: %12.6f%n", fractionZero);

		// Make a histogram of positive discounted payoffs.
		HistogramChart hist = new HistogramChart(
				"Distribution of positive discounted payoff", "Payoff",
				"Frequency", statValuePosMC.getArray(), statValuePosMC.numberObs());
		double[] bounds = { 0, 150, 0, 35000 };  // Range for x and y.
		hist.setManualRange(bounds);
		(hist.getSeriesCollection()).setBins(0, 60, 0, 150); // 60 bins over [0, 150].
		hist.view(800, 500);
		hist.toLatexFile ("asianHist.tex", 10.0, 8.0);  // Stand-alone Latex file.

		// double[] mean = new double[2];
		// double[] variance = new double[2];
		// System.out.println (MonteCarloExperiment.simulateRunsDefaultReportCV 
		//		(asian, n, noise, mean, variance, 0.95, 4, timer));
		
		ListOfTalliesWithCV<Tally> list = ListOfTalliesWithCV.createWithTally (1, 1);
		System.out.println (MonteCarloExperiment.simulateRunsDefaultReportCV 
				(asian, n, noise, list, 0.95, 4, timer));
	}

}
