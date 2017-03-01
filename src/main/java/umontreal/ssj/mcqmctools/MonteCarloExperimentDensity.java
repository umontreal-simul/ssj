package umontreal.ssj.mcqmctools;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.stat.*;
import umontreal.ssj.util.Chrono;

/**
 * Provides generic tools to perform Monte Carlo experiments
 * with a simulation model that implements the MCModelDensity interface,
 * with MC or RQMC.
 */

/**
 * @author Pierre L'Ecuyer
 * 
 */
public class MonteCarloExperimentDensity extends MonteCarloExperiment {

	// Performs n runs using stream and collects statistics in statValue.
	// Note that the stream can be an RQMC stream.
	// Must override because here we must have a TallyStore.
	public static void simulateRuns (MonteCarloModelBounded model, int n,
			RandomStream stream, TallyStore statValue) {
		statValue.init();
		for (int i = 0; i < n; i++) {
			model.simulate(stream);
			statValue.add (model.getPerformance());
			stream.resetNextSubstream();
		}
	}

	// Performs n runs using stream and collects statistics in statValue.
	// Note that the stream can be an RQMC stream.
	public static void simulAndPlotHistogram (MonteCarloModelBounded model, int n, 
			 double a, double b, int numBins) {	
		TallyStore statValue = new TallyStore();
		TallyHistogram hist = new TallyHistogram(a, b, numBins);
		RandomStream stream = new MRG32k3a();
		
		simulateRuns (model, n, stream, statValue); 
		hist.fillFromArray(statValue.getArray());
		ScaledHistogram histDensity = new ScaledHistogram(hist, 1.0);
		// ScaledHistogram histAsh = histDensity.averageShiftedHistogramTrunc (bandwidthASH);
		// histAsh.rescale(1.0);
		HistogramChartToLatex histLatex = new HistogramChartToLatex();
		System.out.println(histLatex.toLatex(histDensity, true, true));
		// histLatex.writeStringTofile ("histmc", histLatex.toLatex(histDensity, true, true));
	}

	// Performs n indep. runs using stream, collects statistics in statValue,
	// and prints a report with a confidence interval.
//	public static void simulateRunsDefaultReport(MonteCarloModelBounded model, int n,
//			RandomStream stream, TallyStore statValue) {
//		Chrono timer = new Chrono();
//		simulateRuns(model, n, stream, statValue);
//		statValue.setConfidenceIntervalStudent();
//		System.out.println(model.toString());
//		System.out.println(statValue.report(0.95, 4));
//		System.out
//				.printf("Variance per run: %9.5g%n", statValue.variance() * n);
//		System.out.println("Total CPU time:      " + timer.format() + "\n");
//	}
	
}
