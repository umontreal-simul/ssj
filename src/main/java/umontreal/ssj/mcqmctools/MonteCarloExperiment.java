package umontreal.ssj.mcqmctools;

import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.*;
import umontreal.ssj.stat.list.lincv.*;
import umontreal.ssj.util.Chrono;
import umontreal.ssj.util.PrintfFormat;

/**
 * Provides generic tools to perform simple Monte Carlo experiments with a simulation model that
 * implements one of the interfaces  
 * @ref MonteCarloModelDouble, @ref MonteCarloModelDoubleArray, or @ref MonteCarloModelCV. 
 * The experiment consists of n independent simulation runs and the results are returned in @ref Tally
 * satistical collectors. The #RandomStream is reset to a new substream after each run.
 */

public class MonteCarloExperiment {

	/**
	 * Performs #n simulation runs of #model using #stream and collects statistics in #statValue.
	 * The #stream is reset to a new substream for each run.
	 */
	public static void simulateRuns(MonteCarloModelDouble model, int n, RandomStream stream,
	        Tally statValue) {
		statValue.init();
		for (int i = 0; i < n; i++) {
			model.simulate(stream);
			statValue.add(model.getPerformance());
			stream.resetNextSubstream();
		}
	}

	/**
	 * Performs n runs of model using #stream and collects statistics for a model with a vector of
	 * control variates. The results are returned in #statWithCV.
	 */
	public static void simulateRunsCV(MonteCarloModelCV model, int n, RandomStream stream,
	        ListOfTalliesWithCV<Tally> statWithCV) {
		statWithCV.init();
		for (int i = 0; i < n; i++) {
			model.simulate(stream);
			statWithCV.add(model.getPerformance(), model.getValuesCV());
			stream.resetNextSubstream();
		}
	}

	/**
	 * Performs n runs using `stream` and collects statistics for a model with a single real-valued
	 * control variate C. The statistics on X and C are collected in `statX` and `statC`.
	 */
	public static void simulateRunsCV (MonteCarloModelCV model, int n, RandomStream stream,
	        TallyStore statX, TallyStore statC) {
		statX.init();
		statC.init();
		for (int i = 0; i < n; i++) {
			model.simulate(stream);
			statX.add(model.getPerformance());
			statC.add(model.getValuesCV()[0]);
			stream.resetNextSubstream();
		}
	}

	/**
	 * Performs n runs using #stream and collects statistics for a model with a single real-valued
	 * control variate C. The statistics are collected, and the mean and variance of the estimators
	 * with and without the control variate are returned in the two-dimensional vectors `mean` and
	 * `variance`, as in `computeMeanVarCV`.
	 */
	public static void simulateRunsCV (MonteCarloModelCV model, int n, RandomStream stream,
	        double[] mean, double[] variance) {
		TallyStore statX = new TallyStore(n);
		TallyStore statC = new TallyStore(n);
		simulateRunsCV(model, n, stream, statX, statC);
		computeMeanVarCV(statX, statC, mean, variance);
	}

	/**
	 * Given statistics collected in `statX` and `statC` as with @ref simulateRunsCV, this method
	 * computes the mean and variance of the estimators with and without the CV and returns them in
	 * the two-dimensional vectors `mean` and `variance` (mean[0] is the value without CV, mean[1] the
	 * value with CV, and similarly for the variance.
	 */
	public static void computeMeanVarCV (TallyStore statX, TallyStore statC, double[] mean,
	        double[] variance) {
		mean[0] = statX.average();
		variance[0] = statX.variance();
		double varC = statC.variance();
		double covCX = statC.covariance(statX);
		double meanC = statC.average();
		double beta = covCX / varC;
		mean[1] = mean[0] - beta * meanC;      // CV has mean 0.
		variance[1] = variance[0] + beta * beta * varC - 2 * beta * covCX;
	}

	/**
	 * Performs `n` simulation runs to estimate the difference in performance between `model2` and
	 * `model1`, divided by `delta`, using common random numbers (CRN) across the two models. One
	 * substream is used for each run and the same `n` substreams are used for the two models. Returns
	 * the statistics on the `n` differences in `statDiff`. By taking `delta` = 1, this just estimates
	 * the difference. By taking `delta` > 0 very small, and if the two models are in fact the same
	 * model but with a parameter that differs by `delta`, this gives a finite-difference estimator of
	 * the derivative of the performance with respect to this parameter.
	 */
	public static void simulFDReplicatesCRN (MonteCarloModelDouble model1,
	        MonteCarloModelDouble model2, double delta, int n, RandomStream stream,
	        Tally statDiff) {
		statDiff.init();
		double value1;
		for (int i = 0; i < n; i++) {
			stream.resetNextSubstream();
			model1.simulate(stream);
			value1 = model1.getPerformance();
			stream.resetStartSubstream();
			model2.simulate(stream);
			statDiff.add((model2.getPerformance() - value1) / delta);
		}
	}

	/**
	 * Similar to @ref simulFDReplicatesCRN, but using independent random numbers (IRN) across the
	 * two models. One substream is used for each run of each model, for a total of 2n substreams.
	 */
	public static void simulFDReplicatesIRN (MonteCarloModelDouble model1,
	        MonteCarloModelDouble model2, double delta, int n, RandomStream stream,
	        Tally statDiff) {
		statDiff.init();
		double value1;
		for (int i = 0; i < n; i++) {
			stream.resetNextSubstream();
			model1.simulate(stream);
			value1 = model1.getPerformance();
			model2.simulate(stream);
			statDiff.add((model2.getPerformance() - value1) / delta);
		}
	}

	/**
	 * Performs n independent runs using n substreams of #stream, collects statistics in #statValue,
	 * and returns a report with a confidence interval of level #level, with #d decimal fractional
	 * digits of precision for the output, computed via a Student distribution.
	 */
	public static String simulateRunsDefaultReportStudent(MonteCarloModelDouble model, int n,
	        RandomStream stream, Tally statValue, double level, int d, Chrono timer) {
		PrintfFormat str = new PrintfFormat();
		timer.init();
		simulateRuns(model, n, stream, statValue);
		statValue.setConfidenceIntervalStudent();
		str.append(model.toString() + "\n");
		str.append(statValue.report(level, d));
		// str.append("Variance per run: " + statValue.variance() + "\n");
		// str.append(7 + 5, 5, 4, statValue.variance());
		str.append("\n");
		str.append("Total CPU time:      " + timer.format() + "\n");
		return str.toString();
	}

	/**
	 * In this version, there is no need to provide a #Chrono; it is created inside.
	 */
	public static String simulateRunsDefaultReportStudent(MonteCarloModelDouble model, int n,
	        RandomStream stream, Tally statValue, double level, int d) {
		Chrono timer = new Chrono();
		return simulateRunsDefaultReportStudent(model, n, stream, statValue, level, d, timer);
	}

	/**
	 * A short-hand equivalent for `simulateRunsDefaultReportStudent (model, n, stream, statValue, 0.95, 4)`
	 */
	public static String simulateRunsDefaultReport (MonteCarloModelDouble model, int n,
	        RandomStream stream, Tally statValue) {
		Chrono timer = new Chrono();
		return simulateRunsDefaultReportStudent (model, n, stream, statValue, 0.95, 4, timer);
	}

	/**
	 * Similar to @ref simulateRunsDefaultReport, but this one uses a vector of control variates. It
	 * returns a report with a confidence interval for the estimator with the CV.
	 */
	public static String simulateRunsDefaultReportCV (MonteCarloModelCV model, int n,
	        RandomStream stream, ListOfTalliesWithCV<Tally> statWithCV, double level, int d,
	        Chrono timer) {
		PrintfFormat str = new PrintfFormat();
		timer.init();
		simulateRunsCV(model, n, stream, statWithCV);
		statWithCV.estimateBeta();  // Computes the variances and covariances!
		// statWithCV.setConfidenceIntervalStudent();
		str.append(model.toString() + "\n");
		// str.append(statWithCV.report(0.95, 4));
		double[] centerAndRadius = new double[2];
		statWithCV.confidenceIntervalStudentWithCV(0, level, centerAndRadius);
		str.append("Average: " + statWithCV.averageWithCV(0) + "\n");
		str.append("Variance per run (Cov[0,0]): " + statWithCV.covarianceWithCV(0, 0) + "\n");
		double[] varCV = new double[2];
		statWithCV.varianceWithCV(varCV);
		str.append("Variance per run with CV: " + varCV[0] + "\n");
		str.append("Center of CI:  " + centerAndRadius[0] + "\n");
		str.append("Radius of CI:  " + centerAndRadius[1] + "\n");
		str.append("Total CPU time:     " + timer.format() + "\n");
		return str.toString();
	}

	/**
	 * This one uses a single real-valued CV, as in @ref simulateRunsCV.
	 */
	public static String simulateRunsDefaultReportCV (MonteCarloModelCV model, int n,
	        RandomStream stream, double[] mean, double[] variance, double level, int d,
	        Chrono timer) {
		PrintfFormat str = new PrintfFormat();
		timer.init();
		simulateRunsCV(model, n, stream, mean, variance);
		// statWithCV.setConfidenceIntervalStudent();
		str.append(model.toString() + "\n");
		// str.append(statWithCV.report(0.95, 4));
		str.append("Average, no CV:  " + mean[0] + "\n");
		str.append("Average with CV:  " + mean[1] + "\n");
		str.append("Variance, no CV:  " + variance[0] + "\n");
		str.append("Variance with CV:  " + variance[1] + "\n");
		str.append("Total CPU time:      " + timer.format() + "\n");
		return str.toString();
	}

}
