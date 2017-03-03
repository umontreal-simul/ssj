package umontreal.ssj.mcqmctools;

import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.*;
import umontreal.ssj.stat.list.lincv.*;
import umontreal.ssj.util.Chrono;

/**
 * Provides generic tools to perform Monte Carlo experiments
 * with a simulation model that implements the MonteCarloModelDouble interface.
 */

/**
 * @author Pierre L'Ecuyer
 * 
 */
public class MonteCarloExperiment {

	// Performs n runs using stream and collects statistics in statValue.
	public static void simulateRuns (MonteCarloModelDouble model, int n,
			RandomStream stream, Tally statValue) {
		statValue.init();
		for (int i = 0; i < n; i++) {
			model.simulate(stream);
			statValue.add(model.getPerformance());
			stream.resetNextSubstream();
		}
	}
	
	/**
	 * Performs <SPAN CLASS="MATH"><I>n</I></SPAN> runs using <TT>stream</TT> and
	 * collects statistics in collectors <TT>statValue</TT> for the performance and 
	 * <TT>statValuesCV</TT> for the control variates.
	 */
	public static void simulateRunsCV (MonteCarloModelCV model, int n, RandomStream stream, 
			  ListOfTalliesWithCV<Tally> statWithCV) {
		// Performs n runs using stream, and collects statistics in
		// statValue for X and in statValuesCV for the CVs.
		statWithCV.init();
		for (int i = 0; i < n; i++) {
			model.simulate(stream);
			statWithCV.add(model.getPerformance(), model.getValuesCV());
			stream.resetNextSubstream();
		}
	}
	
	// Simulation with Monte Carlo with and without CV,  different version....
	public static void simulateRunsCV1 (MonteCarloModelCV model, int n, RandomStream stream, 
			   double[] meanPayoff, double[] varPayoff) 
	   {
	      TallyStore statX = new TallyStore(n);
	      TallyStore statC = new TallyStore(n);
			for (int i = 0; i < n; i++) {
				model.simulate(stream);
				statX.add(model.getPerformance());
				statC.add(model.getValuesCV()[0]);
				stream.resetNextSubstream();
			}
	      computeMeanVar (statX, statC, meanPayoff, varPayoff);
	   }

	
	   // Computes empirical mean and variance for the 
	   // uncontrolled and controlled mean estimators.
	public static void computeMeanVar (TallyStore statX, TallyStore statC, 
			    double[] mean, double[] variance) {
	      mean[0] = statX.average();
	      variance[0] = statX.variance();
	      double varC = statC.variance();
	      double covCX = statC.covariance (statX);
	      double meanC = statC.average();
	      double beta = covCX / varC;
	      mean[1] = mean[0] - beta * meanC;      // CV has mean 0.
	      variance[1] = variance[0] + beta * beta * varC - 2 * beta * covCX;
	   }
	   


    // To estimate a difference between the performance of two models, with CRNs.
	public static void simulFDReplicatesCRN(MonteCarloModelDouble model1, MonteCarloModelDouble model2, double delta,
			int n, RandomStream stream, Tally statDiff) {
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
	
    // To estimate a difference between the performance of two models, with IRNs.
	public static void simulFDReplicatesIRN(MonteCarloModelDouble model1, MonteCarloModelDouble model2, double delta,
			int n, RandomStream stream, Tally statDiff) {
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
	
	// Performs n indep. runs using stream, collects statistics in statValue,
	// and prints a report with a confidence interval.
	public static void simulateRunsDefaultReport(MonteCarloModelDouble model, int n,
			RandomStream stream, Tally statValue) {
		Chrono timer = new Chrono();
		simulateRunsDefaultReport (model, n, stream, statValue, timer);
	}

	public static void simulateRunsDefaultReport (MonteCarloModelDouble model, int n,
			RandomStream stream, Tally statValue, Chrono timer) {
		timer.init();
		simulateRuns(model, n, stream, statValue);
		statValue.setConfidenceIntervalStudent();
		System.out.println(model.toString());
		System.out.println(statValue.report(0.95, 4));
		System.out
				.printf("Variance per run: %9.5g%n", statValue.variance());
		System.out.println("Total CPU time:      " + timer.format() + "\n");
	}
	
	public static void simulateRunsDefaultReportCV (MonteCarloModelCV model, int n,
			RandomStream stream, ListOfTalliesWithCV<Tally> statWithCV, Chrono timer) {
		timer.init();
		simulateRunsCV (model, n, stream, statWithCV);
		statWithCV.estimateBeta();  // Computes the variances and covariances!
		// statWithCV.setConfidenceIntervalStudent();
		System.out.println(model.toString());
		// System.out.println(statWithCV.report(0.95, 4));
		double[] centerAndRadius = new double[2];
        statWithCV.confidenceIntervalStudentWithCV (0, 0.95, centerAndRadius);
		System.out.printf("Average:   %9.5g%n", statWithCV.averageWithCV(0));
		System.out.printf("Variance per run (Cov[0,0]): %9.5g%n", statWithCV.covarianceWithCV(0, 0));
		double[] varCV = new double[2];
		statWithCV.varianceWithCV(varCV);
		System.out.printf("Variance per run with CV: %9.5g%n", varCV[0]);
		System.out.printf("Center of CI:   %9.5g%n", centerAndRadius[0]);
		System.out.printf("Radius of CI:   %9.5g%n", centerAndRadius[1]);
		System.out.println("Total CPU time:      " + timer.format() + "\n");
	}

	public static void simulateRunsDefaultReportCV1 (MonteCarloModelCV model, int n,
			RandomStream stream, double[] meanPayoff, double[] varPayoff, Chrono timer) {
		timer.init();
		simulateRunsCV1 (model, n, stream, meanPayoff, varPayoff);
		// statWithCV.setConfidenceIntervalStudent();
		System.out.println(model.toString());
		// System.out.println(statWithCV.report(0.95, 4));
		System.out.printf("Average, no CV:  %9.5g%n", meanPayoff[0]);
		System.out.printf("Average WITH CV: %9.5g%n", meanPayoff[1]);
		System.out.printf("Variance, no CV:  %9.5g%n", varPayoff[0]);
		System.out.printf("Variance WITH CV: %9.5g%n", varPayoff[1]);
		System.out.println("Total CPU time:      " + timer.format() + "\n");
	}
	
}
