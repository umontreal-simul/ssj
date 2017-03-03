package umontreal.ssj.mcqmctools;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.ArrayList;

import umontreal.ssj.charts.XYLineChart;
import umontreal.ssj.functionfit.LeastSquares;
import umontreal.ssj.gof.GofStat;
import umontreal.ssj.hups.*;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.*;
import umontreal.ssj.stat.density.DensityEstimator;
// import umontreal.ssj.stat.list.ListOfTallies;
import umontreal.ssj.util.Chrono;
import umontreal.ssj.util.Num;
import umontreal.ssj.util.PrintfFormat;

/**
 * Provides generic tools to perform RQMC experiments
 * with a simulation model that implements the MCModelDensity interface.
 */

/**
 * @author Pierre L'Ecuyer
 * 
 */
public class RQMCExperimentDensity extends RQMCExperiment {


	/**
	 * Takes data from previous simulation (m replicates, n points each)
	 * and a density estimator de.
	 * Computes and returns an estimate of the integrated variance (IV) 
	 * for this density estimator, obtained by estimating the variance 
	 * at numEvalPoints equidistant points over [a,b] and summing up.
	 */
	public static double computeDensityVariance (MonteCarloModelBounded model, int n, int m,
			double[][] data, DensityEstimator de, int numEvalPoints) {

		double x, y;
		// TO DO:
		// If the density estimator is a histogram, here we may reset numEvalPoints to 
		// the number of bins of the histogram.
		double evalPoints[] = new double[numEvalPoints];  // Points at which the density will be evaluated
		double estimDens[] = new double[numEvalPoints];   // Value of the density at those points
		double meanDens[] = new double[numEvalPoints];    // Average value over rep replicates
		double varDens[] = new double[numEvalPoints];     // Variance at each evaluation point
		// Arrays.fill(meanDens, 0.0);
		// Arrays.fill(varDens, 0.0);
		for (int rep = 0; rep < m; rep++) {
			// Estimate the density and evaluate it at eval points
			de.constructDensity(data[rep]);
			de.evalDensity(evalPoints, estimDens);
	        // Update the empirical mean and sum of squares of centered observations at each evaluation point.
			for (int j = 0; j < numEvalPoints; j++) {
				x = estimDens[j];
				y = x - meanDens[j];
				meanDens[j] += y / (double) (rep+1);
				varDens[j] += y * (x - meanDens[j]);
			}
		}
		double a = model.getMin();
		double b = model.getMax();
		double sumVar = 0.0;
		for (int i = 0; i < numEvalPoints; ++i)
			sumVar += varDens[i];
		return sumVar * (b - a) / (numEvalPoints * (m - 1));   // Empirical integrated variance.
	}
	

		
	/**
	 * Similar to computeDensityVariance but does it for a list of density estimators,
	 * and returns the results (integrated variance for each DE) in array integVariance.
	 * NOTE: In case we would like to provide more results, we may want to return
	 * them in a list like listDE, but for the results (?).
	 */
	public static void computeDensityVarianceListDE (MonteCarloModelBounded model, int n, int m,
			double[][] data, ArrayList<DensityEstimator> listDE, int numEvalPoints, 
            double[] integVariance) {
        integVariance = new double[listDE.size()];
        int deNumber = 0;
		for (DensityEstimator de : listDE) {
            integVariance[deNumber] = computeDensityVariance (model, n, m,
        			data, de, numEvalPoints);
            deNumber++;
		}
	}
		
}
