package umontreal.ssj.stat.density;

import java.util.ArrayList;


/**
 * This class provides methods to compute density estimators and to estimate their integrated variance (IV),
 * based on data from m independent replications, with n data points each.  
 * The data may have been obtained in any way, e.g., by Monte Carlo or by RQMC.
 * The density estimator is evaluated at #numEvalPoints points either equally spaced or specified explicitly
 * over a given interval [a,b], for each of the m replicates. 
 * By computing the empirical average and variance at each evaluation point, one can estimate the density and the IV.
 */

/**
 * @author Pierre L'Ecuyer
 * 
 */
public class DensityIVEstimation {


	/**
	 * Takes data from previous simulations (m independent replications, with n data points each),
	 * a density estimator de, and an array of evaluation points for this density over the interval [a,b].
	 * Computes and returns the density estimator and the empirical integrated variance (IV) 
	 * for this density estimator, obtained by estimating the variance at the 
	 * numEvalPoints evaluation points provided and taking the average multiplied by (b-a).
	 */
	public static double computeDensityVariance (int m, double[][] data, 
			DensityEstimator de, double a, double b, int numEvalPoints, double[] evalPoints) {
		double x, y;
		// TO DO:
		// If the density estimator is a histogram, here we may reset numEvalPoints to 
		// the number of bins of the histogram.
		// int numEvalPoints = evalPoints.length();
		double estimDens[] = new double[numEvalPoints];   // Value of the density at those points
		double meanDens[] = new double[numEvalPoints];    // Average value over the rep replicates
		double varDens[] = new double[numEvalPoints];     // Variance at each evaluation point
		// Arrays.fill(meanDens, 0.0);
		// Arrays.fill(varDens, 0.0);
		for (int rep = 0; rep < m; rep++) {
			// Estimate the density for this rep and evaluate it at the evaluation points
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
		double sumVar = 0.0;
		for (int i = 0; i < numEvalPoints; ++i)
			sumVar += varDens[i];
		return sumVar * (b - a) / (numEvalPoints * (m - 1));   // Empirical integrated variance.
	}
	

	/**
	 * Similar to the previous method, except that the numEvalPoints evaluation points are equidistant.
	 * We partition [a,b] in numEvalPoints equal intervals and put one point in the middle of each interval.
	 * With delta = (b-a) / numEvalPoints, we put the first point at a + delta/2, and then the successive 
	 * points are spaced by delta. 
	 */
	public static double computeDensityVariance (int m, double[][] data, 
			DensityEstimator de, double a, double b, int numEvalPoints) {
		double evalPoints[] = new double[numEvalPoints];  // Points at which the density will be evaluated
		double delta = (b-a) / (numEvalPoints);
		for (int j = 0; j < numEvalPoints; j++)
            evalPoints[j] = a + delta * (0.5 + j);
		return computeDensityVariance (m, data, de, a, b, numEvalPoints, evalPoints);
	}
	
		
	/**
	 * Similar to computeDensityVariance but does it for a list of density estimators,
	 * and returns the results (integrated variance for each DE) in array #integVariance.
	 * NOTE: In case we would like to provide more results, we may want to return
	 * them in a list like listDE.  
	 * 
	 * We may also want a method somewhere that plots the densities in a listDE ...  
	 */
	public static void computeDensityVarianceListDE (int m, double[][] data, 
			ArrayList<DensityEstimator> listDE, double a, double b, int numEvalPoints, 
			double[] evalPoints, double[] integVariance) {
        integVariance = new double[listDE.size()];
        int deNumber = 0;
		for (DensityEstimator de : listDE) {
            integVariance[deNumber] = computeDensityVariance (m, data, de, a, b, numEvalPoints, evalPoints);
            deNumber++;
		}
	}
	
	/**
	 * Similar to computeDensityVariance but does it for a list of density estimators,
	 * and returns the results (integrated variance for each DE) in array integVariance.
	 * NOTE: In case we would like to provide more results, we may want to return
	 * them in a list like listDE, but for the results (?).
	 */
	public static void computeDensityVarianceListDE (int m, double[][] data, 
			ArrayList<DensityEstimator> listDE, double a, double b, int numEvalPoints, 
			double[] integVariance) {
		double evalPoints[] = new double[numEvalPoints];  // Points at which the density will be evaluated
		double delta = (b-a) / (numEvalPoints);
		for (int j = 0; j < numEvalPoints; j++)
            evalPoints[j] = a + delta * (0.5 + j);
		computeDensityVarianceListDE (m, data, listDE, a, b, numEvalPoints, evalPoints, integVariance); 
	}
		
}
