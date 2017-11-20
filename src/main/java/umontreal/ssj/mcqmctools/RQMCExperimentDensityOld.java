package umontreal.ssj.mcqmctools;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import umontreal.ssj.charts.XYLineChart;
import umontreal.ssj.functionfit.LeastSquares;
import umontreal.ssj.gof.GofStat;
import umontreal.ssj.hups.*;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.*;
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
public class RQMCExperimentDensityOld extends MonteCarloExperimentDouble {

	/*
	 * static int numStats = 12; String[] statNames = { "var of average ",
	 * "histogram ISE, fixed bin 16", "histogram ISE, 16 * 64 ", "histogram ISE, 32 x more bins ",
	 * "ASH, linear weights ", "ASH2: quadratic     ", "kernel dens. Botev  ", };
	 */

	public boolean knownDensity = false;   // Can be changed to "true" if density is known.
	
	static int numStats = 13;
	static int numShiftsASH = 32;  // Value of r for ASH.
	// static int numEvalKDE = 1024 * 16;  // Num of evaluations points of the KDE density.
	String[] statNames = { "var of average ", "histogram ISE, m = 64", "histogram ISE, m = 256",
	        "histogram ISE, m = 1024", "ASH, linear weights, m = 512",
	        "ASH, linear weights, m = 4096", "ASH, linear weights, m = 32768",
	        "ASH, quadratic weights, m = 512", "ASH, quadratic weights, m = 4096",
	        "ASH, quadratic weights, m = 32768", "kernel dens. Botev, h = 1/256",
	        "kernel dens. Botev, h = 1/512", "kernel dens. Botev, h = 1/1024" };
	static int[] numBins = { 1, 64, 256, 1024, 16 * 32, 128 * 32, 1024 * 32, 16 * 32, 128 * 32,
	        1024 * 32, 4096, 4096, 4096 };
	public double[] bandwidthKDE = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.0/256.0, 1.0/512.0, 1.0/1024.0 };
	// public double[] bandwidthKDE = new double[numStats];
	
	public double[] log2n; // Log_2 of n.
	public double[][] log2Stats; // log_2 of each statistic for each n.
	public double[][] log2StatsMISE; // log_2 of each statistic for each n, for known density.
	public double[] regSlope = new double[numStats]; // Regression slopes.



	/**
	 * Computes and returns average of values of state.
	 */
	public void setKnownDensity(boolean known) {
		knownDensity = known;
	}

	/**
	 * Computes and returns average of values of state.
	 */
	public static double meanCentered(double[] x) {
		double sum = 0.0;
		for (int i = 0; i < x.length; ++i)
			sum += x[i] - 0.5;
		return sum / x.length;
	}

	/**
	 * Returns sum of entries of x.
	 */
	public static double sum(double[] x) {
		double sum = 0.0;
		for (int i = 0; i < x.length; ++i)
			sum += x[i];
		return sum;
	}

	/**
	 * Computes and returns sum of squares of values of state.
	 */
	public static double sumSquare(double[] x) {
		double sum2 = 0.0;
		for (int i = 0; i < x.length; ++i)
			sum2 += x[i] * x[i];
		return sum2;
	}

	public static double computeKS(double[] x) {
		double[] ks = GofStat.kolmogorovSmirnov(x);
		return ks[2];
	}

	public static double computeSquareL2StarDim1(double[] x) {
		int n = x.length;
		double discL2Star = 0;
		double[] w = new double[n];
		for (int j = 0; j < n; ++j) {
			w[j] = (j + 0.5) / n;
		}
		// Arrays.sort(x);
		for (int j = 0; j < n; ++j) {
			discL2Star += (w[j] - x[j]) * (w[j] - x[j]);
		}
		discL2Star /= n;
		discL2Star += 1.0 / (12.0 * n * n);
		return discL2Star;
	}

	/**
	 * Computes the exact density at n points equally spaced in [a,b]. 
	 */
	public static void computeExactDensity (MonteCarloModelDensityKnown model, double a, double b, int numBins,
			double[] exact) {
		double delta = (b-a) / numBins;
		double x = a + 0.5 * delta;
		for (int j = 0; j < numBins; j++) {
			// System.out.println(" Exact density, numBins = " + numBins + "\n");
			exact[j] = model.density(x);
			x += delta;
		}
	}
	
	/**
	 * Approximates the ISE by sampling the estimated density and the model exact density at n
	 * points equally spaced in [a,b]. 
	 */
	public static void densityISE (int numBins,
			double[] heights, double[] exact, double[] sumIse) {
	    double y;
		for (int j = 0; j < numBins; j++) {
			y = heights[j] - exact[j];   // model.density(x);
			sumIse[j] += y * y;
		}
	}

	/**
	 * Updates vectors that contain empirical mean and sum of squares of centered observations of
	 * density (or histogram height) on a grid of numBins points.
	 */
	public static void updateDensityMoments(int rep1, int numBins, double[] heights, double[] mean,
	        double sum2c[]) {
		double x, y;
		for (int j = 0; j < numBins; j++) {
			x = heights[j];
			y = x - mean[j];
			mean[j] += y / (double) rep1;
			sum2c[j] += y * (x - mean[j]);
		}
	}

	/**
	 * Updates vectors that contain empirical mean and sum of squares observations of density (or
	 * histogram height) on a grid of numBins points.
	 */
	public static void updateDensityMoments0(int rep, int numBins, double[] heights, double[] mean,
	        double sum2[]) {
		double x, y;
		for (int j = 0; j < numBins; j++) {
			x = heights[j];
			y = x - mean[j];
			mean[j] += y / (double) rep;
			sum2[j] += x * x;
			// sum2c[j] += y*(x - mean[j]);
		}
	}

	// Fits a regression using number values, starting from start.
	public double slope(double[] x, double[] y, int start, int number) {
		double[] x2 = new double[number], y2 = new double[number];
		for (int i = 0; i < number; ++i) {
			x2[i] = x[i+start];
			y2[i] = y[i+start];
		}
		return LeastSquares.calcCoefficients(x2, y2, 1)[1];
	}

    // Fits a regression and prints the coefficients, using log2n[j1..j2-1] 
	// and log2Stats[j1..j2-1].  
	// In each array, it uses number values, starting from start.
	// We add an interaction term 
	// if order = 3 and none if order = 2.  range = b-a.
	// If useBandwidth == true, we use bandwidthKDE[j] for the bandwidth h, otherwise
	// we use (b-a)/numBins[j].
	public void fitPrintRegression (int order, int j1, int j2, boolean useBandwidth, 
			  int start, int number,
			  double range, double alpha, String method, StringBuffer sb) {
		double[][] regDataX = new double[(j2-j1) * number][order];
		double[] regDataY = new double[(j2-j1) * number];
		double[] regDataYMISE = new double[(j2-j1) * number];
		double[] coef;
		double logh;
		for (int j = j1; j < j2; j++) { // For each m.
			if (useBandwidth)
				logh = Num.log2(bandwidthKDE[j]);  // For KDE.
			else
				logh = Num.log2(range / numBins[j]);  // Not true for KDE.
			for (int s = 0; s < number; s++) { // For each cardinality n
				regDataX[(j - j1) * number + s][0] = log2n[start+s];
				regDataX[(j - j1) * number + s][1] = logh;
				if (order > 2) regDataX[(j - j1) * number + s][2] = logh * log2n[start+s];
				regDataY[(j - j1) * number + s] = log2Stats[j][start+s];
				regDataYMISE[(j - j1) * number + s] = log2StatsMISE[j][start+s];
			}
		}
		coef = LeastSquares.calcCoefficients0(regDataX, regDataY);
		// System.out.println(" coef computed 1.\n");
		sb.append("  Regression coefficients for " + method + ".\n");
		sb.append("  C     = " + Math.exp(coef[0]) + "\n");
		sb.append("  beta  = " + -coef[1] + "\n");
		sb.append("  delta = " + -coef[2] + "\n");
		if (order > 2) sb.append("  inter = " + -coef[3] + "\n");
		sb.append("  gamma = " + (-coef[1])/(alpha - coef[2]) + "\n");	
		sb.append("  nu    = " + (-alpha * coef[1])/(alpha - coef[2]) + "\n\n");	

		if (knownDensity) {
			coef = LeastSquares.calcCoefficients0(regDataX, regDataYMISE);
			sb.append("  C for MISE     = " + Math.exp(coef[0]) + "\n");
			sb.append("  beta for MISE  = " + -coef[1] + "\n");
			sb.append("  delta for MISE = " + -coef[2] + "\n\n");
		}
	}

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


	/**
	 * Simulate <tt>m/tt> replications and return an array of stat collectors that contain the
	 * <tt>m/tt> realizations of the KS and CvM test statistics, the centered mean, etc., for the
	 * states of the chains at the last step.
	 */
	public void simulReplicatesRQMCFixedh (MonteCarloModelBounded model, PointSet p, PointSetRandomization rand,
	        int m, double[] meanReps, double[] meanRepsMISE) {
		// double stats[] = new double[numStats];
		// Mean and variance of density (w.r.t. rep) in each bin or at each point.
		double meanDens[][] = new double[numStats][];
		double varDens[][] = new double[numStats][];
		double miseDens[][] = new double[numStats][];  // MSE of density in each bin.
		double exactDens[][] = new double[numStats][];  // Exact density in each bin.

		Tally statReps = new Tally();   // For average performance
		int n = p.getNumPoints();
		// To store the n outputs X for each replication.
		TallyStore statValue = new TallyStore(n);
		PointSetIterator stream = p.iterator();
		double a = model.getMin();
		double b = model.getMax();

		TallyHistogram hist;
		ScaledHistogram histDensity;
		KernelDensityEstimator1d kde1d = new KernelDensityEstimator1d();

		// For ASH.
		double[] weights2 = new double[numShiftsASH];   // Quadratic weights for ASH
		for (int j = 0; j < numShiftsASH; j++)
			weights2[j] = 1.0 - (double) (j * j) / ((double) (numShiftsASH * numShiftsASH));

		for (int j = 1; j < numStats; j++) { // For each statistic
			Arrays.fill(meanDens[j] = new double[numBins[j]], 0.0);
			Arrays.fill(varDens[j] = new double[numBins[j]], 0.0);
			Arrays.fill(miseDens[j] = new double[numBins[j]], 0.0);
			if (knownDensity) 
				computeExactDensity ((MonteCarloModelDensityKnown)model, a, b, numBins[j], exactDens[j]);
		}

		for (int rep = 0; rep < m; rep++) {
			p.randomize(rand);
			stream.resetStartStream();
			simulateRuns(model, n, stream, statValue);
			statValue.quickSort();
			statReps.add(statValue.average());   // For the estimator of the mean.

			for (int j = 1; j < 4; j++) {
				hist = new TallyHistogram(a, b, numBins[j]);  //
				hist.fillFromArray(statValue.getArray());
				histDensity = new ScaledHistogram(hist, 1.0);
				updateDensityMoments(rep + 1, numBins[j], histDensity.getHeights(), meanDens[j],
				        varDens[j]);
				if (knownDensity) 
					densityISE(numBins[j], histDensity.getHeights(), exactDens[j], miseDens[j]); 
			}

			for (int j = 4; j < 7; j++) {
				hist = new TallyHistogram(a, b, numBins[j]);  //
				hist.fillFromArray(statValue.getArray());
				histDensity = new ScaledHistogram(hist, 1.0);
				ScaledHistogram histAsh = histDensity.averageShiftedHistogramTrunc(numShiftsASH);
				histAsh.rescale(1.0);
				updateDensityMoments(rep + 1, numBins[j], histAsh.getHeights(), meanDens[j],
				        varDens[j]);
				if (knownDensity) 
					densityISE(numBins[j], histAsh.getHeights(), exactDens[j], miseDens[j]); 
			}

			for (int j = 7; j < 10; j++) {
				hist = new TallyHistogram(a, b, numBins[j]);  //
				hist.fillFromArray(statValue.getArray());
				histDensity = new ScaledHistogram(hist, 1.0);
				ScaledHistogram histAsh = histDensity.averageShiftedHistogramTrunc(numShiftsASH,
				        weights2);
				histAsh.rescale(1.0);
				updateDensityMoments(rep + 1, numBins[j], histAsh.getHeights(), meanDens[j],
				        varDens[j]);
				if (knownDensity) 
					densityISE(numBins[j], histAsh.getHeights(), exactDens[j], miseDens[j]); 
			}

			for (int j = 10; j < 13; j++) {
				kde1d.kde(statValue.getArray(), numBins[j], a, b, bandwidthKDE[j]);
				kde1d.getDensity();
				updateDensityMoments(rep + 1, numBins[j], kde1d.getDensity(), meanDens[j],
				        varDens[j]);
				if (knownDensity) 
					densityISE(numBins[j], kde1d.getDensity(), exactDens[j], miseDens[j]); 
			}
		}
		meanReps[0] = statReps.variance();
		for (int j = 1; j < numStats; j++) { // For each statistic
			meanReps[j] = sum(varDens[j]) * (b - a) / (numBins[j] * (m - 1));
			if (knownDensity) 
			   meanRepsMISE[j] = sum(miseDens[j]) * (b - a) / (numBins[j] * (m - 1));
			// System.out.println(" IV estimate = " + meanReps[j] + "\n");
		} // ****
	}


	
	/**
	 * Simulate <tt>m/tt> replications and return a list of stat collectors that contain the
	 * <tt>m/tt> realizations of the KS and CvM test statistics, the centered mean, etc., for the
	 * states of the chains at the last step.
	 */
	public void simulReplicatesRQMC(MonteCarloModelBounded model, PointSet p, PointSetRandomization rand, int m,
	        double[] meanReps, double[] meanRepsMISE) {
		// double stats[] = new double[numStats];
		// Mean and variance of density (w.r.t. rep) in each bin or at each point.
		double meanDens[][] = new double[numStats][];
		double varDens[][] = new double[numStats][];
		int numBins[] = new int[numStats];

		Tally statReps = new Tally();   // For average performance
		int n = p.getNumPoints();
		// To store the n outputs X for each replication.
		TallyStore statValue = new TallyStore(n);
		PointSetIterator stream = p.iterator();
		double a = model.getMin();
		double b = model.getMax();

		// For histograms.
		double factorHist = 16;
		int bandwidthASH = 32;
		double factorASH = factorHist * bandwidthASH;
		double factorKDE = 32;
		TallyHistogram hist;
		ScaledHistogram histDensity;
		KernelDensityEstimator1d kde1d = new KernelDensityEstimator1d();
		numBins[1] = (int) factorHist;
		// numBins[2] = (int) (factorHist * Math.pow((double)n, 0.20)); // c1 * n^{1/5}
		numBins[2] = (int) factorHist * 64;
		numBins[3] = (int) (factorASH * Math.pow((double) n, 0.20));    //
		numBins[5] = numBins[4] = numBins[3];
		numBins[6] = 1024 * 16;
		double bandwidthKDE = Math.pow((double) n, -0.20) / factorKDE;
		double[] weights2 = new double[bandwidthASH];   // Quadratic weights for ASH
		for (int j = 0; j < bandwidthASH; j++)
			weights2[j] = 1.0 - (double) (j * j) / ((double) (bandwidthASH * bandwidthASH));

		// System.out.println(" numBins2 = " + numBins2 + " numBins3 = " + numBins3 + "\n");
		for (int j = 1; j < numStats; j++) { // For each statistic
			Arrays.fill(meanDens[j] = new double[numBins[j]], 0.0);
			Arrays.fill(varDens[j] = new double[numBins[j]], 0.0);
		}

		for (int rep = 0; rep < m; rep++) {
			rand.randomize(p);
			// p.randomize(rand);
			stream.resetStartStream();
			simulateRuns(model, n, stream, statValue);
			statValue.quickSort();
			statReps.add(statValue.average());   // For the estimator of the mean.

			hist = new TallyHistogram(a, b, numBins[1]);  //
			hist.fillFromArray(statValue.getArray());
			histDensity = new ScaledHistogram(hist, 1.0);
			updateDensityMoments(rep + 1, numBins[1], histDensity.getHeights(), meanDens[1],
			        varDens[1]);

			hist = new TallyHistogram(a, b, numBins[2]);  //
			hist.fillFromArray(statValue.getArray());
			histDensity = new ScaledHistogram(hist, 1.0);
			updateDensityMoments(rep + 1, numBins[2], histDensity.getHeights(), meanDens[2],
			        varDens[2]);

			hist = new TallyHistogram(a, b, numBins[3]);
			hist.fillFromArray(statValue.getArray());
			histDensity = new ScaledHistogram(hist, 1.0);
			ScaledHistogram histAsh = histDensity.averageShiftedHistogramTrunc(bandwidthASH);
			histAsh.rescale(1.0);
			updateDensityMoments(rep + 1, numBins[3], histDensity.getHeights(), meanDens[3],
			        varDens[3]);
			updateDensityMoments(rep + 1, numBins[4], histAsh.getHeights(), meanDens[4],
			        varDens[4]);

			histAsh = histDensity.averageShiftedHistogramTrunc(bandwidthASH, weights2);
			histAsh.rescale(1.0);
			updateDensityMoments(rep + 1, numBins[5], histAsh.getHeights(), meanDens[5],
			        varDens[5]);

			kde1d.kde(statValue.getArray(), numBins[6], a, b, bandwidthKDE);
			kde1d.getDensity();
			updateDensityMoments(rep + 1, numBins[6], kde1d.getDensity(), meanDens[6], varDens[6]);
		}
		meanReps[0] = statReps.variance();
		for (int j = 1; j < numStats; j++) { // For each statistic
			meanReps[j] = sum(varDens[j]) * (b - a) / (numBins[j] * (m - 1));
			// meanReps[j] = sum(varDens[j]) / (double)((m-1)*n); // - sumSquare(meanDens[j]);
			// System.out.println(" IV estimate = " + meanReps[j] + "\n");
		} // ****
	}

	/**
	 * Applies <tt>simulReplicatesArrayRQMC</tt> to run an experiment with a list of <tt>numSets/tt>
	 * point sets of various sizes, which are passed in the array <tt>pointSets/tt>. The descriptor
	 * string is used to describe this experiment in reports.
	 */

	//  RandPointSetsResults ???   Results for vector of point sets of different sizes.
	//  RandPointSetsVector  ???   Vector of point sets of different sizes.
	//  ListOfRandPointSetsVectors 
	
	public String simulReplicatesRQMCManySets(MonteCarloModelBounded model, int numSets, 
			int numSkipReg, PointSet[] pointSets,
	        PointSetRandomization rand, int m, String descPoints) {
		// ListOfTallies<TallyStore> statReps = ListOfTallies.createWithTallyStore(numStats);
		double[] meanReps = new double[numStats];     // mean for each statistic in statReps.
		double[] meanRepsMISE = new double[numStats]; // same, for MISE, for known density.
		int n;
		log2n = new double[numSets];
		log2Stats = new double[numStats][numSets]; // log_2 of each
		                                           // statistic for
		                                           // each n.
		log2StatsMISE = new double[numStats][numSets]; // log_2 of each for MISE.
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + m + "\n");
		sb.append(" Point sets: " + descPoints + "\n\n");
		System.out.println(sb); // ****
		Chrono timer = new Chrono();
		for (int s = 0; s < numSets; s++) { // For each cardinality n
			n = pointSets[s].getNumPoints();
			log2n[s] = Num.log2(n);
			// System.out.println(" n = " + n + ", Lg n = " + log2n[s] + "\n"); // ****
			System.out.println("  n = " + n);

			simulReplicatesRQMCFixedh(model, pointSets[s], rand, m, meanReps, meanRepsMISE);  // ******
			// sb.append(" Average state = " + meanReps[0] + "\n");
			// sb.append(" Average square state = " + meanReps[1] + "\n");
			for (int j = 0; j < numStats; j++) { // Each stat is now a discrepancy.
				log2Stats[j][s] = Num.log2(meanReps[j]);
			    if (knownDensity) 
			    	log2StatsMISE[j][s] = Num.log2(meanRepsMISE[j]);
			}
		}
		for (int j = 0; j < numStats; j++) { // For each statistic (uniformity
		                                     // measure).
			sb.append(statNames[j] + "\n");
			sb.append("    n       log2 mean square \n");
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				n = pointSets[s].getNumPoints();
				sb.append( n + "   " + PrintfFormat.f (7, 2, log2Stats[j][s]) + "\n");
			}
			// n = pointSets[numSets - 1].getNumPoints();
			// sb.append(statNames[j] + ", n = " + n + ", mean square = " + meanReps[j] + "\n");
			// regSlope[j] = LeastSquares.calcCoefficients(log2n, log2Stats[j],
			// 1)[1];
			regSlope[j] = slope(log2n, log2Stats[j], numSkipReg, numSets-numSkipReg);
			sb.append("    Slope of log2 = " + PrintfFormat.f (8, 5, regSlope[j]) + "\n\n");
			if (knownDensity & (j > 1)) {
   			   for (int s = 0; s < numSets; s++) { // For each cardinality n
				   n = pointSets[s].getNumPoints();
				   sb.append(n + "   " + PrintfFormat.f (7, 2, log2StatsMISE[j][s]) + "\n");
			   }
			   double regSlopeMISE = slope(log2n, log2StatsMISE[j], numSkipReg, numSets-numSkipReg);
			   sb.append("    Slope of log2 for MISE = " + PrintfFormat.f (8, 5, regSlopeMISE) + "\n\n");
			}
		}
		// System.out.println(" Experiments done.\n");
	    makePlotsIV (numSets, m, (model.toString()).split(" ")[0], descPoints);

		n = pointSets[numSets-1].getNumPoints();
		sb.append("***  Point sets: " + descPoints + " Summary for n = " + n + " \n \\hline \n");
		sb.append(" Number & log Var  & Var rate  &  name  \\\\ \\hline\n");
		for (int j = 0; j < numStats; j++) { // For mean and cdf estimators.
			sb.append(j + "       & ");
			sb.append(PrintfFormat.f (7, 2, log2Stats[j][numSets-1]));
			sb.append("  & " + PrintfFormat.f (7, 3, regSlope[j]) + "   & " + statNames[j] + " \\\\ \n");
		}
		sb.append("\n");

		double a = model.getMin();
		double b = model.getMax();

		fitPrintRegression (2, 1, 4, false,  numSkipReg, numSets-numSkipReg, b-a, 2.0, "Histogram", sb);
		fitPrintRegression (2, 4, 7, false,  numSkipReg, numSets-numSkipReg, b-a, 4.0, "ASH linear", sb);
		fitPrintRegression (2, 7, 10, false,  numSkipReg, numSets-numSkipReg, b-a, 4.0, "ASH quadratic", sb);
		fitPrintRegression (2, 10, 13, true,  numSkipReg, numSets-numSkipReg, b-a, 4.0, "KDE Gaussian", sb);

		sb.append("\n Total CPU Time = " + timer.format() + "\n");
		sb.append("-----------------------------------------------------\n");
		return sb.toString();
		// May want to return log2n, log2Stats[][], and regSlope[] as well.
	}

	// Method that produces LaTeX plots for the experiment just made.
	public void makePlotsIV (int numSteps, int m, String descModel, String descPoints) {
		// makeGraph();
		try {
			XYLineChart chart = new XYLineChart();
			for (int j = 1; j < numStats; j++)
				chart.add(log2n, log2Stats[j], PrintfFormat.d(j), " ");
			FileWriter file = new FileWriter(
			        descModel + "_" + descPoints + "_IV.tex");
			file.write(chart.toLatex(12, 8));
			file.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// Method that does this for several types of point sets.
	public void TestRQMCManyPointTypes(MonteCarloModelBounded model, int dim, int mink, int numSets, int numSkipReg, int[] N, int[] a,
	        int m) {
		PointSet[] pointSets = new PointSet[numSets];
		RandomStream noise = new MRG32k3a();
		PointSetRandomization rand = new RandomShift(noise);
		// int dim = model.getDimension();
		int i;
		int k;

		// Independent points (Monte Carlo)
		for (i = 0; i < numSets; ++i)
			pointSets[i] = new IndependentPointsCached(N[i], dim);
		System.out.println(
		        simulReplicatesRQMCManySets(model, numSets, numSkipReg, pointSets, rand, m, "IndepPoints"));

		// // Latin hypercube
		// for (i = 0; i < numSets; ++i)
		// pointSets[i] = new LatinHypercube(N[i], 2);
		// System.out.println(simulReplicatesArrayRQMCManySets(moel, numSets, startSetReg, 
		// pointSets, rand, 1, numSteps, m, descModel,
		// "LatinHypercube"));
		//
		// Stratification
		// for (i = 0; i < numSets; ++i) {
		// 	  k = (int) Math.round(Math.pow(Num.TWOEXP[i + mink], 1.0/dim));
		// n = k * k;
		// pointSets[i] = new StratifiedUnitCube(k, 2);
		// pointSets[i] = new CachedPointSet2 (new StratifiedUnitCube (k,
		// 2), n, 1, 1);
		// }
		// System.out.println(simulReplicatesRQMCManySets(model, numSets, startSetReg, pointSets, rand, m,
		// "Stratification"));
		//
		
		 // Stratification + antithetic
/*		 for (i = 0; i < numSets; ++i) {
			 k = (int) Math.round(Math.pow(Num.TWOEXP[i + mink], 1.0/dim));
		     pointSets[i] = new StratifiedUnitCubeAnti(k, dim);
		 }
		 System.out.println(simulReplicatesRQMCManySets(model, numSets, numSkipReg, 
		    pointSets, rand, m, "StratificationAnti"));
*/
		
		// Korobov lattice
		for (i = 0; i < numSets; ++i)
			pointSets[i] = new KorobovLattice(N[i], a[i], dim, 0);
		System.out.println(
		        simulReplicatesRQMCManySets(model, numSets, numSkipReg, pointSets, rand, m, "Korobov"));

		// Korobov lattice + baker
		// a[9] = 50673;
		for (i = 0; i < numSets; ++i)
			pointSets[i] = new BakerTransformedPointSet(new KorobovLattice(N[i], a[i], dim));
		System.out.println(
		        simulReplicatesRQMCManySets(model, numSets, numSkipReg, pointSets, rand, m, "Korobov+Baker"));

		// Sobol + LMS + shift
		for (i = 0; i < numSets; ++i)
			pointSets[i] = (new SobolSequence(i + mink, 31, dim)).toNetShiftCj();
		rand = new LMScrambleShift(new MRG32k3a());
		System.out.println(
		        simulReplicatesRQMCManySets(model, numSets, numSkipReg, pointSets, rand, m, "Sobol+LMS"));

		 // Sobol + Baker
		 for (i = 0; i < numSets; ++i)
		    pointSets[i] = new BakerTransformedPointSet((new SobolSequence(i + mink, 31,
		 dim)).toNetShiftCj());
		 System.out.println(
		    simulReplicatesRQMCManySets(model, numSets, numSkipReg, pointSets, rand, m, "Sobol+LMS+baker"));

   }

}
