package umontreal.ssj.mcqmctools;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import umontreal.ssj.charts.XYLineChart;
import umontreal.ssj.functionfit.LeastSquares;
import umontreal.ssj.gof.GofStat;
import umontreal.ssj.hups.*;
import umontreal.ssj.hups.CachedPointSet;
import umontreal.ssj.hups.IndependentPointsCached;
import umontreal.ssj.hups.KorobovLattice;
import umontreal.ssj.hups.LMScrambleShift;
import umontreal.ssj.hups.NestedUniformScrambling;
import umontreal.ssj.hups.PointSet;
import umontreal.ssj.hups.PointSetIterator;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.RandomShift;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.hups.StratifiedUnitCube;
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
public class RQMCExperimentDensityU01 extends RQMCExperimentDensityKnown {
	
	static int numStats = 15;
	String[] statNames = { "var of average  ", 
			"square KS stat      ", "square L2star disc    ",
	        "histogram ISE, $m = 4 n^{1/3}$ ", "histogram ISE, $m = 16 n^{1/3}$ ", 
	        "histogram ISE, $m = 64 n^{1/3}$ ", "ASH, linear weights, $m = 256 n^{1/5}$",
	        "ASH, linear weights, $m = 512 n^{1/5}$", "ASH, linear weights, $m = 1024 n^{1/5}$",
	        "ASH, quadratic weights, $m = 256 n^{1/5}$", "ASH, quadratic weights, $m = 512 n^{1/5}$", 
	        "ASH, quadratic weights, $m = 1024 n^{1/5}$", "kernel dens. Botev, $h = n^{-1/5} /64$",
	        "kernel dens. Botev, $h = n^{-1/5} / 256$", "kernel dens. Botev, $h = n^{-1/5}/1024$" };
	static int[] numBins = new int[numStats];  // Should be inside simulReplicatesRQMC.
	
	/*	String[] statNames = { "var of average ", 
			"square KS stat      ", "square L2star disc    ",
			"histogram ISE, m = 16", "histogram ISE, m = 16 * 8",
	        "histogram ISE, m = 16 * 64", "ASH, linear weights, m = 16 * 64",
	        "ASH, linear weights, m = 128 * 64", "ASH, linear weights, m = 1024 * 64",
	        "ASH, quadratic weights, m = 16 * 64", "ASH, quadratic weights, m = 128 * 64",
	        "ASH, quadratic weights, m = 1024 * 64", "kernel dens. Botev, m = 256",
	        "kernel dens. Botev, m = 256 * 4", "kernel dens. Botev, m = 256 * 16" };
	static int[] numBins = { 1, 1, 1, 16, 128, 1024, 16 * 64, 128 * 64, 1024 * 64, 16 * 64, 128 * 64,
	        1024 * 64, 4096, 4096, 4096 };
	*/
	public double[] bandwidthKDE = new double[numStats];
	  //   {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.0/64.0, 1.0/256.0, 1.0/1024.0 };
	public double[] regSlope = new double[numStats]; // Regression slopes.

	
	// Performs n runs using stream and collects statistics in statValue.
	// Note that the stream can be an RQMC stream.
	public static void mapRunsTo01 (MonteCarloModelDensityU01 model, 
			TallyStore statValue, TallyStore statValue01) {
		statValue01.init();
		for (int i = 0; i < statValue.numberObs(); i++) 
			statValue01.add(model.cdf(statValue.getArray()[i]));
	}

    /**
	 * Computes and returns average of values of x-0.5.
	 */
	public static double meanCentered(double[] x) {
		double sum = 0.0;
		for (int i = 0; i < x.length; ++i)
			sum += x[i] - 0.5;
		return sum / x.length;
	}

	/**
	 * Computes and returns average of squares of values of x, minus 1/3.
	 */
	public static double meanSquareCentered(double[] x) {
		double sum2 = 0.0;
		for (int i = 0; i < x.length; ++i)
			sum2 += x[i] * x[i];
		return (sum2 / x.length) - 1.0 / 3.0;
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


	/**
	 * Simulate <tt>m/tt> replications and return a list of stat collectors that
	 * contain the <tt>m/tt> realizations of the KS and CvM test statistics, the
	 * centered mean, etc., for the states of the chains
	 * at the last step.
	 */
	public void simulReplicatesRQMC (MonteCarloModelDensityU01 model,
			PointSet p, PointSetRandomization rand, int m,
			double[] meanReps, double[] meanRepsMISE) {
		// Mean and variance of density (w.r.t. rep) in each bin or at each point.
		double meanDens[][] = new double[numStats][];
		double varDens[][] = new double[numStats][];
		double miseDens[][] = new double[numStats][];  // MSE of density in each bin.
		double exactDens[][] = new double[numStats][];  // Exact density in each bin.
		Tally statReps = new Tally();   // For average performance
		Tally statKS  = new Tally();   // For KS test
		Tally statCVM = new Tally();   // For Cramer von Mises

		int n = p.getNumPoints();
		// To store the n outputs X for each replication.
		TallyStore statValue = new TallyStore(n);
		TallyStore statValue01 = new TallyStore(n);
		PointSetIterator stream = p.iterator();
		double a = model.getMin();
		double b = model.getMax();

		TallyHistogram hist;
		ScaledHistogram histDensity;
		KernelDensityEstimator1d kde1d = new KernelDensityEstimator1d();

		// This is good for Normal2Sum.
        // Eventually, numBins should also depend on the RQMC method.
		double gamma = 0.20;
		numBins[0] = 1;
		numBins[1] = 1;
		numBins[2] = 1;
		numBins[3] = (int) (8 * Math.pow((double) n, 0.33333333));  // Histogram
		numBins[4] = (int) (16 * Math.pow((double) n, 0.33333333)); 
		numBins[5] = (int) (32 * Math.pow((double) n, 0.33333333));
		numBins[6] = (int) (256 * Math.pow((double) n, gamma));   // ASH
		numBins[7] = (int) (512 * Math.pow((double) n, gamma)); 
		numBins[8] = (int) (1024 * Math.pow((double) n, gamma));
		numBins[9] = (int) (256 * Math.pow((double) n, gamma));  // ASH quad
		numBins[10] = (int) (512 * Math.pow((double) n, gamma)); 
		numBins[11] = (int) (1024 * Math.pow((double) n, gamma));
		bandwidthKDE[12] = (1.0 / (64.0 * Math.pow((double) n, gamma)));  // Gaussian kernel
		bandwidthKDE[13] = (2.0 / (256.0 * Math.pow((double) n, gamma)));
		bandwidthKDE[14] = (4.0 / (1024.0 * Math.pow((double) n, gamma)));
		numBins[12] = numBins[13] = numBins[14] = 4096;  // Num eval. points for KDE.
		
		
		// For ASH.
		double[] weights2 = new double[numShiftsASH];   // Quadratic weights for ASH
		for (int j = 0; j < numShiftsASH; j++)
			weights2[j] = 1.0 - (double) (j * j) / ((double) (numShiftsASH * numShiftsASH));

		for (int j = 3; j < numStats; j++) { // For each statistic on density.
			Arrays.fill(meanDens[j] = new double[numBins[j]], 0.0);
			Arrays.fill(varDens[j] = new double[numBins[j]], 0.0);
			Arrays.fill(miseDens[j] = new double[numBins[j]], 0.0);
			exactDens[j] = new double[numBins[j]];
		    computeExactDensity (model, a, b, numBins[j], exactDens[j]);
		}
		for (int rep = 0; rep < m; rep++) {
			// p.randomize(rand);
			rand.randomize(p);
			stream.resetStartStream();
			simulateRuns(model, n, stream, statValue);
			statValue.quickSort();
			statReps.add(statValue.average());   // For the estimator of the mean.
			mapRunsTo01 (model, statValue, statValue01);
			double ks = computeKS(statValue01.getArray());
			statKS.add(ks * ks);   // Squared KS statistic.
            statCVM.add(computeSquareL2StarDim1(statValue01.getArray()));  // CVM stat.

			for (int j = 3; j < 6; j++) {
				hist = new TallyHistogram(a, b, numBins[j]);  //
				hist.fillFromArray(statValue.getArray());
				histDensity = new ScaledHistogram(hist, 1.0);
				updateDensityMoments(rep + 1, numBins[j], histDensity.getHeights(), meanDens[j],
				        varDens[j]);
				densityISE(numBins[j], histDensity.getHeights(), exactDens[j], miseDens[j]); 
			}

			for (int j = 6; j < 9; j++) {
				hist = new TallyHistogram(a, b, numBins[j]);  //
				hist.fillFromArray(statValue.getArray());
				histDensity = new ScaledHistogram(hist, 1.0);
				ScaledHistogram histAsh = histDensity.averageShiftedHistogramTrunc(numShiftsASH);
				histAsh.rescale(1.0);
				updateDensityMoments(rep + 1, numBins[j], histAsh.getHeights(), meanDens[j],
				        varDens[j]);
				densityISE(numBins[j], histAsh.getHeights(), exactDens[j], miseDens[j]); 
			}

			for (int j = 9; j < 12; j++) {
				hist = new TallyHistogram(a, b, numBins[j]);  //
				hist.fillFromArray(statValue.getArray());
				histDensity = new ScaledHistogram(hist, 1.0);
				ScaledHistogram histAsh = histDensity.averageShiftedHistogramTrunc(numShiftsASH,
				        weights2);
				histAsh.rescale(1.0);
				updateDensityMoments(rep + 1, numBins[j], histAsh.getHeights(), meanDens[j],
				        varDens[j]);
				densityISE(numBins[j], histAsh.getHeights(), exactDens[j], miseDens[j]); 
			}

			for (int j = 12; j < 15; j++) {
				kde1d.kde(statValue.getArray(), numBins[j], a, b, bandwidthKDE[j]);
				// kde1d.getDensity();
				updateDensityMoments(rep + 1, numBins[j], kde1d.getDensity(), meanDens[j],
				        varDens[j]);
			    densityISE(numBins[j], kde1d.getDensity(), exactDens[j], miseDens[j]); 
			}
		}
		meanReps[0] = statReps.variance();
		meanReps[1] = statKS.average();
		meanReps[2] = statCVM.average();
		for (int j = 3; j < numStats; j++) { // For each statistic on density estim.
			meanReps[j] = sum(varDens[j]) * (b - a) / (numBins[j] * (m - 1));
			meanRepsMISE[j] = sum(miseDens[j]) * (b - a) / (numBins[j] * (m - 1));
			// meanReps[j] = sum(varDens[j]) / (double)((m-1)*n); // - sumSquare(meanDens[j]);
			// System.out.println(" IV estimate = " + meanReps[j] + "\n");
		} // ****
	}

	/**
	 * Applies <tt>simulReplicatesArrayRQMC</tt> to run an experiment with a
	 * list of <tt>numSets/tt> point sets of various sizes, which are passed in
	 * the array <tt>pointSets/tt>. The descriptor string is used to describe
	 * this experiment in reports.
	 */

	public String simulReplicatesRQMCManySets(MonteCarloModelDensityU01 model, int numSets, int numSkipReg, PointSet[] pointSets, 
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
	    knownDensity = true;

	    StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for density estimation: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + m + "\n");
		sb.append(" Point sets: " + descPoints + "\n");
		sb.append(" Number of small point sets skipped for the regressions = " + numSkipReg + "\n\n");
		System.out.println(sb); // ****
		Chrono timer = new Chrono();
		System.out.println("    n        CPU time  ");
		for (int s = 0; s < numSets; s++) { // For each cardinality n
			n = pointSets[s].getNumPoints();
			log2n[s] = Num.log2(n);
			// System.out.println(" n = " + n + ", Lg n = " + log2n[s] + "\n"); // ****
			System.out.println("  " + n + "     " + timer.format());
			simulReplicatesRQMC(model, pointSets[s], rand, m, meanReps, meanRepsMISE);
			// sb.append(" Average centered state = " + meanReps[0] + "\n");
			// sb.append("\n Average centered square state = " + meanReps[1]);
			for (int j = 0; j < numStats; j++) {  // Each stat is now a discrepancy.
				log2Stats[j][s] = Num.log2(meanReps[j]);
	    	    log2StatsMISE[j][s] = Num.log2(meanRepsMISE[j]);
			}
		}
		for (int j = 0; j < numStats; j++) { // For each statistic (uniformity measure).
			sb.append(statNames[j] + "\n");
			sb.append("Variance (IV):\n");
			sb.append("    n       log2 mean square \n");
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				n = pointSets[s].getNumPoints();
				sb.append( n + "   " + PrintfFormat.f (7, 2, log2Stats[j][s]) + "\n");
			}
			// n = pointSets[numSets - 1].getNumPoints();
			// sb.append(statNames[j] + ", n = " + n + ", mean square = " + meanReps[j] + "\n");
			// regSlope[j] = LeastSquares.calcCoefficients(log2n, log2Stats[j], 1)[1];
			regSlope[j] = slope(log2n, log2Stats[j], numSkipReg, numSets-numSkipReg);
			sb.append("    Slope of log2 = " + PrintfFormat.f (8, 5, regSlope[j]) + "\n\n");
			if (j > 2) {
				sb.append("MISE:\n");
				for (int s = 0; s < numSets; s++) { // For each cardinality n
					n = pointSets[s].getNumPoints();
					sb.append(n + "   " + PrintfFormat.f (7, 2, log2StatsMISE[j][s]) + "\n");
				}
				regSlope[j] = slope(log2n, log2StatsMISE[j], numSkipReg, numSets-numSkipReg);
				sb.append("    Slope of log2 for MISE = " + PrintfFormat.f (8, 5, regSlope[j]) + "\n\n");
			}
		}
	    makePlotsMISE (numSets, m, (model.toString()).split(" ")[0], descPoints);

		// System.out.println(" Experiments done.\n");

		n = pointSets[numSets-1].getNumPoints();
		sb.append("***  Point sets: " + descPoints + " Summary for n = " + n + " \n \\hline \n");
		sb.append(" Number &  log Var  &          &  Var rate  &  name  \\\\ \\hline\n");
		for (int j = 0; j < 3; j++) { // For mean and cdf estimators.
			sb.append(j + "       & ");
			sb.append(PrintfFormat.f (7, 2, log2Stats[j][numSets-1]));
			sb.append(" &          & " + PrintfFormat.f (7, 3, regSlope[j]) + " & " + statNames[j] + " \\\\ \n");
		}
		sb.append("\\hline \n");
		sb.append(" Number &  log Var  & log MISE &  MISE rate &  name  \\\\ \n \\hline \n");
		for (int j = 3; j < numStats; j++) { // For each density estimator.
			sb.append(j + "       & ");
			sb.append(PrintfFormat.f (7, 2, log2Stats[j][numSets-1]) + " & ");
			sb.append(PrintfFormat.f (7, 2, log2StatsMISE[j][numSets-1]) + " & ");
			sb.append(PrintfFormat.f (7, 3, regSlope[j]) + " & " + statNames[j] + " \\\\ \n");
		}
		sb.append("\n");
		
	    double a = model.getMin();
		double b = model.getMax();
		fitPrintRegression (2, 3, 6, false,  numSkipReg, numSets-numSkipReg, b-a, 2.0, "Histogram", sb);
		fitPrintRegression (2, 6, 9, false,  numSkipReg, numSets-numSkipReg, b-a, 4.0, "ASH linear", sb);
		fitPrintRegression (2, 9, 12, false,  numSkipReg, numSets-numSkipReg, b-a, 4.0, "ASH quadratic", sb);
		fitPrintRegression (2, 12, 15, true,  numSkipReg, numSets-numSkipReg, b-a, 4.0, "KDE Gaussian", sb);

		sb.append("\n Total CPU Time = " + timer.format() + "\n");
		sb.append("-----------------------------------------------------\n");
		return sb.toString();
		// May want to return log2n, log2Stats[][], and regSlope[] as well.
	}
	
	// Method that produces a report on the experiment just made and make LaTeX
	// plots.
	public void makePlotsMISE (int numSets, int m, String descModel, String descPoints) {
		// makeGraph();
		try {
			// String title = descModel + "; " + descPoints;
			// double [][] bidon = new double[0][0];
			// XYLineChart chart = new XYLineChart (title, "lg n", "lg MISE", bidon);
			// chart.init (title, "lg n", "lg MISE");
		    XYLineChart chart = new XYLineChart ();
			for (int j = 3; j < numStats; j++)
				// chart.add(log2n, log2StatsMISE[j], statNames[j], " ");
			    chart.add(log2n, log2StatsMISE[j], PrintfFormat.d(j), " ");
			FileWriter file = new FileWriter(descModel + "_" + descPoints + "_MISE.tex");
			file.write(chart.toLatex(12, 8));
			file.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	

	// Method that does this for several types of point sets.
	public void TestRQMCManyPointTypes(MonteCarloModelDensityU01 model, int dim, int mink, 
			 int numSets, int numSkipReg, int[] N, int[] a, int m) {
		PointSet[] pointSets = new PointSet[numSets];
		RandomStream noise = new MRG32k3a();
		PointSetRandomization rand = new RandomShift(noise);
		// int dim = model.getDimension();
		int i, k;
	    knownDensity = true;

		// Independent points (Monte Carlo)
		for (i = 0; i < numSets; ++i)
			pointSets[i] = new IndependentPointsCached(N[i], dim);
		System.out.println(
			simulReplicatesRQMCManySets(model, numSets, numSkipReg, pointSets, rand, m, "IndepPoints"));

		// // Latin hypercube
		// for (i = 0; i < numSets; ++i)
		// pointSets[i] = new LatinHypercube(N[i], 2);
		// System.out.println(simulReplicatesArrayRQMCManySets(model, numSets, startSetReg, 
		// pointSets, rand, 1, numSteps, m, descModel,
		// "LatinHypercube"));
		//
		// Stratification
		for (i = 0; i < numSets; ++i) {
			k = (int) Math.round(Math.pow(Num.TWOEXP[i + mink], 1.0/dim));
			// n = k * k;
			pointSets[i] = new StratifiedUnitCube(k, dim);
			// pointSets[i] = new CachedPointSet2 (new StratifiedUnitCube (k,
			// 2), n, 1, 1);
		}
		System.out.println(simulReplicatesRQMCManySets(model, numSets, numSkipReg, pointSets, rand, m, 
				"Stratification"));
		//
		// // Stratification + antithetic
		// for (i = 0; i < numSets; ++i) {
		// k = (int) Math.round(Math.sqrt(Num.TWOEXP[i + mink] / 2.0));
		// // n = 2 * k * k;
		// pointSets[i] = new StratifiedUnitCubeAnti(k, 2);
		// }
		// System.out.println(simulReplicatesArrayRQMCManySets(model, numSets, startSetReg, 
		// pointSets, rand, 1, numSteps, m, descModel,
		// "StratificationAnti"));

		// Korobov lattice
		for (i = 0; i < numSets; ++i)
			pointSets[i] = new KorobovLattice(N[i], a[i], dim, 0);
		System.out.println(
				simulReplicatesRQMCManySets(model, numSets, numSkipReg, pointSets, rand, m, "Korobov2"));

		// Korobov lattice + baker
		// a[9] = 50673;
		for (i = 0; i < numSets; ++i)
			pointSets[i] = new BakerTransformedPointSet(new KorobovLattice(N[i], a[i], dim));
		System.out.println(
				simulReplicatesRQMCManySets(model, numSets, numSkipReg, pointSets, rand, m, "KorobovBaker2"));

		// Sobol
		for (i = 0; i < numSets; ++i)
			pointSets[i] = (new SobolSequence(i + mink, 31, dim)).toNetShiftCj();
		rand = new LMScrambleShift(new MRG32k3a());
		System.out.println(
				simulReplicatesRQMCManySets(model, numSets, numSkipReg, pointSets, rand, m, "Sobol2"));

		// Sobol + NUS
		for (i = 0; i < numSets; ++i) {
			DigitalNetBase2 sobolNet = (new SobolSequence(i + mink, 31, dim)).toNetShiftCj();
		    pointSets[i] = new CachedPointSet(sobolNet);
		}
		rand = new NestedUniformScrambling (new MRG32k3a(), 31);  // Scramble all 31 bits.
		System.out.println(
				simulReplicatesRQMCManySets(model, numSets, numSkipReg, pointSets, rand, m, "Sobol2 + NUS"));
	}

}
