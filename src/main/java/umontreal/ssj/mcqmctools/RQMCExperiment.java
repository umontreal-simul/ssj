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
import umontreal.ssj.stat.list.lincv.ListOfTalliesWithCV;
// import umontreal.ssj.stat.list.ListOfTallies;
import umontreal.ssj.util.Chrono;
import umontreal.ssj.util.Num;
import umontreal.ssj.util.PrintfFormat;

/**
 * Provides generic tools to perform RQMC experiments
 * with a simulation model that implements the MonteCarloModelDouble interface.
 */

/**
 * @author Pierre L'Ecuyer
 * 
 */
public class RQMCExperiment extends MonteCarloExperiment {

	/*
	 * static int numStats = 12; String[] statNames = { "var of average ",
	 * "histogram ISE, fixed bin 16", "histogram ISE, 16 * 64 ", "histogram ISE, 32 x more bins ",
	 * "ASH, linear weights ", "ASH2: quadratic     ", "kernel dens. Botev  ", };
	 */
	
	public double[] log2n; // Log_2 of n.
	public double[] log2Stats; // log_2 of variance.
	public double regSlope; // Regression slopes.

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

	// Fits a regression using number values, starting from start.
	public double slope(double[] x, double[] y, int start, int number) {
		double[] x2 = new double[number], y2 = new double[number];
		for (int i = 0; i < number; ++i) {
			x2[i] = x[i+start];
			y2[i] = y[i+start];
		}
		return LeastSquares.calcCoefficients(x2, y2, 1)[1];
	}

 
	/**
	 * Simulate <tt>m/tt> replications and return the RQMC variance.
	 */
	public Tally simulReplicatesRQMC (MonteCarloModelDouble model, PointSet p, PointSetRandomization rand,
	        int m) {
		Tally statReps = new Tally();   // For average performance
		int n = p.getNumPoints();
		// To store the n outputs X for each replication.
		Tally statValue = new Tally();
		PointSetIterator stream = p.iterator();
		for (int rep = 0; rep < m; rep++) {
			rand.randomize(p);
			stream.resetStartStream();
			simulateRuns(model, n, stream, statValue);
			statReps.add(statValue.average());   // For the estimator of the mean.
		}
		return statReps;
	}

	/**
	 * Simulate <tt>m/tt> replications and return the RQMC variance.
	 */
	public static void simulReplicatesRQMC (MonteCarloModelDouble model, PointSet p, PointSetRandomization rand,
	        int m, Tally statReps) {
		// Tally statReps = new Tally();   // For average performance
		statReps.init();
		int n = p.getNumPoints();
		// To store the n outputs X for each replication.
		Tally statValue = new Tally();
		PointSetIterator stream = p.iterator();
		for (int rep = 0; rep < m; rep++) {
			rand.randomize(p);
			stream.resetStartStream();
			simulateRuns(model, n, stream, statValue);
			statReps.add(statValue.average());   // For the estimator of the mean.
		}
		// return statReps;
	}

	/**
	 * Makes <SPAN CLASS="MATH"><I>m</I></SPAN> independent randomizations of the point set
	 * <SPAN CLASS="MATH"><I>p</I></SPAN> using stream <TT>noise</TT>. For each of them, performs
	 * one simulation run for each point of <SPAN CLASS="MATH"><I>p</I></SPAN>, and adds the
	 * averages of the performance (value) and of the CVs over these points to 
	 * <TT>statWithCV</TT>.  This implements the replication method for RQMC with CVs,
	 * as in Section 5 of Hickernell, Lemieux and Owen (2005).
	 * CV estimates, variances, and covariances can be recovered from  <TT>statWithCV</TT>.
	 */
	public void simulReplicatesRQMCCV (MonteCarloModelCV model, PointSet p, 
			PointSetRandomization rand, int m, ListOfTalliesWithCV<Tally> statWithCV) {
		// Makes m independent randomizations of the point set p using stream
		// noise. For each of them, performs one simulation run for each point
		// of p, and adds the averages over these points to the statWithCV.

		statWithCV.init();
		int n = p.getNumPoints();
		int numCV = model.getNumberCV();
		PointSetIterator stream = p.iterator();
		for (int rep = 0; rep < m; rep++) {
			rand.randomize(p);
			stream.resetStartStream();
			double sumValues = 0.0;
			double[] sumValuesCV = new double[numCV]; 
			double[] curValuesCV = new double[numCV];  // Current value of the CV
			for (int i = 0; i < n; i++) {
				model.simulate(stream);
				sumValues += model.getPerformance();
				curValuesCV = model.getValuesCV();
				for (int k=0; k < numCV; k++) sumValuesCV[k] += curValuesCV[k];
				stream.resetNextSubstream();
			}
			for (int k=0; k < numCV; k++) 
				sumValuesCV[k] /= (double)n;
			statWithCV.add(sumValues/(double)n, sumValuesCV);
		}
	}
	

	//  RQMCPointSetsResults ???   Results for vector of point sets of different sizes.
	//  RQMCPointSetsArray  ???   Vector of point sets of different sizes.
	//  ListOfRQMCPointSetArrays 
	//  ListOfRQMCPointSets 

	/**
	 * Applies <tt>simulReplicatesArrayRQMC</tt> to run an experiment with a list of <tt>numSets/tt>
	 * point sets of various sizes, which are passed in the array <tt>pointSets/tt>. The descriptor
	 * string is used to describe this experiment in reports.
	 */

	public String simulReplicatesRQMCManySets(MonteCarloModelDouble model, int numSets, 
			int numSkipReg, PointSet[] pointSets,
	        PointSetRandomization rand, int m, String descPoints) {
		// ListOfTallies<TallyStore> statReps = ListOfTallies.createWithTallyStore(numStats);
		int n;
		double mean;  //  moyenne
		double var;  //  variance

		Tally statReps = new Tally();
		log2n = new double[numSets];
		log2Stats = new double[numSets]; // log_2 of each
		                                           // statistic for
		                                           // each n.
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for mean estimation: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + m + "\n");
		sb.append(" Point sets: " + descPoints + "\n\n");
		// System.out.println(sb); // ****
		Chrono timer = new Chrono();
		// System.out.println("    n        CPU time  ");
		sb.append("RQMC variance \n");
		sb.append("    n      mean       log2 var \n");
		System.out.println(sb);
		for (int s = 0; s < numSets; s++) { // For each cardinality n
			n = pointSets[s].getNumPoints();
			log2n[s] = Num.log2(n);
			// System.out.println(" n = " + n + ", Lg n = " + log2n[s] + "\n"); // ****
			// System.out.println("  " + n + "     " + timer.format());
			simulReplicatesRQMC (model, pointSets[s], rand, m, statReps);
			mean = statReps.average();
			var = statReps.variance();
			// sb.append(" Average state = " + meanReps[0] + "\n");
			// sb.append(" Average square state = " + meanReps[1] + "\n");
		    log2Stats[s] = Num.log2(var);
			n = pointSets[s].getNumPoints();
			System.out.println("  " + n + "     " + timer.format() + 
			              "   " + PrintfFormat.f(10, 5, mean) + 
					      "   " + PrintfFormat.f(7, 2, log2Stats[s]));
			// sb.append("  " + n + "   " + PrintfFormat.f(10, 5, mean) + 
			//		      "   " + PrintfFormat.f(7, 2, log2Stats[s]) + "\n");
		}
		double regSlope = slope(log2n, log2Stats, numSkipReg, numSets - numSkipReg);
		sb.delete(0, 1000);
		sb.append("    Slope of log2 = " + PrintfFormat.f(8, 5, regSlope) + "\n\n");
		// System.out.println(" Experiments done.\n");
	    // makePlotsIV (numSets, m, (model.toString()).split(" ")[0], descPoints);

		sb.append("\n Total CPU Time = " + timer.format() + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
		// May want to return log2n, log2Stats[][], and regSlope[] as well.
	}

	/**
	 * Applies <tt>simulReplicatesArrayRQMC</tt> to run an experiment with a list of <tt>numSets/tt>
	 * point sets of various sizes, which are passed in the array <tt>pointSets/tt>. The descriptor
	 * string is used to describe this experiment in reports.
	 */

	public String simulReplicatesRQMCManySetsCV (MonteCarloModelCV model, int numSets, 
			int numSkipReg, PointSet[] pointSets,
	        PointSetRandomization rand, int m, String descPoints) {
		int numCV = model.getNumberCV();
		ListOfTalliesWithCV<Tally> statWithCV = ListOfTalliesWithCV.createWithTally(1, numCV);
		statWithCV.setExpectedValue (0, 0.0);  // The CV is centered to 0.
		int n;
		double mean;  //  moyenne
		double var;  //  variance

		log2n = new double[numSets];
		log2Stats = new double[numSets]; // log_2 of each
		                                           // statistic for
		                                           // each n.
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for mean estimation: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + m + "\n");
		sb.append(" Point sets: " + descPoints + "\n\n");
		// System.out.println(sb); // ****
		Chrono timer = new Chrono();
		// System.out.println("    n        CPU time  ");
		sb.append("RQMC variance \n");
		sb.append("    n     cpu time     mean       log2 var \n");
		System.out.println(sb);
		for (int s = 0; s < numSets; s++) { // For each cardinality n
			n = pointSets[s].getNumPoints();
			log2n[s] = Num.log2(n);
			// System.out.println(" n = " + n + ", Lg n = " + log2n[s] + "\n"); // ****
			// System.out.println("  " + n + "     " + timer.format());
			simulReplicatesRQMCCV (model, pointSets[s], rand, m, statWithCV);
			statWithCV.estimateBeta();    // This is where the var. and covar. are computed!
			mean = statWithCV.averageWithCV(0);
			var = statWithCV.covarianceWithCV(0, 0);
			// sb.append(" Average state = " + meanReps[0] + "\n");
			// sb.append(" Average square state = " + meanReps[1] + "\n");
		    log2Stats[s] = Num.log2(var);
			n = pointSets[s].getNumPoints();
			System.out.println("  " + n + "     " + timer.format() + 
			              "   " + PrintfFormat.f(10, 5, mean) + 
					      "   " + PrintfFormat.f(7, 2, log2Stats[s]));
			// sb.append("  " + n + "   " + PrintfFormat.f(10, 5, mean) + 
			//		      "   " + PrintfFormat.f(7, 2, log2Stats[s]) + "\n");
		}
		double regSlope = slope(log2n, log2Stats, numSkipReg, numSets - numSkipReg);
		sb.delete(0, 1000);
		sb.append("    Slope of log2 = " + PrintfFormat.f(8, 5, regSlope) + "\n\n");
		// System.out.println(" Experiments done.\n");
	    // makePlotsIV (numSets, m, (model.toString()).split(" ")[0], descPoints);

		sb.append("\n Total CPU Time = " + timer.format() + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
		// May want to return log2n, log2Stats[][], and regSlope[] as well.
	}

	
	// Method that does this for several types of point sets.
	public void TestRQMCManyPointTypes(MonteCarloModelDouble model, int dim, int mink, int numSets, int numSkipReg, int[] N, int[] a,
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
		// System.out.println(
		 //       simulReplicatesRQMCManySets(model, numSets, numSkipReg, pointSets, rand, m, "Korobov"));

		// Korobov lattice + baker
		// a[9] = 50673;
		for (i = 0; i < numSets; ++i)
			pointSets[i] = new BakerTransformedPointSet(new KorobovLattice(N[i], a[i], dim));
		//System.out.println(
		 //       simulReplicatesRQMCManySets(model, numSets, numSkipReg, pointSets, rand, m, "Korobov+Baker"));

		// Sobol + nested uniform scramble
		for (i = 0; i < numSets; ++i) {
			DigitalNetBase2 sobolNet = (new SobolSequence(i + mink, 31, dim)).toNetShiftCj();
	        pointSets[i] = new CachedPointSet(sobolNet);
		}
		rand = new NestedUniformScrambling (new MRG32k3a(), 31);
		System.out.println(simulReplicatesRQMCManySets(model, numSets, numSkipReg, pointSets, rand,
		        m, "Sobol+NUS"));

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
		 // System.out.println(
		  //  simulReplicatesRQMCManySets(model, numSets, numSkipReg, pointSets, rand, m, "Sobol+LMS+baker"));

   }

	// Method that does this for several types of point sets, with CVs.
	public void TestRQMCManyPointTypesCV (MonteCarloModelCV model, int dim, int mink, int numSets, int numSkipReg, int[] N, int[] a,
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
		       simulReplicatesRQMCManySetsCV (model, numSets, numSkipReg, pointSets, rand, m, "IndepPoints"));

		// // Latin hypercube
		// for (i = 0; i < numSets; ++i)
		// pointSets[i] = new LatinHypercube(N[i], 2);
		// System.out.println(simulReplicatesArrayRQMCManySetsCV (model, numSets, startSetReg, 
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
		// System.out.println(simulReplicatesRQMCManySetsCV (model, numSets, startSetReg, pointSets, rand, m,
		// "Stratification"));
		//
		
		 // Stratification + antithetic
/*		 for (i = 0; i < numSets; ++i) {
			 k = (int) Math.round(Math.pow(Num.TWOEXP[i + mink], 1.0/dim));
		     pointSets[i] = new StratifiedUnitCubeAnti(k, dim);
		 }
		 System.out.println(simulReplicatesRQMCManySetsCV (model, numSets, numSkipReg, 
		    pointSets, rand, m, "StratificationAnti"));
*/
		
		// Korobov lattice
		for (i = 0; i < numSets; ++i)
			pointSets[i] = new KorobovLattice(N[i], a[i], dim, 0);
		// System.out.println(
		 //       simulReplicatesRQMCManySetsCV (model, numSets, numSkipReg, pointSets, rand, m, "Korobov"));

		// Korobov lattice + baker
		// a[9] = 50673;
		for (i = 0; i < numSets; ++i)
			pointSets[i] = new BakerTransformedPointSet(new KorobovLattice(N[i], a[i], dim));
		//System.out.println(
		 //       simulReplicatesRQMCManySetsCV (model, numSets, numSkipReg, pointSets, rand, m, "Korobov+Baker"));

		// Sobol + nested uniform scramble
		for (i = 0; i < numSets; ++i) {
			DigitalNetBase2 sobolNet = (new SobolSequence(i + mink, 31, dim)).toNetShiftCj();
	        pointSets[i] = new CachedPointSet(sobolNet);
		}
		rand = new NestedUniformScrambling (new MRG32k3a(), 31);
		System.out.println(simulReplicatesRQMCManySetsCV (model, numSets, numSkipReg, pointSets, rand,
		        m, "Sobol+NUS"));

		// Sobol + LMS + shift
		for (i = 0; i < numSets; ++i)
			pointSets[i] = (new SobolSequence(i + mink, 31, dim)).toNetShiftCj();
		rand = new LMScrambleShift(new MRG32k3a());
		System.out.println(
		        simulReplicatesRQMCManySetsCV (model, numSets, numSkipReg, pointSets, rand, m, "Sobol+LMS"));

		 // Sobol + Baker
		 for (i = 0; i < numSets; ++i)
		    pointSets[i] = new BakerTransformedPointSet((new SobolSequence(i + mink, 31,
		 dim)).toNetShiftCj());
		 // System.out.println(
		  //  simulReplicatesRQMCManySetsCV (model, numSets, numSkipReg, pointSets, rand, m, "Sobol+LMS+baker"));

   }
}
