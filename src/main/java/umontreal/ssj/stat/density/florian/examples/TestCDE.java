package umontreal.ssj.stat.density.florian.examples;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import umontreal.ssj.functionfit.LeastSquares;
import umontreal.ssj.hups.BakerTransformedPointSet;
import umontreal.ssj.hups.CachedPointSet;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.DigitalNetBase2FromFile;
import umontreal.ssj.hups.IndependentPointsCached;
import umontreal.ssj.hups.LMScrambleShift;
import umontreal.ssj.hups.NestedUniformScrambling;
import umontreal.ssj.hups.PointSet;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.RQMCPointSet;
import umontreal.ssj.hups.RandomShift;
import umontreal.ssj.hups.Rank1Lattice;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.hups.StratifiedUnitCube;
import umontreal.ssj.latnetbuilder.DigitalNetSearch;
import umontreal.ssj.latnetbuilder.Search;
import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.mcqmctools.MonteCarloModelDoubleArray;
import umontreal.ssj.mcqmctools.RQMCExperiment;
import umontreal.ssj.mcqmctools.florian.examples.BucklingStrengthVars;
import umontreal.ssj.mcqmctools.florian.examples.LookbackOptionGBMVars;
import umontreal.ssj.mcqmctools.florian.examples.MultiNormalIndependent;
import umontreal.ssj.mcqmctools.florian.examples.San13Vars;
import umontreal.ssj.mcqmctools.florian.examples.San13VarsCDE;
import umontreal.ssj.mcqmctools.florian.examples.ShortColumnFunctionVars;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.PgfDataTable;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.density.ConditionalDensityEstimator;
import umontreal.ssj.stat.density.DensityEstimator;
import umontreal.ssj.stat.list.ListOfTallies;
import umontreal.ssj.stochprocess.BrownianMotion;
import umontreal.ssj.stochprocess.GeometricBrownianMotion;
import umontreal.ssj.util.Num;
import umontreal.ssj.util.PrintfFormat;

/**
 * This class provides all necessary tools as well as a main() function to run
 * an experiment for estimating densities with a conditional density estimator
 * (CDE) or a log-likelihood estimator (LR). Typically, this will require that a
 * #MonteCarloModelArray is provided together with a pertinent
 * #ConditionalDensityEstimator.
 * 
 * @author florian
 *
 */
public class TestCDE {
	private double a, b;
	private boolean displayExec;
	private boolean producePlots;

	private double[] logIV;
	private double[] logN;
	private double logOfBase = Math.log(2.0);
	private double baseOfLog = 2.0;

	private String[] tableFields = { "logN", "logIV" };

	public TestCDE(double a, double b) {
		this.a = a;
		this.b = b;
	}

	private void preprocess(RQMCPointSet[] rqmcPts) {
		logN = new double[rqmcPts.length];
		for (int i = 0; i < rqmcPts.length; i++)
			logN[i] = Math.log((double) rqmcPts[i].getNumPoints()) / logOfBase;
		logIV = new double[logN.length];
	}

	/**
	 * Constructs a formatted String summarizing the input for the experiment.
	 * 
	 * @param pointLabel     short string identifying the point sets.
	 * @param estimatorLabel short string describing the estimator.
	 * @param numEvalPoints  number of evaluation points.
	 * @param m              number of independent replications
	 * @return the input parameters as a formatted string.
	 */
	public String formatHead(String pointLabel, String estimatorLabel, String numEvalPoints, int m) {
		StringBuffer sb = new StringBuffer("");
		sb.append("Estimation of IV rate over [" + a + ", " + b + "]\n");
		sb.append("----------------------------------------------------------------\n\n");
		sb.append("Estimator: " + estimatorLabel + "\n");
		sb.append("Point set used: " + pointLabel + "\n");
		sb.append("Number of repititions: m = " + m + "\n");
		sb.append("Evaluation points: " + numEvalPoints + "\n");
		sb.append("----------------------------------------------------------------\n\n");
		if (displayExec)
			System.out.print(sb.toString());
		return sb.toString();
	}

	/**
	 * For the given \a model, this function estimates the IV for various RQMC point
	 * sets \a rqmcPts based on \a m independently constructed replications of the
	 * estimator and returns a formatted string to summarize the results in a table.
	 * 
	 * @param model   the underlying model
	 * @param rqmcPts the used RQMC points
	 * @param m       number of independent replications
	 * @param cde     the conditional density estimator.
	 * @param evalPts the evaluation points.
	 * @return
	 */
	public String estimateIVComputeTable(MonteCarloModelDoubleArray model, RQMCPointSet[] rqmcPts, int m,
			ConditionalDensityEstimator cde, double[] evalPts) {
		StringBuffer sb = new StringBuffer("");
		String str;

		str = "log(n)\t  log(IV)\n\n";

		sb.append(str);
		if (displayExec)
			System.out.print(str);

		double[][][] data;
		double[][] density;
		double[] variance;

		ListOfTallies<Tally> statRepsList = ListOfTallies.createWithTally(model.getPerformanceDim());

		for (int i = 0; i < rqmcPts.length; i++) {// rqmc point sets indexed by i
			data = new double[m][rqmcPts[i].getNumPoints()][model.getPerformanceDim()];
			RQMCExperiment.simulReplicatesRQMC(model, rqmcPts[i], m, statRepsList, data);

			density = new double[m][evalPts.length];
			for (int rep = 0; rep < m; rep++) {
				density[rep] = cde.evalDensity(evalPts, data[rep]);
//				System.out.println("TEST:\t" + rep + "\t" + density[rep][0]);
			}
			variance = new double[evalPts.length];

			logIV[i] = Math.log(DensityEstimator.computeIV(density, a, b, variance)) / logOfBase;

			str = PrintfFormat.f(3, 1, logN[i]) + "\t " + PrintfFormat.f(8, 6, logIV[i]) + "\n";

			sb.append(str);
			if (displayExec)
				System.out.print(str);
		}

		str = "\n\n";
		sb.append(str);
		if (displayExec)
			System.out.print(str);

		return sb.toString();
	}

	/**
	 * Based on @f$m@f$ replications of the experiment for various values of @f$n@f$
	 * that are internally stored, this function estimates the convergence rate of
	 * the MISE via linear regression. The result is returned in a formatted string.
	 * 
	 * @return a formatted string.
	 */
	public String estimateIVSlopes() {
		double[] regCoeffs = new double[2];

		String str = "Regression data:\n";
		str += "********************************************\n\n";
		str += "IV:\n";

		regCoeffs = LeastSquares.calcCoefficients(logN, logIV);
		str += "Slope:\t" + regCoeffs[1] + "\n";
		str += "Const.:\t" + regCoeffs[0] + "\n\n";

		if (displayExec)
			System.out.print(str);
		return str;
	}

	/**
	 * Carries out the simulation for one type of RQMC points of different
	 * sizes @f$n@f$ provided in \a rqmcPts. This function then estimates the MISE
	 * for every given @f$n@f$, estimates the convergence rate of the MISE and then
	 * returns a summary of the experiment as a formatted String.
	 * 
	 * @param model      the underlying model
	 * @param rqmcPts    the RQMC point sets
	 * @param m          the number of independent replications of the estimator.
	 * @param cde        the conditional density estimator
	 * @param evalPoints the evaluation points.
	 * @return a formatted String containing the main result of the experiment.
	 */
	public String testIVRate(MonteCarloModelDoubleArray model, RQMCPointSet[] rqmcPts, int m,
			ConditionalDensityEstimator cde, double[] evalPoints) {
		StringBuffer sb = new StringBuffer("");
		sb.append(formatHead(rqmcPts[0].getLabel(), cde.toString(), Integer.toString(evalPoints.length), m));
		preprocess(rqmcPts);
		sb.append(estimateIVComputeTable(model, rqmcPts, m, cde, evalPoints));
		sb.append(estimateIVSlopes());

		return sb.toString();
	}

	/**
	 * Carries out the simulation for several types of RQMC point sets of different
	 * sizes @f$n@f$ provided in \a rqmcPtsList. For each type of point set, this
	 * function then estimates the MISE for every given @f$n@f$, estimates the
	 * convergence rate of the MISE, produces a plot, and then returns a summary of
	 * the experiment as a formatted String.
	 * 
	 * @param model       the underlying model
	 * @param rqmcPtsList a list of different types of RQMC point sets for
	 *                    various @f$n@f$.
	 * @param m           number of independent repetitions
	 * @param cde         the density estimator.
	 * @param evalPoints  the evaluation points
	 * @return a formatted String containing the main result of the experiment.
	 * @throws IOException
	 */
	public String testIVRate(MonteCarloModelDoubleArray model, ArrayList<RQMCPointSet[]> rqmcPtsList, int m,
			ConditionalDensityEstimator cde, double[] evalPoints) throws IOException {
		StringBuffer sb = new StringBuffer("");
		ArrayList<PgfDataTable> pgfTblList = new ArrayList<PgfDataTable>();
		for (RQMCPointSet[] rqmcPts : rqmcPtsList) {
			sb.append(testIVRate(model, rqmcPts, m, cde, evalPoints));
			if (producePlots)
				pgfTblList.add(genPgfDataTable(rqmcPts[0].getLabel(), rqmcPts[0].getLabel()));
		}
		if (producePlots)
			genPlots(cde.toString(), pgfTblList);
		return sb.toString();
	}

	/**
	 * Creates the PgfDataTable used for the plot.
	 * 
	 * @param tableName
	 * @param tableLabel
	 * @return
	 */
	public PgfDataTable genPgfDataTable(String tableName, String tableLabel) {
		int len = logN.length;
		double[][] pgfData = new double[len][tableFields.length];
		for (int i = 0; i < len; i++) {
			pgfData[i][0] = logN[i];
			pgfData[i][1] = logIV[i];
		}
		return new PgfDataTable(tableName, tableLabel, tableFields, pgfData);
	}

	/**
	 * Creates a plot of the IV vs. @f$n@f$ in a log-log scale.
	 * 
	 * @param cdeDescr
	 * @param pgfTblList
	 * @throws IOException
	 */
	public void genPlots(String cdeDescr, ArrayList<PgfDataTable> pgfTblList) throws IOException {
		FileWriter fw;
		String plotBody;

		plotBody = PgfDataTable.drawPgfPlotManyCurves("log(IV) vs log(n)", "axis", 0, 1, pgfTblList, (int) baseOfLog,
				"", " ");
		fw = new FileWriter(cdeDescr + "_IV.tex");
		fw.write(PgfDataTable.pgfplotFileHeader() + plotBody + PgfDataTable.pgfplotEndDocument());
		fw.close();
	}

	private static double[] genEvalPoints(int numPts, double a, double b, RandomStream stream) {
		double[] evalPts = new double[numPts];
		double invNumPts = 1.0 / ((double) numPts);
		for (int i = 0; i < numPts; i++)
			evalPts[i] = a + (b - a) * ((double) i + stream.nextDouble()) * invNumPts;
		return evalPts;
	}

	/**
	 * Returns flag, whether to print output during runtime or not
	 * 
	 * @return
	 */
	public boolean getDisplayExec() {
		return displayExec;
	}

	/**
	 * Sets flag, whether to print output during runtime or not
	 */
	public void setDisplayExec(boolean displayExec) {
		this.displayExec = displayExec;
	}

	/**
	 * Returns flag, whether to produce the IV vs. @f$n@f$ plot or not.
	 * 
	 * @return
	 */
	public boolean getProducePlots() {
		return producePlots;
	}

	/**
	 * Sets flag, whether to produce the IV vs. @f$n@f$ plot or not.
	 * 
	 * @param producePlots the producePlots to set
	 */
	public void setProducePlots(boolean producePlots) {
		this.producePlots = producePlots;
	}

	public double getLogOfBase() {
		return logOfBase;
	}

	public void setLogOfBase(double logOfBase) {
		this.logOfBase = logOfBase;
	}

	public double getBaseOfLog() {
		return baseOfLog;
	}

	public void setBaseOfLog(double baseOfLog) {
		this.baseOfLog = baseOfLog;
	}

	public static void main(String[] args) throws IOException {

		/*
		 * ************************ UTIL PARAMETERS
		 ****************************************/

		RandomStream noise = new MRG32k3a();
//		((MRG32k3a)noise).setSeed(new long[] {1606215560,1140697538,1523809004,1292913007,1858992010,665386629,1173670500});
		int mink = 14; // first log(N) considered
		int i;
		int m = 100; // m= 20;// Number of RQMC randomizations.
//		int[] N = { 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152 }; 
		int[] N = { 16384, 32768, 65536, 131072, 262144, 524288 };

		int numSets = N.length; // Number of sets in the series.
		int numEvalPts = 128; // numEvalPts = 8;// normal: 128

		// 0.6^k
		int[] aMult = { 1, 103259, 511609, 482163, 299529, 491333, 30987, 286121, 388189, 39885, 413851, 523765, 501705,
				93009, 44163, 325229, 345483, 168873, 376109, 146111 };

		// 0.05^k
		int[][] aa = {
				{ 1, 3455, 1899, 2921, 3663, 2823, 3977, 2761, 255, 845, 3029, 3831, 2089, 3691, 1771, 3907, 337, 3735,
						1373, 1795 }, // 13

				{ 1, 6915, 4877, 7479, 1203, 3941, 2159, 3225, 5219, 6307, 2643, 633, 7139, 869, 7239, 7019, 8151, 3853,
						8019, 5731 }, // 14

				{ 1, 12033, 3801, 5023, 10647, 14127, 12751, 7461, 11901, 1167, 14349, 1951, 2209, 7397, 2505, 5675,
						12195, 1801, 7707, 13443 }, // 15

				{ 1, 25015, 11675, 7425, 3289, 17821, 5649, 32161, 10285, 12031, 26337, 13403, 14547, 18661, 7993, 1299,
						15111, 12735, 13129, 12655 }, // 16

				{ 1, 38401, 48799, 17301, 59639, 20297, 26805, 53109, 4365, 14055, 5023, 48499, 37937, 5155, 44255,
						61671, 11409, 38529, 61887, 19183 }, // 17

				{ 1, 96407, 36479, 31333, 63411, 80945, 24597, 41083, 70179, 42983, 62013, 48035, 80011, 105415, 108151,
						68869, 104973, 20719, 72257, 59193 }, // 18

				{ 1, 154805, 243089, 211205, 258913, 18107, 174117, 67287, 3585, 155767, 31401, 154275, 35513, 36509,
						162377, 51021, 88413, 190981, 145989, 257551 }, // 19

				{ 1, 387275, 457903, 282967, 117983, 355873, 439959, 109733, 382437, 297385, 267803, 68841, 343399,
						171303, 420841, 136437, 423733, 355591, 415917, 406205 } // 20
		};

		/*
		 * ************************ MODEL
		 ****************************************/
		// LOOKBACK OPTION

		int dim = 12;
		MonteCarloModelDoubleArray model = new LookbackOptionGBMVars(dim);

//		double[] mus = {2.9E7,500.0,1000.0}; //Canti
//		double[] sigmas = {1.45E6,100.0,100.0}; //Canti
//		int dim = mus.length;

		// CANTI, NORMALS
//		int dim = 11;
//		double[] mus = new double[dim];
//		Arrays.fill(mus,0.0);
//		double[] sigmas = new double[dim];
////		Arrays.fill(sigmas,1.0);
//		sigmas[0] = 1.0;
//		for(int j = 1; j<dim; j++)
//			sigmas[j] = sigmas[j-1] / Math.sqrt(2.0);
//
//		MonteCarloModelDoubleArray model = new MultiNormalIndependent(mus,sigmas);

		// BUCKLING
//		double[] mus = { 0.992 * 24.0, 1.05 * 0.5, 1.3 * 34.0, 0.987 * 29.0E3,0.35, 5.25};
//		double[] covs = { 0.028, 0.044, 0.1235, 0.076,0.05,0.07 };
//		int dim = mus.length;
//		double[] sigmas = new double[dim];
//		for(int j = 0; j < dim; j++)
//			sigmas[j] = mus[j] * covs[j];
//		MonteCarloModelDoubleArray model = new BucklingStrengthVars(mus, sigmas);

		// SAN
//		int dim = 13;

//		MonteCarloModelDoubleArray model = new San13Vars("san13a.dat"); //LR

//		MonteCarloModelDoubleArray model = new San13VarsCDE("san13a.dat"); //CDE

		// SHORT COLUMN
//		double h = 15.0;
//		double bb = 5.0;
//		double muY = 5.0;
//		double muM = 2000.0;
//		double muP = 500.0;
//		double sigmaY = 0.5;
//		double[][] sigma = { { 400.0, 50.0 }, { 0.0, 86.60254037844386 } };
//
//		MonteCarloModelDoubleArray model = new ShortColumnFunctionVars(muY, muM, muP, sigmaY, sigma);
//		int dim = 3;

		/*
		 * ************************ DENSITY ESTIMATOR
		 ****************************************/

		double strike = 101.0;
		double s0 = 100.0;
		double sigma = 0.12136;
		double r = 0.1;

		double a = strike;
		double b = strike + 34.4; // Lookback; cuts 0.08 left and 0.05 right --> 87% of mass

		ConditionalDensityEstimator cde = new LRLookbackOptionGBM(dim, s0, strike, r, sigma);
		String descr = "LRLookback";

		// CANTILEVER
//		double a = 0.407;
//		double b = 1.515;	
//		double D0 = 2.2535;
//		a = (a+1) * D0;
//		b = (b+1) * D0;
//		
//		double L = 100.0;
//		double t = 2.0;
//		double w = 4.0;
////		double[] weights = {0.25, -1.0 , 0.75};
//		double[] weights = {-1.0,1.0, -1.0};
//
//		ConditionalDensityEstimator cde = new CDECantilever(L, t, w, mus[0], sigmas[0], mus[1], sigmas[1], mus[2],
//				sigmas[2], weights);
//		String descr = "cdeCantiY2";

//		double pp = 0.76171875;

//		ConditionalDensityEstimator cde = new LRCantilever(L, t, w, mus[0],  sigmas[0],  mus[1],  sigmas[1], mus[2],  sigmas[2],pp);
//		String descr = "lrCantiOpt";

		// NORMALS
//		double a = -2.0;
//		double b = 2.0;
//		
//		double[] weights = new double[dim];
////		Arrays.fill(weights,1.0/((double) dim)); //set the unused ones negative, then they will be omitted!
//		
//		Arrays.fill(weights,-1.0); //set the unused ones negative, then they will be omitted!
//		weights[0] = 1.0;
//		
//		ConditionalDensityEstimator cde = new CDENormalizedSumOfNormals(dim,weights,sigmas);
//		String descr = "cdeSumOfNormals" + dim;

		// BUCKLING
//		double a = 0.5169;
//		double b = 0.6511;
//
////		double pp = 0.146171875; //ll opt
//
////		ConditionalDensityEstimator cde = new LRBucklingStrength(mus[4],sigmas[4],mus[5],sigmas[5],pp);
////		String descr = "LRBucklingStrength";
////		
//
////		double pp = 0.00244140625; //CDE opt
//
//		ConditionalDensityEstimator cde = new CDEBucklingStrength(mus[4],sigmas[4],mus[5],sigmas[5],pp);
//		String descr = "CDEBucklingStrength";

		// SAN
//		double a =22.0;
//		double b = 106.24; //95% non-centralized
//
////		ConditionalDensityEstimator cde = new LRSan13();
////		String descr = "LRSan13";
//
//		ConditionalDensityEstimator cde = new CDESan13("san13a.dat");
//		String descr = "CDESan";

		// SHORT COLUMN
//		double a = -5.338; // ShortColumn
//		double b = -0.528; // 99% centered
//
//		ConditionalDensityEstimator cde = new CDEShortColumnFunction(bb,h,muY,sigmaY);
//		ConditionalDensityEstimator cde = new LRShortColumnFunction(bb, h, muY, sigmaY, 0.5);
//		String descr = "LRSHortColumn";
//		String descr = "CDEShortColumn";

		double[] evalPoints = genEvalPoints(numEvalPts, a, b, noise);

		/*
		 * ************************ POINT SETS
		 ****************************************/
		// Create a list of series of RQMC point sets.
		ArrayList<RQMCPointSet[]> listRQMC = new ArrayList<RQMCPointSet[]>();
		PointSet p;
		PointSetRandomization rand;
		RQMCPointSet[] rqmcPts;

		// Independent points (Monte Carlo)
		rqmcPts = new RQMCPointSet[numSets];
		for (i = 0; i < numSets; ++i) {
			p = new IndependentPointsCached(N[i], dim);
			rand = new RandomShift(noise);
			rqmcPts[i] = new RQMCPointSet(p, rand);
		}
		rqmcPts[0].setLabel("Independent points");
		listRQMC.add(rqmcPts);
//		
//		// Stratification
//		rqmcPts = new RQMCPointSet[numSets];
//		int k;
//		for (i = 0; i < numSets; ++i) {
//			k = (int) Math.round(Math.pow(Num.TWOEXP[i + mink], 1.0 / (double) (dim)));
//			p = new StratifiedUnitCube(k, dim);
//
//			rand = new RandomShift(noise);
//			rqmcPts[i] = new RQMCPointSet(p, rand);
//		}
//		rqmcPts[0].setLabel("Stratification");
//		listRQMC.add(rqmcPts);

//		 //lattice+shift
//		 rqmcPts = new RQMCPointSet[numSets];
//		 for (i = 0; i < numSets; ++i) {
//		
////		 p = new Rank1Lattice(N[i],aa[mink-13+i],dim);
//			 p = new Rank1Lattice(N[i],aMult,dim);
//
//		
//		 rand = new RandomShift(noise);
//		 rqmcPts[i] = new RQMCPointSet(p, rand);
//		 }
//		 rqmcPts[0].setLabel("Lattice+Shift");
//		 listRQMC.add(rqmcPts);
//		
//		// lattice+baker
//				 rqmcPts = new RQMCPointSet[numSets];
//				 for (i = 0; i < numSets; ++i) {
//				
//				 p =  new BakerTransformedPointSet(new Rank1Lattice(N[i],aa[mink-13+i],dim));
////					 p =  new BakerTransformedPointSet(new Rank1Lattice(N[i],aMult,dim));
//
//				
//				 rand = new RandomShift(noise);
//				 rqmcPts[i] = new RQMCPointSet(p, rand);
//				 }
//				 rqmcPts[0].setLabel("Lattice+Baker");
//				 listRQMC.add(rqmcPts);

//		 Sobol + LMS
//		rqmcPts = new RQMCPointSet[numSets];
//		for (i = 0; i < numSets; ++i) {
//
//			p = new SobolSequence(i + mink, 31, dim);
//
//			rand = new LMScrambleShift(noise);
//			rqmcPts[i] = new RQMCPointSet(p, rand);
//		}
//		rqmcPts[0].setLabel("Sobol+LMS");
//		listRQMC.add(rqmcPts);
//
//		// Sobol+NUS
//		rqmcPts = new RQMCPointSet[numSets];
//		for (i = 0; i < numSets; ++i) {
//			CachedPointSet cp = new CachedPointSet(new SobolSequence(i + mink, 31, dim));
//			cp.setRandomizeParent(false);
//			p = cp;
//
//			rand = new NestedUniformScrambling(noise);
//			rqmcPts[i] = new RQMCPointSet(p, rand);
//		}
//		rqmcPts[0].setLabel("Sobol+NUS");
//		listRQMC.add(rqmcPts);

		TestCDE test = new TestCDE(a, b);

		test.setDisplayExec(true);
		test.setProducePlots(true);

		FileWriter fw = new FileWriter(descr + ".txt");
		String str = test.testIVRate(model, listRQMC, m, cde, evalPoints);

		fw.write(str);
		fw.close();

		System.out.println(str);
	}

}
