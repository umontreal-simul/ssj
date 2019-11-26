package umontreal.ssj.mcqmctools.florian.examples;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import umontreal.ssj.hups.BakerTransformedPointSet;
import umontreal.ssj.hups.CachedPointSet;
import umontreal.ssj.hups.FaureSequence;
import umontreal.ssj.hups.IndependentPointsCached;
import umontreal.ssj.hups.LMScrambleShift;
import umontreal.ssj.hups.NestedUniformScrambling;
import umontreal.ssj.hups.PointSet;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.RQMCPointSet;
import umontreal.ssj.hups.RandomShift;
import umontreal.ssj.hups.Rank1Lattice;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.hups.SortedAndCutPointSet;
import umontreal.ssj.hups.StratifiedUnitCube;
import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.mcqmctools.RQMCExperimentSeries;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.PgfDataTable;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.Num;

/**
 * Class containing a main() function to run an experiment for RQMC mean estimation.
 * @author florian
 *
 */

public class TestRQMCExperimentDouble {

	public static void main(String[] args) throws IOException {
		MonteCarloModelDouble model;
		String modelDescr;
		RandomStream noise = new MRG32k3a();
		int dim;
		String outdir = "";
		FileWriter fw;
		File file;
		StringBuffer sb = new StringBuffer("");


//		double a = 0.5;
//		dim = 5;
//		model = new GFunction(a, dim);
//		sb.append(model.toString());
//		modelDescr = "GFunction";
//		System.out.println(sb.toString());
		
		double a = 2;
		dim = 5;
		double u = 0.5;
		model = new GenzGaussianPeak(a, u,dim);
		sb.append(model.toString());
		modelDescr = "GenzGaussianPeak";
		System.out.println(sb.toString());


		// Define the RQMC point sets to be used in experiments.
		int basis = 2; // Basis for the loglog plots.
		int numSkipReg = 0; // Number of sets skipped for the regression.
		int mink = 13; // first log(N) considered
		int i;
		int m = 500; // Number of RQMC randomizations.
//		int[] N = { 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152 }; // 13
		int[] N = { 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576 };
		int numSets = N.length; // Number of sets in the series.
		
		//Lattice Parameters:
		// 0.6^k
				int[][] aa = {
						{ 1, 3455, 1967, 1029, 2117, 3871, 533, 2411, 1277, 2435, 1723, 3803, 1469, 569, 1035, 3977, 721, 797,
								297, 1659 }, // 13

						{ 1, 6915, 3959, 7743, 3087, 5281, 6757, 3369, 7107, 6405, 7753, 1641, 3613, 1819, 5827, 2087, 4417,
								6909, 5623, 4739 }, // 14

						{ 1, 12031, 14297, 677, 6719, 15787, 10149, 7665, 1017, 2251, 12105, 2149, 16273, 14137, 8179, 6461,
								15051, 6593, 12763, 8497 }, // 15

						{ 1, 19463, 8279, 14631, 12629, 26571, 30383, 1337, 6431, 3901, 12399, 20871, 5175, 3111, 26857, 15111,
								22307, 30815, 25901, 27415 }, // 16

						{ 1, 38401, 59817, 33763, 32385, 2887, 45473, 48221, 3193, 63355, 40783, 37741, 54515, 11741, 10889,
								17759, 6115, 18687, 19665, 26557 }, // 17

						{ 1, 100135, 28235, 46895, 82781, 36145, 36833, 130557, 73161, 2259, 3769, 2379, 80685, 127279, 45979,
								66891, 8969, 56169, 92713, 67743 }, // 18

						{ 1, 154805, 242105, 171449, 27859, 76855, 183825, 38785, 178577, 18925, 260553, 130473, 258343, 79593,
								96263, 36291, 2035, 198019, 15473, 148703 }, // 19

						{ 1, 387275, 314993, 50301, 174023, 354905, 303021, 486111, 286797, 463237, 211171, 216757, 29831,
								155061, 315509, 193933, 129563, 276501, 395079, 139111 } // 20
				};

				// 0.6^k
				int[] aMult = { 1, 103259, 511609, 482163, 299529, 491333, 30987, 286121, 388189, 39885, 413851, 523765, 501705,
						93009, 44163, 325229, 345483, 168873, 376109, 146111 };


		// Create a list of series of RQMC point sets.
		ArrayList<RQMCPointSet[]> listRQMC = new ArrayList<RQMCPointSet[]>();
		PointSet p;
		PointSetRandomization rand;
		RQMCPointSet[] rqmcPts;
		String[] merits;
		String path;

//		// Independent points (Monte Carlo)
//		rqmcPts = new RQMCPointSet[numSets];
//		for (i = 0; i < numSets; ++i) {
//			p = new IndependentPointsCached(N[i], dim);
//			rand = new RandomShift(noise);
//			rqmcPts[i] = new RQMCPointSet(p, rand);
//		}
//		rqmcPts[0].setLabel("Independent points");
//		listRQMC.add(rqmcPts);

		// Stratification
//		 rqmcPts = new RQMCPointSet[numSets];
//		 int k;
//		 for (i = 0; i < numSets; ++i) {
//		 k = (int) Math.round(Math.pow(Num.TWOEXP[i + mink], 1.0 / (double) (dim)));
//		 p = new StratifiedUnitCube(k, dim);
//		
//		 rand = new RandomShift(noise);
//		 rqmcPts[i] = new RQMCPointSet(p, rand);
//		 }
//		 rqmcPts[0].setLabel("Stratification");
//		 listRQMC.add(rqmcPts);

		// lattice+shift
		rqmcPts = new RQMCPointSet[numSets];
		for (i = 0; i < numSets; ++i) {

//			p = new Rank1Lattice(N[i], aa[i], dim);
			p = new Rank1Lattice(N[i], aMult, dim);

			rand = new RandomShift(noise);
			rqmcPts[i] = new RQMCPointSet(p, rand);
		}
		rqmcPts[0].setLabel("Lattice+Shift");
		listRQMC.add(rqmcPts);

		// lattice+baker
		rqmcPts = new RQMCPointSet[numSets];
		for (i = 0; i < numSets; ++i) {

//		 p =  new BakerTransformedPointSet(new Rank1Lattice(N[i],aa[i],dim));
			p = new BakerTransformedPointSet(new Rank1Lattice(N[i], aMult, dim));

			rand = new RandomShift(noise);
			rqmcPts[i] = new RQMCPointSet(p, rand);
		}
		rqmcPts[0].setLabel("Lattice+Baker");
		listRQMC.add(rqmcPts);

//		// Sobol + LMS
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
//			CachedPointSet cp = new CachedPointSet(new SobolSequence(N[i], dim));
//
//			cp.setRandomizeParent(false);
//			p = cp;
//
//			rand = new NestedUniformScrambling(noise, i + mink + 1);
//			rqmcPts[i] = new RQMCPointSet(p, rand);
//		}
//		rqmcPts[0].setLabel("Sobol+NUS");
//		listRQMC.add(rqmcPts);

		boolean makePgfTable = true;
		boolean printReport = true;
		boolean details = true;

		// Perform an experiment with the list of series of RQMC point sets.
		// This list contains two series: lattice and Sobol.
		ArrayList<PgfDataTable> listCurves = new ArrayList<PgfDataTable>();
		RQMCExperimentSeries experSeries = new RQMCExperimentSeries(listRQMC.get(0), basis);
		experSeries.setExecutionDisplay(details);
		file = new File(outdir + "reportRQMC.txt");
//		file.getParentFile().mkdirs();
		fw = new FileWriter(file);

		fw.write(modelDescr + experSeries.testVarianceRateManyPointTypes(model, listRQMC, m, numSkipReg, makePgfTable,
				printReport, details, listCurves));
		fw.close();

		// Produces LaTeX code to draw these curves with pgfplot.
		sb = new StringBuffer("");
		sb.append(PgfDataTable.pgfplotFileHeader());
		sb.append(PgfDataTable.drawPgfPlotManyCurves(modelDescr + ": Mean values", "axis", 3, 1, listCurves, basis, "",
				" "));
		sb.append(PgfDataTable.pgfplotEndDocument());
		file = new File(outdir + "plotMean.tex");
		fw = new FileWriter(file);
		fw.write(sb.toString());
		fw.close();

		sb = new StringBuffer("");
		sb.append(PgfDataTable.pgfplotFileHeader());
		sb.append(PgfDataTable.drawPgfPlotManyCurves(modelDescr + ": Variance", "axis", 3, 4, listCurves, basis, "",
				" "));
		sb.append(PgfDataTable.pgfplotEndDocument());
		file = new File(outdir + "plotVariance.tex");
		fw = new FileWriter(file);
		fw.write(sb.toString());
		fw.close();

	}

}
