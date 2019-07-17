package umontreal.ssj.markovchainrqmc.florian.chemicalReactionNetworks;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import umontreal.ssj.functionfit.LeastSquares;
import umontreal.ssj.hups.BakerTransformedPointSet;
import umontreal.ssj.hups.CachedPointSet;
import umontreal.ssj.hups.IndependentPointsCached;
import umontreal.ssj.hups.KorobovLattice;
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
import umontreal.ssj.markovchainrqmc.ArrayOfComparableChains;
import umontreal.ssj.markovchainrqmc.MarkovChainComparable;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.Num;
import umontreal.ssj.util.sort.BatchSort;
import umontreal.ssj.util.sort.HilbertCurveBatchSort;
import umontreal.ssj.util.sort.HilbertCurveSort;
import umontreal.ssj.util.sort.MultiDimSort;
import umontreal.ssj.util.sort.SplitSort;

/**
 * Class that provides a main() method to simulate a #ChemicalReactionNetwork
 * with array-RQMC when one uses one and the same sorting algorithm in every
 * step.
 * 
 * @author florian
 *
 */

public class TestChemicalReactionNetwork {

	public static void main(String[] args) throws IOException {

		ChemicalReactionNetwork model;

		/*
		 * ******************* REVERSIBLE ISO
		 **********************/
//		double epsInv = 1E2;
//		double alpha = 1E-4;
//		double[]c = {1.0,alpha};
//		double[] x0 = {epsInv};
//		double N0 = epsInv + epsInv/alpha;
//		double T = 1.6;
//		double tau = T/8.0;
//		
//		 model = new ReversibleIsomerization(c,x0,tau,T,N0);

		/*
		 * ******************* SCHLOEGL System
		 **********************/
//		double[] c = { 3E-7, 1E-4, 1E-3, 3.5 };
//		double[] x0 = { 250.0, 1E5};
//		double N0 = 2E5 + 1E5 + 250.0;
//		double T = 4.0;
////		double tau = T/20.0;
//		double tau =T/15.0;
//
//		model = new SchloeglSystem(c, x0, tau, T,N0);

		/*
		 * ******************* PKA
		 **********************/
		double[] c = { 8.696E-5, 0.02, 1.154E-4, 0.02, 0.016, 0.0017 };// Nano: 1E-9
		double[] x0 = { 33000.0, 33030.0, 1100.0, 1100.0, 1100.0, 1100.0 };
		double T = 0.00005;
		double tau = T / 15.0;

		model = new PKA(c, x0, tau, T);

		System.out.println(model.toString());

		String modelDescription = "pka-laststep";

		//whether the importance function should contain a constant term or not
		boolean bias = true;

		String dataFolder = "data/PKA/";
		model.init();

		ArrayOfComparableChains chain = new ArrayOfComparableChains(model);

		int[] N = { 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576 };
		int[] logN = { 13, 14, 15, 16, 17, 18, 19, 20 };
		int mink = 13;
		int numSets = N.length;
		int m = 100;

		// generating vector for lattice. order dep. weights 0.6^k
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

		int numSteps = (int) (T / tau);

		ArrayList<Integer> sortCoordPtsList = new ArrayList<Integer>();
		ArrayList<MultiDimSort> sortList = new ArrayList<MultiDimSort>();
		ArrayList<MultiDimSort> sortPts = new ArrayList<MultiDimSort>();

		// number of sample chains to read in and use for fitting importance function
		int rows = 262144; // rows = 1048576;

		// number of variables, i.e. species, plus 1. This equals the number of rows in
		// the data-files.
		int cols = 7; // PKA
//		int cols = 3; //Schloegl

		Scanner sc;
		double[][] vars = new double[rows][cols - 1];
		double[] response = new double[rows];
		double[] reg;

// location of data-file.
	

		sc = new Scanner(new BufferedReader(new FileReader(
				"/u/puchhamf/misc/jars/chemical/PKA/PKA/data/" + "MCData_Step_" + (numSteps - 1) + ".csv")));

		// Monomials that should be considered in the importance function. The integers
		// in each tuple indicate the coordinate indices of the variables contained in the monmial.
		int[][] reducedCols = { { 2 }, { 3 }, { 4 }, { 0, 1 }, { 1, 2 }, { 4, 5 }, { 0, 1, 1 }, { 1, 1, 2 },
				{ 4, 5, 5 } }; // PKA, linear
//		int[][] reducedCols = {{5}};   //PKA4


//		 //Schloegl, linear
//		int[][] reducedCols = { { 0 }, { 1 }, { 0, 0 }, { 0, 1 }, { 0, 0, 0 }, { 0, 0, 1 } };
//		


		//Fit the importance  function
		double[][] varsReduced = new double[rows][reducedCols.length];

		for (int i = 0; i < rows; i++) {
			String[] line = sc.nextLine().trim().split(",");
			response[i] = Double.parseDouble(line[cols - 1]);
			for (int j = 0; j < cols - 1; j++) {
				vars[i][j] = Double.parseDouble(line[j]);
			}
			int j = 0;
			for (int[] tuples : reducedCols) {
				varsReduced[i][j] = 1.0;
				for (int col : tuples) {
					varsReduced[i][j] *= vars[i][col];
				}
				++j;
			}
		}
		sc.close();

//		 reg = LeastSquares.calcCoefficients(vars, response);
		if (bias)
			reg = LeastSquares.calcCoefficients0(varsReduced, response);
		else
			reg = LeastSquares.calcCoefficients(varsReduced, response);

//		Add the desired sort
		sortList.add(new PKASort(reg,reducedCols,bias));
//		sortList.add(new SchloeglSystemSort(reg,reducedCols,bias));

		//Point will only be sorted by 1st coordinate.
		sortPts.add(new SplitSort(1));
		sortCoordPtsList.add(1);

		//standard sorts
		
		/* BATCH SORT */

		double[] batchExp = { 0.16, 0.16, 0.16, 0.16, 0.16, 0.2 };
//		sortList.add(new BatchSort<MarkovChainComparable>(batchExp));
//		sortPts.add(new BatchSort<MarkovChainComparable>(batchExp));
//		sortCoordPtsList.add(model.dimension());
//		modelDescription = "PKAReoredered-batch-sort";

		/* SPLIT SORT */
//		sortList.add(new SplitSort<MarkovChainComparable>(model.dimension()));
//		sortPts.add(new SplitSort<MarkovChainComparable>(model.dimension()));
//		sortCoordPtsList.add(model.dimension());
//		modelDescription = "reviso-split-sort";

		/* HILBERT BATCH SORT */
		sortList.add(new HilbertCurveBatchSort<MarkovChainComparable>(batchExp, 20));
		sortPts.add(new BatchSort(new double[] { 1.0 }));
		sortCoordPtsList.add(1);
		modelDescription = "pka-hilbert-batch";

		/* HILBERT CURVE SORT */
//		sortList.add(new HilbertCurveSort(model.dimension(), 20));
//		sortPts.add(new BatchSort(new double[] {1.0}));
//		sortCoordPtsList.add(1);
//		modelDescription = "pka-hilbert-curve";

		StringBuffer sb = new StringBuffer("");
		String str;
		String outFile = modelDescription + ".txt";

		RandomStream stream = new MRG32k3a();
		RQMCPointSet[] rqmcPts;
		PointSet[] pointSets = new PointSet[numSets];
		PointSetRandomization rand;
		RQMCPointSet prqmc;
		int i, s;

		int nMC = (int) 1E6; // nMC = 100;// n to estimate MC variance.
		Tally statMC = new Tally();
		statMC.init();
		// model.simulRunsWithSubstreams(nMC, model.numSteps, stream, statMC);
		model.simulRuns(nMC, model.numSteps, stream, statMC);
		double varMC = statMC.variance();
		str = "\n\n --------------------------\n";
		str += "MC average  = " + statMC.average() + "\n";
		str += "MC variance = " + varMC + "\n\n";
		sb.append(str);
		System.out.println(str);

		i = 0; // Sorts indexed by i
		for (MultiDimSort sort : sortList) {
			str = "****************************************************\n";
			str += "*\t" + sort.toString() + "\n";
			str += "****************************************************\n\n";
			sb.append(str);
			System.out.println(str);
			ArrayList<RQMCPointSet[]> listP = new ArrayList<RQMCPointSet[]>();

			// Independent points (Monte Carlo)
//			 rqmcPts = new RQMCPointSet[numSets];
//			 for (s = 0; s < numSets; ++s) {
//			 pointSets[s] = new IndependentPointsCached(N[s], model.K + model.N);
//			 rand = new RandomShift(stream);
//			 prqmc = new RQMCPointSet(pointSets[s], rand);
//			 rqmcPts[s] = prqmc;
//			 }
//			 rqmcPts[0].setLabel("Independent points");
//			 listP.add(rqmcPts);
//
//			// Stratification
//				rqmcPts = new RQMCPointSet[numSets];
//				int k;
//				for (s = 0; s < numSets; ++s) {
//					k = (int) Math.round(Math.pow(Num.TWOEXP[s + mink], 1.0 / (double) (sortCoordPtsList.get(i) + model.K)));
//					pointSets[s] = new StratifiedUnitCube(k, sortCoordPtsList.get(i) + model.K);
//					// Here the points must be sorted at each step, always.
//					// In the case of Hilbert map, the points should be 2d and sorted
//					// based on one coordinate,
//					// whereas the states are 2d and sorted by the Hilbert sort.
//					rand = new RandomShift(stream);
//					prqmc = new RQMCPointSet(pointSets[s], rand);
//					rqmcPts[s] = prqmc;
//				}
//				rqmcPts[0].setLabel("Stratification");
//				listP.add(rqmcPts);

//			 Lattice + Shift
			rqmcPts = new RQMCPointSet[numSets];
			for (s = 0; s < numSets; ++s) {

				pointSets[s] = new SortedAndCutPointSet(
						new Rank1Lattice(N[s], aa[s], sortCoordPtsList.get(i) + model.K), sortPts.get(i));

				rand = new RandomShift(stream);
				prqmc = new RQMCPointSet(pointSets[s], rand);
				rqmcPts[s] = prqmc;
			}
			rqmcPts[0].setLabel("lattice+shift");
			listP.add(rqmcPts);

			// Rank1Lattice +baker
			rqmcPts = new RQMCPointSet[numSets];
			for (s = 0; s < numSets; ++s) {

				// The points are sorted here, but only once.
				pointSets[s] = new SortedAndCutPointSet(
						new BakerTransformedPointSet(new Rank1Lattice(N[s], aa[s], sortCoordPtsList.get(i) + model.K)),
						sortPts.get(i));
//							 
				rand = new RandomShift(stream);
				prqmc = new RQMCPointSet(pointSets[s], rand);
				rqmcPts[s] = prqmc;
			}
			rqmcPts[0].setLabel("lattice+ baker ");
			listP.add(rqmcPts);

			// Korobov + Shift
//			rqmcPts = new RQMCPointSet[numSets];
//			 for (s = 0; s < numSets; ++s){
//				
//					pointSets[s] = new SortedAndCutPointSet (new KorobovLattice(N[s],korA[s],sortCoordPtsList.get(i)+model.K),sortPts.get(i));
//
//					 rand = new RandomShift(stream);
//				 prqmc = new RQMCPointSet(pointSets[s], rand);
//				 rqmcPts[s] = prqmc;
//			 }
//			 rqmcPts[0].setLabel("korobov+shift");
//			 listP.add(rqmcPts);

			// Korobov +baker
//			rqmcPts = new RQMCPointSet[numSets];
//			for (s = 0; s < numSets; ++s) {
//				
//					// The points are sorted here, but only once.
//					pointSets[s] = new SortedAndCutPointSet(new BakerTransformedPointSet(
//							new KorobovLattice(N[s],korA[s],sortCoordPtsList.get(i)+model.K)), sortPts.get(i));
////							 
//				rand = new RandomShift(stream);
//				prqmc = new RQMCPointSet(pointSets[s], rand);
//				rqmcPts[s] = prqmc;
//			}
//			rqmcPts[0].setLabel("korobov+baker ");
//			listP.add(rqmcPts);

//			 Sobol + LMS
			rqmcPts = new RQMCPointSet[numSets];
			for (s = 0; s < numSets; ++s) {

				pointSets[s] = new SortedAndCutPointSet(
						new SobolSequence(s + mink, 31, sortCoordPtsList.get(i) + model.K), sortPts.get(i));

				rand = new LMScrambleShift(stream);
				prqmc = new RQMCPointSet(pointSets[s], rand);
				rqmcPts[s] = prqmc;
			}
			rqmcPts[0].setLabel("Sobol+LMS");
			listP.add(rqmcPts);

			// Sobol+NUS
			rqmcPts = new RQMCPointSet[numSets];
			for (s = 0; s < numSets; ++s) {

				CachedPointSet p = new CachedPointSet(
						new SobolSequence(s + mink, 31, sortCoordPtsList.get(i) + model.K));
				p.setRandomizeParent(false);
				// The points are sorted here, but only once.
				pointSets[s] = new SortedAndCutPointSet(p, sortPts.get(i));

				rand = new NestedUniformScrambling(stream);
				prqmc = new RQMCPointSet(pointSets[s], rand);
				rqmcPts[s] = prqmc;
			}
			rqmcPts[0].setLabel("Sobol+NUS");
			listP.add(rqmcPts);

			// Sobol + LMS + Baker
//			rqmcPts = new RQMCPointSet[numSets];
//			for (s = 0; s < numSets; ++s) {
//				if (sortCoordPtsList.get(i) == 1)
//					pointSets[s] = new BakerTransformedPointSet(new SobolSequence(s + mink, 31, 1 + model.K));
//				else
//					pointSets[s] = new SortedAndCutPointSet(
//							new BakerTransformedPointSet(new SobolSequence(s + mink, 31, sortCoordPtsList.get(i) + model.K)),
//							sortList.get(i));
//
//				rand = new LMScrambleShift(stream);
//				prqmc = new RQMCPointSet(pointSets[s], rand);
//				rqmcPts[s] = prqmc;
//			}
//			rqmcPts[0].setLabel("Sobol+LMS+baker");
//			listP.add(rqmcPts);

			for (RQMCPointSet[] ptSeries : listP) {
				String label = ptSeries[0].getLabel();
				str = label;
				str += "\n-----------------------------\n";
				sb.append(str + "\n");
				System.out.println(str);
				// If Stratification, then we need to sort point set in every step
				int sortedCoords = label.startsWith("St") ? sortCoordPtsList.get(i) : 0;
				str = (chain.testVarianceRateFormat(ptSeries, sort, sortedCoords, model.numSteps, m, varMC,
						modelDescription + "-" + sort.toString() + "-" + label, label));
				System.out.println(str);
				sb.append(str + "\n");

			}
			i++;

		}
		FileWriter file = new FileWriter(outFile);
		file.write(sb.toString());
		file.close();

	}

}
