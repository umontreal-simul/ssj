package latnetbuilder;

import umontreal.ssj.latnetbuilder.PolynomialLatticeSearch;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.stat.Tally;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.Set;

public class RunExperiments {
	private static PolynomialLatticeSearch initializePLRSearch(int s, int k, String pathToLNB, String pointSetType) {
		PolynomialLatticeSearch search = new PolynomialLatticeSearch(pointSetType);
		search.setPathToLatNetBuilder(pathToLNB);
		search.setDimension("" + s);
		search.setSizeParameter("2^" + k);
		return search;
	}

	private static void createExperimentFiles(Path experimentDirPath, String pathToLNB, String pathToJoeKuo,
			String experimentName, Set<Integer> experimentHashes) throws IOException {

		int s = 3;
		double[] c = {0.7, 0.2, 0.5};
		int r = -1;
		String w = "product:0:0.7,0.2,0.5";
		for (int k = 5; k <= 18; k += 1) {

			SobolTestFunc integrand = new SobolTestFunc(c, s, r);

			int n = 10000;
			Tally statValue = new Tally("MC method");
			integrand.simulateRuns(n, new MRG32k3a(), statValue);
			double varMC = statValue.variance();

			for (String randomization : new String[]{"LMS - Shift", "NUS - Interlace"})
			{
				PolynomialLatticeSearch search = initializePLRSearch(s, k, pathToLNB, "lattice");
				search.setFigureOfMerit("CU:P2");
				search.setExplorationMethod("fast-CBC");
				search.setNormType("2");
				search.addWeight(w);
				Experiment exp = new Experiment(experimentName, search, integrand, randomization);
				exp.setVarianceMC(varMC);
				
				exp.setSearch(search);
				exp.toJson(experimentDirPath, experimentHashes);
			}

			for (String randomization : new String[]{"LMS - Interlace - Shift", "NUS - Interlace"})
			{
				PolynomialLatticeSearch search = initializePLRSearch(s, k, pathToLNB, "lattice");
				search.setFigureOfMerit("CU:IC2");
				search.setInterlacing("2");
				search.setExplorationMethod("fast-CBC");
				search.setNormType("1");
				search.addWeight(w);
				Experiment exp = new Experiment(experimentName, search, integrand, randomization);
				exp.setSearch(search);
				exp.toJson(experimentDirPath, experimentHashes);
			}

			for (String randomization : new String[]{"LMS - Interlace - Shift", "NUS - Interlace"})
			{
				PolynomialLatticeSearch search = initializePLRSearch(s, k, pathToLNB, "lattice");
				search.setFigureOfMerit("CU:IC3");
				search.setInterlacing("3");
				search.setExplorationMethod("fast-CBC");
				search.setNormType("1");
				search.addWeight(w);
				Experiment exp = new Experiment(experimentName, search, integrand, randomization);
				exp.setSearch(search);
				exp.toJson(experimentDirPath, experimentHashes);
			}

		}
	}

	private static void runExperiments(Path experimentDirPath, Experiment[] constructedExperiments, boolean redoRandomization) throws IOException {
		for (int i = 0; i < constructedExperiments.length; i++){
			constructedExperiments[i].runExperiment(experimentDirPath, i, redoRandomization);
		}
	}

	/**
	 * Runs batches of RQMC experiments with LatNet Builder. One experiment corresponds to the choice of:
	 *   - an integrand.
	 *   - a point set (either pre-defined, e.g. Joe & Kuo's Sobol, or through a call to LatNet Builder).
	 *   - a randomization method.
	 * 
	 * The idea is that each experiment is saved in a JSON file. First, the JSON files are created. Then
	 * for each file, the experiment is computed.
	 * 
	 * Must be provided the following arguments:
	 * args[0]: name of the batch of experiments. (e.g. "myNewExperiment").
	 * args[1]: path to LatNet Builder executable (in string format).
	 * args[2]: path to result folder (in string format).
	 * args[3]: path to Joe & Kuo direction numbers file (in string format).
	 * args[4]: boolean, whether to create the batch of experiments.
	 * args[5]: boolean, whether to run the batch of experiments.
	 */
	public static void main(String[] args) throws IOException {
		String experimentName = args[0];
		String pathToLNB = args[1];
		String pathToResults = args[2];
		String pathToJoeKuo = args[3];
		boolean createExperimentFiles = Boolean.parseBoolean(args[4]);
		boolean runNewExperiments = Boolean.parseBoolean(args[5]);
		boolean runAllExperiments = Boolean.parseBoolean(args[6]);


		Path experimentDirPath = Paths.get(pathToResults, experimentName);
		if (!experimentDirPath.toFile().exists()) {
			experimentDirPath.toFile().mkdirs();
		}
		Experiment[] constructedExperiments = Experiment.fromJSONs(experimentDirPath);
		Set<Integer> experimentHashes = new HashSet<Integer>();
		for (final Experiment exp : constructedExperiments){
			experimentHashes.add(exp.getHash());
		}

		if (createExperimentFiles){
			createExperimentFiles(experimentDirPath, pathToLNB, pathToJoeKuo, experimentName, experimentHashes);
		}
		
		if (runNewExperiments){
			runExperiments(experimentDirPath, constructedExperiments, false);
		}

		if (runAllExperiments){
			runExperiments(experimentDirPath, constructedExperiments, true);
		}
	}
}
