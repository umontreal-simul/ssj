package latnetbuilder;

import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.SobolSequence;
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

		for (int effectDim = 1; effectDim <= 1; effectDim += 3) {
			for (int s = 3; s <= 3; s += 10) {
				for (int k = 3; k <= 5; k += 1) {

					double c = (double) effectDim / s;
					double w = c * c;
					SobolTestFunc integrand = new SobolTestFunc(c, s);

					int n = 10000;
					Tally statValue = new Tally("MC method");
					integrand.simulateRuns(n, new MRG32k3a(), statValue);
					double varMC = statValue.variance();

					{
						PolynomialLatticeSearch search = initializePLRSearch(s, k, pathToLNB, "net");
						search.setFigureOfMerit("CU:P2");
						search.setNormType("2");
						search.addWeight("product:" + w);
						Experiment exp = new Experiment(experimentName, search, integrand, "LMS - Shift");
						exp.setVarianceMC(varMC);

						search.setExplorationMethod("random:1");
						exp.setSearch(search);
						exp.toJson(experimentDirPath, experimentHashes);

						search.setExplorationMethod("random:100");
						exp.setSearch(search);
						exp.toJson(experimentDirPath, experimentHashes);

						search.changePointSetTypeView("lattice");
						search.setExplorationMethod("fast-CBC");
						exp.setSearch(search);
						exp.toJson(experimentDirPath, experimentHashes);
					}

					{
						PolynomialLatticeSearch search = initializePLRSearch(s, k, pathToLNB, "net");
						search.setFigureOfMerit("CU:IC2");
						search.setNormType("1");
						search.addWeight("product:" + w);
						search.setExplorationMethod("random:100");
						Experiment exp = new Experiment(experimentName, search, integrand, "NUS - Interlace");
						exp.setVarianceMC(varMC);

						search.setInterlacing("2");
						exp.setSearch(search);
						exp.toJson(experimentDirPath, experimentHashes);

						search.changePointSetTypeView("lattice");
						search.setExplorationMethod("fast-CBC");
						search.setInterlacing("2");
						exp.setSearch(search);
						exp.toJson(experimentDirPath, experimentHashes);
					}

					{
						PolynomialLatticeSearch search = initializePLRSearch(s, k, pathToLNB, "net");
						search.setFigureOfMerit("CU:IC2");
						search.setNormType("1");
						search.addWeight("product:" + w);
						search.setExplorationMethod("random:100");
						Experiment exp = new Experiment(experimentName, search, integrand, "LMS - Interlace - Shift");
						exp.setVarianceMC(varMC);

						search.setInterlacing("2");
						exp.setSearch(search);
						exp.toJson(experimentDirPath, experimentHashes);

						search.setInterlacing("3");
						exp.setSearch(search);
						exp.toJson(experimentDirPath, experimentHashes);

						search.changePointSetTypeView("lattice");
						search.setExplorationMethod("fast-CBC");
						search.setInterlacing("2");
						exp.setSearch(search);
						exp.toJson(experimentDirPath, experimentHashes);

						search.setInterlacing("3");
						exp.setSearch(search);
						exp.toJson(experimentDirPath, experimentHashes);
					}

					{
						DigitalNetBase2 q = new SobolSequence(pathToJoeKuo, k, DigitalNetBase2.getMaxBits(), s);
						Experiment exp = new Experiment(experimentName, q, integrand, "LMS - Shift");
						exp.setVarianceMC(varMC);
						exp.toJson(experimentDirPath, experimentHashes);
					}

					{
						DigitalNetBase2 q = new SobolSequence(pathToJoeKuo, k, DigitalNetBase2.getMaxBits(), 2 * s);
						q.setInterlacing(2);
						Experiment exp = new Experiment(experimentName, q, integrand, "LMS - Interlace - Shift");
						exp.setVarianceMC(varMC);
						exp.toJson(experimentDirPath, experimentHashes);
					}

					{
						DigitalNetBase2 q = new SobolSequence(pathToJoeKuo, k, DigitalNetBase2.getMaxBits(), 3 * s);
						q.setInterlacing(3);
						Experiment exp = new Experiment(experimentName, q, integrand, "LMS - Interlace - Shift");
						exp.setVarianceMC(varMC);
						exp.toJson(experimentDirPath, experimentHashes);
					}
				}
			}
		}
	}

	private static void runExperiments(Path experimentDirPath, Experiment[] constructedExperiments) throws IOException {
		for (int i = 0; i < constructedExperiments.length; i++){
			constructedExperiments[i].runExperiment(experimentDirPath, i);
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
		boolean runExperiments = Boolean.parseBoolean(args[5]);


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
		
		if (runExperiments){
			runExperiments(experimentDirPath, constructedExperiments);
		}
	}
}
