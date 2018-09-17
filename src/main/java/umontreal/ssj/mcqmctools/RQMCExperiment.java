package umontreal.ssj.mcqmctools;

import umontreal.ssj.hups.*;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.*;
import umontreal.ssj.stat.list.ListOfTallies;
import umontreal.ssj.stat.list.lincv.ListOfTalliesWithCV;
import umontreal.ssj.util.Chrono;
import umontreal.ssj.util.PrintfFormat;

/**
 * Provides basic generic tools to perform RQMC experiments
 * with a simulation model that implements the MonteCarloModelDouble interface.
 */

/**
 * @author Pierre L'Ecuyer
 * 
 */
public class RQMCExperiment extends MonteCarloExperiment {

	/**
	 * Simulate <tt>m</tt> replications with #prqmc and return the #m RQMC average
	 * observations in statReps. These <tt>m</tt> observations are usually
	 * independent if the randomizations used in the RQMC point set #prqmc are
	 * independent.
	 */
	public static void simulReplicatesRQMC(MonteCarloModelDouble model, RQMCPointSet prqmc, int m, Tally statReps) {
		simulReplicatesRQMC(model, prqmc.getPointSet(), prqmc.getRandomization(), m, statReps);
	}

	/**
	 * Simulate <tt>m</tt> replications and return the #m RQMC average observations
	 * in statReps. Here the QMC point set #p and its randomization #rand are
	 * specified directly.
	 */
	public static void simulReplicatesRQMC(MonteCarloModelDouble model, PointSet p, PointSetRandomization rand, int m,
			Tally statReps) {
		statReps.init();
		int n = p.getNumPoints();
		// Internal collector for stats on the n outputs X, for each replication.
		Tally statValue = new Tally();
		PointSetIterator stream = p.iterator();
		for (int rep = 0; rep < m; rep++) {
			rand.randomize(p);
			stream.resetStartStream();
			simulateRuns(model, n, stream, statValue);
			statReps.add(statValue.average()); // For the estimator of the mean.
		}
	}

	/**
	 * Same as @ref simulReplicatesRQMC, except that all the n observations for each
	 * the m replications are saved and returned in a new two-dimensional
	 * array @data, which is an array of <tt>m</tt> arrays of size <tt>n</tt>, i.e.,
	 * <tt>double[m][n]</tt>, created inside this method. Each array of size
	 * <tt>n</tt> is sorted by increasing order. This is useful for density and cdf
	 * estimation, or for further processing of the data, for example.
	 */
	public static void simulReplicatesRQMC(MonteCarloModelDouble model, RQMCPointSet prqmc, int m, Tally statReps,
			double[][] data) {
		simulReplicatesRQMC(model, prqmc.getPointSet(), prqmc.getRandomization(), m, statReps, data);
	}

	/**
	 * Here the QMC point set #p and its randomization #rand are specified directly.
	 */
	public static void simulReplicatesRQMC(MonteCarloModelDouble model, PointSet p, PointSetRandomization rand, int m,
			Tally statReps, double[][] data) {
		int n = p.getNumPoints();
		data = new double[m][]; // ????
		// Internal collector for stats on the n outputs X, for each replication.
		TallyStore statSave = new TallyStore(n);
		PointSetIterator stream = p.iterator();
		for (int rep = 0; rep < m; rep++) {
			rand.randomize(p);
			stream.resetStartStream();
			simulateRuns(model, n, stream, statSave);
			statReps.add(statSave.average()); // For the estimator of the mean.
			statSave.quickSort();
			data[rep] = statSave.getArray(); // Instead of copy, just exchange pointers! Good?
			statSave = new TallyStore(n);
		}
	}
	

	public static String simulReplicatesRQMCDefaultReport (MonteCarloModelDouble model, 
			PointSet p, PointSetRandomization rand, int m, RandomStream noise,
			Tally statRQMC) {
		PrintfFormat str = new PrintfFormat();
		Chrono timer = new Chrono();
		simulReplicatesRQMC(model, p, rand, m, statRQMC);
		statRQMC.setConfidenceIntervalStudent();
		str.append (model.toString());
		str.append (p.toString());
		str.append (rand.toString());
		str.append (statRQMC.report(0.95, 4));
		str.append ("Total CPU time:      " + timer.format() + "\n");
		return str.toString();
	}

	public static String simulReplicatesRQMCDefaultReportCompare (MonteCarloModelDouble model,
			PointSet p, PointSetRandomization rand, int m,
			Tally statRQMC, double varianceMC, double secondsMC) {
		PrintfFormat str = new PrintfFormat();
		Chrono timer = new Chrono();
		simulReplicatesRQMC(model, p, rand, m, statRQMC);
		double secondsRQMC = timer.getSeconds() / (m * p.getNumPoints());
		double varianceRQMC = p.getNumPoints() * statRQMC.variance();
		statRQMC.setConfidenceIntervalStudent();
		str.append (model.toString());
		str.append (p.toString());
		str.append (rand.toString());
		str.append (statRQMC.report(0.95, 4));
		str.append ("Total CPU time:      " + timer.format() + "\n");
		str.append ("Variance per run: ");
		str.append(12, 5, 4, varianceRQMC);
		str.append ("Variance ratio:   ");
		str.append(12, 5, 4, varianceMC / varianceRQMC);
		str.append ("Efficiency ratio: ");
		str.append(12, 5, 4, (varianceMC * secondsMC) / (varianceRQMC * secondsRQMC));
		str.append ("-------------------------------------------\n");
		return str.toString();
	}
	
	/**
	 * Similar to
	 * {@link #simulReplicatesRQMC(MonteCarloModelDouble, RQMCPointSet, int, Tally)}
	 * but for a model of type @ref MonteCarloModelDoubleArray. Consequently, the
	 * statistics are collected in a @ref ListOfTallies. The \f$t\f$-th element of
	 * #statValueList collects the statistics for the \f$t\f$-th coordinate of the
	 * performance vector of #model.
	 * 
	 * @param model
	 *            the underlying model which is simulated.
	 * @param prqmc
	 *            the RQMC-point set used.
	 * @param m
	 *            number of independent replications.
	 * @param statReps
	 *            statistical container collecting the obtained estimates.
	 */
	public static void simulReplicatesRQMC(MonteCarloModelDoubleArray model, RQMCPointSet prqmc, int m,
			ListOfTallies<Tally> statRepsList) {
		simulReplicatesRQMC(model, prqmc.getPointSet(), prqmc.getRandomization(), m, statRepsList);
	}

	/**
	 * Same as
	 * {@link #simulReplicatesRQMC(MonteCarloModelArrayOfDoubles, RQMCPointSet, int, ListOfTallies)}
	 * but with the point set #p and its randomization #rand specified directly.
	 * 
	 * @param model
	 *            the underlying model which is simulated.
	 * @param p
	 *            the point set used.
	 * @param rand
	 *            the point set randomization used.
	 * @param m
	 *            number of independent replications.
	 * @param statRepsList
	 *            statistical container collecting the obtained estimates.
	 */
	public static void simulReplicatesRQMC(MonteCarloModelDoubleArray model, PointSet p, PointSetRandomization rand,
			int m, ListOfTallies<Tally> statRepsList) {
		statRepsList.init();
		int n = p.getNumPoints();
		int t = model.getPerformanceDim();
		// Internal collector for stats on the n outputs X, for each
		// replication.
		ListOfTallies<Tally> statValue = ListOfTallies.createWithTally(t);
		PointSetIterator stream = p.iterator();
		for (int rep = 0; rep < m; rep++) {
			rand.randomize(p);
			stream.resetStartStream();
			simulateRuns(model, n, stream, statValue);
			double[] means = new double[t];
			statValue.average(means);
			statRepsList.add(means); // For the estimator of the mean.
		}
	}

	/**
	 * Same as
	 * {@link #simulReplicatesRQMC(MonteCarloModelDoubleArray, RQMCPointSet, int, ListOfTallies)},
	 * except that all the \f$n\f$ observations of dimension \f$t\f$ for each of
	 * the @f$m@f$ replications are saved and returned in a new three-dimensional
	 * array <tt>data<\tt>, which is 
	 * an array of size \f$m\times n\times t\f$, i.e., <tt>double[m][n][t]</tt>. 
	 * 
	 * Note that the \f$n\times t\f$ submatrices of <tt>data</tt> are not sorted.
	 * 
	 * @param model
	 *            the underlying model which is simulated.
	 * @param prqmc
	 *            the RQMC-point set used.
	 * @param m
	 *            number of independent replications.
	 * @param statRepsList
	 *            statistical container collecting the obtained estimates.
	 * @param data
	 *            array in which all the observations are stored.
	 */
	public static void simulReplicatesRQMC(MonteCarloModelDoubleArray model, RQMCPointSet prqmc, int m,
			ListOfTallies<Tally> statRepsList, double[][][] data) {
		simulReplicatesRQMC(model, prqmc.getPointSet(), prqmc.getRandomization(), m, statRepsList, data);

	}

	/**
	 * Same as {@link #simulReplicatesRQMC(MonteCarloModelDoubleArray, RQMCPointSet, int, ListOfTallies, double[][][])},
	 * but with the point set and the randomization  passed separately.
	 * 
	 * @param model
	 *            the underlying model which is simulated.
	 * @param p
	 *            the QMC-point set used.
	 * @param rand
	 *            the point set randomization used.
	 * @param m
	 *            number of independent replications.
	 * @param statRepsList
	 *            statistical container collecting the obtained estimates.
	 * @param data
	 *            array in which all the observations are stored.
	 */

	public static void simulReplicatesRQMC(MonteCarloModelDoubleArray model, PointSet p, PointSetRandomization rand,
			int m, ListOfTallies<Tally> statRepsList, double[][][] data) {
		int n = p.getNumPoints();
		int t = model.getPerformanceDim();
		// Internal collector for stats on the n outputs X, for each
		// replication.
		ListOfTallies<TallyStore> statSave = ListOfTallies.createWithTallyStore(n, t);
		PointSetIterator stream = p.iterator();
		for (int rep = 0; rep < m; rep++) {		
				p.randomize(rand);
			stream.resetStartStream();
			simulateRuns(model, n, stream, statSave);
			// TODO: should we keep the possibility to sort the data?
			double[] means = new double[t];
			statSave.average(means);
			statRepsList.add(means); // For the estimator of the mean.
			for (int i = 0; i < t; i++)
				// This allows to just exchange pointers instead of copying
				// entries. But the array will be [m][t][n]
				// which is less practical than [m][n][t]
				// data[rep][i] = statSave.get(i).getArray();
				for (int k = 0; k < n; k++)
					// creates [m][n][t]-array
					data[rep][k][i] = statSave.get(i).getArray()[k];
			statSave = ListOfTallies.createWithTallyStore(n, t);

		}

	}

	/**
	 * Same as @ref simulReplicatesRQMC, except that this one uses control variates.
	 * It returns in <TT>statWithCV</TT> the statistics for m observations which
	 * corresponds to the m RQMC replications. Each observation is a vector that
	 * contains the average over the n simulation runs of the performance value and
	 * of the control variates. This implements the replication method for RQMC with
	 * CVs, as in Section 5 of Hickernell, Lemieux and Owen (2005). The CV
	 * estimates, variances, and covariances can be recovered from
	 * <TT>statWithCV</TT>.
	 */
	public static void simulReplicatesRQMCCV(MonteCarloModelCV model, RQMCPointSet prqmc, int m,
			ListOfTalliesWithCV<Tally> statWithCV) {
		simulReplicatesRQMCCV(model, prqmc.getPointSet(), prqmc.getRandomization(), m, statWithCV);
	}

	/**
	 * Same as @ref simulReplicatesRQMCCV above, except that here the point set #p
	 * and the randomization #rand are specified directly.
	 */
	public static void simulReplicatesRQMCCV(MonteCarloModelCV model, PointSet p, PointSetRandomization rand, int m,
			ListOfTalliesWithCV<Tally> statWithCV) {
		statWithCV.init();
		int n = p.getNumPoints();
		int numCV = model.getNumberCV();
		PointSetIterator stream = p.iterator();
		for (int rep = 0; rep < m; rep++) {
			rand.randomize(p);
			stream.resetStartStream();
			double sumValues = 0.0;
			double[] sumValuesCV = new double[numCV];
			double[] curValuesCV = new double[numCV]; // Current value of the CV
			for (int i = 0; i < n; i++) {
				model.simulate(stream);
				sumValues += model.getPerformance();
				curValuesCV = model.getValuesCV();
				for (int k = 0; k < numCV; k++)
					sumValuesCV[k] += curValuesCV[k];
				stream.resetNextSubstream();
			}
			for (int k = 0; k < numCV; k++)
				sumValuesCV[k] /= (double) n;
			statWithCV.add(sumValues / (double) n, sumValuesCV);
		}
	}
	
	/**
	 * To estimate a derivative via a finite difference.
	 */	
	public static void simulFDReplicatesRQMC (MonteCarloModelDouble model1, MonteCarloModelDouble model2, double delta,
			PointSet p, PointSetRandomization rand, int m, Tally statDiffRQMC) {
		Tally statValue = new Tally("stat on value");
		statDiffRQMC.init();
		double average1; 
		PointSetIterator stream = p.iterator();
		for (int j = 0; j < m; j++) {
			rand.randomize(p);
			stream.resetStartStream();
			simulateRuns(model1, p.getNumPoints(), stream, statValue);
			average1 = statValue.average();
			stream.resetStartStream();
			simulateRuns(model2, p.getNumPoints(), stream, statValue);
			statDiffRQMC.add((statValue.average() - average1) / delta);
		}
	}
}
