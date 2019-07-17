package umontreal.ssj.markovchainrqmc.florian.chemicalReactionNetworks;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import umontreal.ssj.hups.CachedPointSet;
import umontreal.ssj.hups.PointSet;
import umontreal.ssj.hups.PointSetIterator;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.RQMCPointSet;
import umontreal.ssj.markovchainrqmc.ArrayOfComparableChains;
import umontreal.ssj.markovchainrqmc.MarkovChainComparable;
import umontreal.ssj.stat.PgfDataTable;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.Num;
import umontreal.ssj.util.sort.MultiDimSort;
import umontreal.ssj.util.sort.florian.MultiDim;

/**
 * This class offers a subset of the methods of #MarkovChainComparable  but allows to have a different sorting algorithm in every step.
 * @remark **Florian:** It really does not do anything very new. So, I skipped the documentation for this one.
 * @author florian
 *
 * @param <T>
 */
public class ArrayOfComparableChainsMultipleSorts<T extends MarkovChainComparable> extends ArrayOfComparableChains<T> {
	
	/**
	 * List containing one sort for each step.
	 */
	List<MultiDimSort<MultiDim>> sortList;
	
	public ArrayOfComparableChainsMultipleSorts(T baseChain, List<MultiDimSort<MultiDim>> sortList) {
		super(baseChain);
		this.sortList = sortList;
	}
	
	public ArrayOfComparableChainsMultipleSorts(T baseChain, PointSetRandomization rand, List<MultiDimSort<MultiDim>> sortList) {
		super(baseChain);
		this.sortList = sortList;
		this.randomization = rand;
	}
	
	//Same as usual, but removed sort
		public double simulArrayRQMC(PointSet p, PointSetRandomization rand,  int sortCoordPts,
				int numSteps) throws IOException {
		
			initialStates();
			int numNotStopped = n;
			int dim = chains[0].stateDim;
			int step = 0;
			MultiDimSort<MultiDim> sort;
			while (step < numSteps && numNotStopped > 0) {
				sort = sortList.get(step);
//				System.out.println("Network for step " + step + "loaded!"); //test
				if (numNotStopped == n) {

					sort.sort((MultiDim[]) chains, 0, n); 
//					sort.sort((MultiDim[]) chains);

					
				}
				else
					sortNotStoppedChains(sort); // Sort the numNotStopped first chains.
				
				p.randomize(rand); // Randomize the point set.
				if (sortCoordPts > 0) {
					if (!(p instanceof CachedPointSet))
						throw new IllegalArgumentException("p is not a CachedPointSet.");
					if (sortCoordPts > 1)
						((CachedPointSet) p).sort(sort); // Sort points using first sortCoordPts coordinates.
					else
						((CachedPointSet) p).sortByCoordinate(0); // Sort by first coordinate.
				}
				PointSetIterator stream = p.iterator();
				stream.resetCurPointIndex(); // Go to first point.
				int i = 0;
				for (T mc : chains) { // Assume the chains are sorted
					if (mc.hasStopped()) {
						numNotStopped--;
					} else {
						stream.setCurCoordIndex(sortCoordPts); // Skip first sortCoordPts coord.
						mc.nextStep(stream); // simulate next step of the chain.
						stream.resetNextSubstream(); // Go to next point.
						if (mc.hasStopped())
							numNotStopped--;
					}
					performances[i] = mc.getPerformance();
					++i;
				}
				performancePerRun[step].add( Arrays.stream(performances).average().orElse(Double.NaN) );
				++step;
			}
			// System.out.println("calcMeanPerf"+calcMeanPerf());
			return calcMeanPerf();
			}
		
		//NEW BY FLORIAN
		public double simulArrayRQMC(PointSet p, PointSetRandomization rand,  int sortCoordPts, int coordsForSort,
				int numSteps) throws IOException {
		
			initialStates();
			int numNotStopped = n;
			int dim = chains[0].stateDim;
			int step = 0;
			MultiDimSort<MultiDim> sort;
			while (step < numSteps && numNotStopped > 0) {
				sort = sortList.get(step);
//				System.out.println("Network for step " + step + "loaded!"); //test
				if (numNotStopped == n) {

					sort.sort((MultiDim[]) chains, 0, n); 
//					sort.sort((MultiDim[]) chains);

					
				}
				else
					sortNotStoppedChains(sort); // Sort the numNotStopped first chains.
				
				p.randomize(rand); // Randomize the point set.
				if (sortCoordPts > 0) {
					if (!(p instanceof CachedPointSet))
						throw new IllegalArgumentException("p is not a CachedPointSet.");
					if (sortCoordPts > 1)
						((CachedPointSet) p).sort(sort); // Sort points using first sortCoordPts coordinates.
					else
						((CachedPointSet) p).sortByCoordinate(0); // Sort by first coordinate.
				}
				PointSetIterator stream = p.iterator();
				stream.resetCurPointIndex(); // Go to first point.
				int i = 0;
				for (T mc : chains) { // Assume the chains are sorted
					if (mc.hasStopped()) {
						numNotStopped--;
					} else {
						stream.setCurCoordIndex(coordsForSort); // Skip first sortCoordPts coord.
						mc.nextStep(stream); // simulate next step of the chain.
						stream.resetNextSubstream(); // Go to next point.
						if (mc.hasStopped())
							numNotStopped--;
					}
					performances[i] = mc.getPerformance();
					++i;
				}
				++step;
			}
			// System.out.println("calcMeanPerf"+calcMeanPerf());
			return calcMeanPerf();
			}
		
		//Same as usual, but removed sort
		public void simulReplicatesArrayRQMC(PointSet p, PointSetRandomization rand, int sortCoordPts,
				int numSteps, int m, Tally statReps) throws IOException {
			makeCopies(p.getNumPoints());
			statReps.init();
			for (int rep = 0; rep < m; rep++) {
				statReps.add(simulArrayRQMC(p, rand, sortCoordPts, numSteps));
			}
		}
		
		public void simulReplicatesArrayRQMC(PointSet p, PointSetRandomization rand, int sortCoordPts, int coordsForSort,
				int numSteps, int m, Tally statReps) throws IOException {
			makeCopies(p.getNumPoints());
			statReps.init();
			for (int rep = 0; rep < m; rep++) {
				statReps.add(simulArrayRQMC(p, rand, sortCoordPts, coordsForSort, numSteps));
			}
		}
		
		public String testVarianceRateFormat(RQMCPointSet[] rqmcPts, 
				int sortCoordPts, int numSteps, int m, double varMC, String filenamePlot, String methodLabel) throws IOException {
			int numSets = rqmcPts.length; // Number of point sets.
			Tally statPerf = new Tally("Performance");
			double[] logn = new double[numSets];
			double[] variance = new double[numSets];
			double[] logVariance = new double[numSets];
			long initTime; // For timings.

			StringBuffer str = new StringBuffer("\n\n --------------------------");
			str.append(methodLabel + "\n  MC Variance : " + varMC + "\n\n");

			// Array-RQMC experiment with each pointSet.
			for (int i = 0; i < numSets; ++i) {
				for(int j = 0; j < numSteps; ++j)
					performancePerRun[j] = new Tally();
				initTime = System.currentTimeMillis();
				n = rqmcPts[i].getNumPoints();
				str.append("n = " + n + "\n");
				simulReplicatesArrayRQMC(rqmcPts[i].getPointSet(), rqmcPts[i].getRandomization(), sortCoordPts, numSteps, m, statPerf);
				logn[i] = Num.log2(n);
				variance[i] = statPerf.variance();
				logVariance[i] = Num.log2(variance[i]);
				str.append("  Average = " + statPerf.average() + "\n");
				str.append(" RQMC Variance : " +  variance[i] + "\n\n");
				str.append("  VRF =  " + varMC / (n * variance[i]) + "\n");
				str.append(formatTime((System.currentTimeMillis() - initTime) / 1000.) + "\n");
			}
			// Estimate regression slope and print plot and overall results.
			double regSlope = slope(logn, logVariance, numSets);
			str.append("Regression slope (log) for variance = " + regSlope + "\n\n");

			String[] tableField = { "log(n)", "log(Var)" };
			double[][] data = new double[numSets][2];
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				data[s][0] = logn[s];
				data[s][1] = logVariance[s];
			}

			// Print plot and overall results in files.
			if (filenamePlot != null)
				try {
					PgfDataTable pgf = new PgfDataTable(filenamePlot, rqmcPts[0].getLabel(), tableField, data);
					String pVar = pgf.drawPgfPlotSingleCurve(filenamePlot, "axis", 0, 1, 2, "", "");
					String plotIV = (PgfDataTable.pgfplotFileHeader() + pVar + PgfDataTable.pgfplotEndDocument());

					FileWriter fileIV = new FileWriter(filenamePlot + "_" + "VAr.tex");
					fileIV.write(plotIV);
					fileIV.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			
			tableField = new String[] {"step", "average", "variance"};
			data = new double[numSteps][3];
			for (int s = 0; s < numSteps; s++) { // For each cardinality n
				data[s][0] = s+1;
				data[s][1] = performancePerRun[s].average();
				data[s][2] = performancePerRun[s].variance();
			}
			
			if (filenamePlot != null)
				try {

					PgfDataTable pgf = new PgfDataTable(filenamePlot, rqmcPts[0].getLabel(), tableField, data);
					String pMean = pgf.drawPgfPlotSingleCurve(filenamePlot, "axis", 0, 1, 2, "", "");
					String pVar = pgf.drawPgfPlotSingleCurve(filenamePlot, "axis", 0, 2, 2, "", "");
					
					String plotMean = (PgfDataTable.pgfplotFileHeader() + pMean + PgfDataTable.pgfplotEndDocument());
					String plotVar = (PgfDataTable.pgfplotFileHeader() + pVar + PgfDataTable.pgfplotEndDocument());

					FileWriter fileIV = new FileWriter(filenamePlot + "_" + "MeanPerStep.tex");
					fileIV.write(plotMean);
					fileIV.close();
					
					fileIV = new FileWriter(filenamePlot + "_" + "VariancePerStep.tex");
					fileIV.write(plotVar);
					fileIV.close();
				} catch (IOException e) {
					e.printStackTrace();
				}


			return str.toString();
		}
		
		//NEW BY FLORIAN
		public String testVarianceRateFormat(RQMCPointSet[] rqmcPts, 
				int sortCoordPts, int coordsForSort, int numSteps, int m, double varMC, String filenamePlot, String methodLabel) throws IOException {
			int numSets = rqmcPts.length; // Number of point sets.
			Tally statPerf = new Tally("Performance");
			double[] logn = new double[numSets];
			double[] variance = new double[numSets];
			double[] logVariance = new double[numSets];
			long initTime; // For timings.

			StringBuffer str = new StringBuffer("\n\n --------------------------");
			str.append(methodLabel + "\n  MC Variance : " + varMC + "\n\n");

			// Array-RQMC experiment with each pointSet.
			for (int i = 0; i < numSets; ++i) {
				for(int j = 0; j < numSteps; ++j)
					performancePerRun[j] = new Tally();
				initTime = System.currentTimeMillis();
				n = rqmcPts[i].getNumPoints();
				str.append("n = " + n + "\n");
				simulReplicatesArrayRQMC(rqmcPts[i].getPointSet(), rqmcPts[i].getRandomization(), sortCoordPts, coordsForSort, numSteps, m, statPerf);
				logn[i] = Num.log2(n);
				variance[i] = statPerf.variance();
				logVariance[i] = Num.log2(variance[i]);
				str.append("  Average = " + statPerf.average() + "\n");
				str.append(" RQMC Variance : " +  variance[i] + "\n\n");
				str.append("  VRF =  " + varMC / (n * variance[i]) + "\n");
				str.append(formatTime((System.currentTimeMillis() - initTime) / 1000.) + "\n");
			}
			// Estimate regression slope and print plot and overall results.
			double regSlope = slope(logn, logVariance, numSets);
			str.append("Regression slope (log) for variance = " + regSlope + "\n\n");

			String[] tableField = { "log(n)", "log(Var)" };
			double[][] data = new double[numSets][2];
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				data[s][0] = logn[s];
				data[s][1] = logVariance[s];
			}

			// Print plot and overall results in files.
			if (filenamePlot != null)
				try {
					PgfDataTable pgf = new PgfDataTable(filenamePlot, rqmcPts[0].getLabel(), tableField, data);
					String pVar = pgf.drawPgfPlotSingleCurve(filenamePlot, "axis", 0, 1, 2, "", "");
					String plotIV = (PgfDataTable.pgfplotFileHeader() + pVar + PgfDataTable.pgfplotEndDocument());

					FileWriter fileIV = new FileWriter(filenamePlot + "_" + "VAr.tex");
					fileIV.write(plotIV);
					fileIV.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			
			tableField = new String[] {"step", "average", "variance"};
			data = new double[numSteps][3];
			for (int s = 0; s < numSteps; s++) { // For each cardinality n
				data[s][0] = s+1;
				data[s][1] = performancePerRun[s].average();
				data[s][2] = performancePerRun[s].variance();
			}
			
			if (filenamePlot != null)
				try {

					PgfDataTable pgf = new PgfDataTable(filenamePlot, rqmcPts[0].getLabel(), tableField, data);
					String pMean = pgf.drawPgfPlotSingleCurve(filenamePlot, "axis", 0, 1, 2, "", "");
					String pVar = pgf.drawPgfPlotSingleCurve(filenamePlot, "axis", 0, 2, 2, "", "");
					
					String plotMean = (PgfDataTable.pgfplotFileHeader() + pMean + PgfDataTable.pgfplotEndDocument());
					String plotVar = (PgfDataTable.pgfplotFileHeader() + pVar + PgfDataTable.pgfplotEndDocument());

					FileWriter fileIV = new FileWriter(filenamePlot + "_" + "MeanPerStep.tex");
					fileIV.write(plotMean);
					fileIV.close();
					
					fileIV = new FileWriter(filenamePlot + "_" + "VariancePerStep.tex");
					fileIV.write(plotVar);
					fileIV.close();
				} catch (IOException e) {
					e.printStackTrace();
				}


			return str.toString();
		}
		
		public String testVarianceRateFormat(PointSet[] pointSets, PointSetRandomization rand, 
				int sortCoordPts, int numSteps, int m, double varMC, String filenamePlot, String methodLabel) throws IOException {
			int numSets = pointSets.length; // Number of point sets.
			Tally statPerf = new Tally("Performance");
			double[] logn = new double[numSets];
			double[] variance = new double[numSets];
			double[] logVariance = new double[numSets];
			long initTime; // For timings.

			StringBuffer str = new StringBuffer("\n\n --------------------------");
			str.append(methodLabel + "\n  MC Variance : " + varMC + "\n\n");

			// Array-RQMC experiment with each pointSet.
			for (int i = 0; i < numSets; ++i) {
				for(int j = 0; j < numSteps; ++j)
					performancePerRun[j] = new Tally();
				initTime = System.currentTimeMillis();
				n = pointSets[i].getNumPoints();
				str.append("n = " + n + "\n");
				simulReplicatesArrayRQMC(pointSets[i], rand, sortCoordPts, numSteps, m, statPerf);
				logn[i] = Num.log2(n);
				variance[i] = statPerf.variance();
				logVariance[i] = Num.log2(variance[i]);
				str.append("  Average = " + statPerf.average() + "\n");
				str.append(" RQMC Variance : " + n * variance[i] + "\n\n");
				str.append("  VRF =  " + varMC / (n * variance[i]) + "\n");
				str.append(formatTime((System.currentTimeMillis() - initTime) / 1000.) + "\n");
			}
			// Estimate regression slope and print plot and overall results.
			double regSlope = slope(logn, logVariance, numSets);
			str.append("Regression slope (log) for variance = " + regSlope + "\n\n");

			String[] tableField = { "log(n)", "log(Var)" };
			double[][] data = new double[numSets][2];
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				data[s][0] = logn[s];
				data[s][1] = logVariance[s];
			}

			// Print plot and overall results in files.
			if (filenamePlot != null)
				try {
					/*
					 * Writer file = new FileWriter (filenamePlot + ".tex"); XYLineChart chart = new
					 * XYLineChart(); // ("title", "$log_2(n)$", "$log_2 Var[hat mu_{rqmc,s,n}]$");
					 * chart.add (logn, logVariance); file.write (chart.toLatex(12, 8));
					 * file.close();
					 */
					PgfDataTable pgf = new PgfDataTable(filenamePlot, pointSets.toString(), tableField, data);
					String pVar = pgf.drawPgfPlotSingleCurve(filenamePlot, "axis", 0, 1, 2, "", "");
					String plotIV = (PgfDataTable.pgfplotFileHeader() + pVar + PgfDataTable.pgfplotEndDocument());

					FileWriter fileIV = new FileWriter(filenamePlot + "_" + "VAr.tex");
					fileIV.write(plotIV);
					fileIV.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			
			tableField = new String[] {"step", "average", "variance"};
			data = new double[numSteps][3];
			for (int s = 0; s < numSteps; s++) { // For each cardinality n
				data[s][0] = s+1;
				data[s][1] = performancePerRun[s].average();
				data[s][2] = performancePerRun[s].variance();
			}
			
			if (filenamePlot != null)
				try {

					PgfDataTable pgf = new PgfDataTable(filenamePlot, methodLabel, tableField, data);
					String pMean = pgf.drawPgfPlotSingleCurve(filenamePlot, "axis", 0, 1, 2, "", "");
					String pVar = pgf.drawPgfPlotSingleCurve(filenamePlot, "axis", 0, 2, 2, "", "");
					
					String plotMean = (PgfDataTable.pgfplotFileHeader() + pMean + PgfDataTable.pgfplotEndDocument());
					String plotVar = (PgfDataTable.pgfplotFileHeader() + pVar + PgfDataTable.pgfplotEndDocument());

					FileWriter fileIV = new FileWriter(filenamePlot + "_" + "MeanPerStep.tex");
					fileIV.write(plotMean);
					fileIV.close();
					
					fileIV = new FileWriter(filenamePlot + "_" + "VariancePerStep.tex");
					fileIV.write(plotVar);
					fileIV.close();
				} catch (IOException e) {
					e.printStackTrace();
				}


			return str.toString();
		}
}
