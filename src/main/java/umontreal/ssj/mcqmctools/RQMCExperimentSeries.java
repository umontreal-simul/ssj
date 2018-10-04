/*
 * Class:        RQMCExperimentSeries
 * Description:  randomized quasi-Monte Carlo simulations
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author
 * @since
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package umontreal.ssj.mcqmctools;

import umontreal.ssj.functionfit.LeastSquares;
import umontreal.ssj.hups.*;
import umontreal.ssj.stat.PgfDataTable;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.list.lincv.ListOfTalliesWithCV;
import umontreal.ssj.util.Chrono;
import umontreal.ssj.util.PrintfFormat;
import java.util.ArrayList;


/**
 * This class offers facilities to perform experiments on the convergence
 * of the variance when estimating a mean (expectation) with a series of RQMC 
 * point sets usually of the same type, but different sizes @f$n@f$.
 * The series of RQMC point sets of different sizes can be passed in an array 
 * to the constructor. The method @ref testVarianceRate performs an experiment 
 * with a given @ref MonteCarloModelDouble and the series of point sets. 
 * In this experiment, for each size @f$n@f$ of point set, @f$m@f$ independent replicates 
 * of the RQMC estimator (which is an average over RQMC @f$n@f$ points) are computed.
 * One can then recover the average and the empirical variance of these @f$m@f$ replicates
 * of the @f$n@f$-point average, their logs in any base, etc., in arrays, 
 * as well as the estimated linear regression of log(variance) as a function of log(n). 
 * 
 * One can also perform experiments with many such series of point sets for the same model, 
 * and display the results for the different series on the same plot.
 * Likewise, one can also perform an experiment with many different estimators 
 * for the same model and same point sets.
 */

public class RQMCExperimentSeries {
	int numSets = 0;   // Number of point sets in the series.
    RQMCPointSet[] theSets;   
    double base = 2.0;    // Base for the logs (in base 2 by default)
    double logOfBase;       // Math.log(base)
	double[] size = new double[numSets];    // values of n
	double[] mean = new double[numSets];    // average performance for each point set 
	double[] variance = new double[numSets]; // variance of the average for each point set
	double[] logn = new double[numSets];   // log_base n 
	double[] logVar = new double[numSets]; // log_base (variance)
	String[] tableFields = {"n", "mean", "variance", "log(n)", "log(variance)"};
	                                       // Names of fields for table.
	boolean displayExec = false;   // When true, prints a display of execution in real time
	int numReplicates;    // last value of m
	MonteCarloModelDouble model;
	// int numSkipRegression = 0; // Number of values of n that are skipped for the regression
	String cpuTime;       // time for last experiment\
    String title;

   /**
    * Constructor with a give series of RQMC point sets.
    *  @param theSets      the RQMC point sets.
    *  @param base 		the base used for all logarithms.
    */
   public RQMCExperimentSeries (RQMCPointSet[] theSets, double base) {
	   init(theSets, base);
   }

   /**
    * Resets the array of RQMC point sets for this object, and initializes 
    * (re-creates) the arrays that will contain the results.
    */
   public void init(RQMCPointSet[] theSets, double base) {
	   this.base = base;
	   this.logOfBase = Math.log(base);
	     this.numSets = theSets.length;
	     this.theSets = theSets;
		 size = new double[numSets]; //  n for each point set
	     mean = new double[numSets]; //  average for each point set
	     variance = new double[numSets]; // variance for each point set
	     logn = new double[numSets];    // log n
	     logVar = new double[numSets];  // log (variance)
   }
   
   /**
    * When set to true, a real-time display of execution results and CPU times 
    * will be printed on the default output.
    */
   public void setExecutionDisplay (boolean display) {
	      this.displayExec = display;
   }
   
   /**
    * Sets the base for the logs to b.
    */
   public void setBase (double b) {
      base = b;
	  logOfBase = Math.log(base);
   }

   /**
    * Returns the base used for the logs.
    */
   public double getBase() {
      return base;
   }

   /**
    * Returns the point set number i associated to this object (starts at 0).
    *  @return the ith point set associated to this object
    */
   public RQMCPointSet getSet(int i) {
      return theSets[i];
   }

   /**
    * Returns the vector of values of n, after an experiment.
    */
   public double[] getValuesn() {
      return size;
   }

   /**
    * Returns the vector of log_base(n), after an experiment.
    */
   public double[] getLogn() {
      return logn;
   }

  /**
    * Returns the vector of means from last experiment, after an experiment.
    */
   public double[] getMeans() {
      return mean;
   }

   /**
    * Returns the vector of variances from the last experiment.  
    * Each variance in the vector is the empirical variance of the m replicates
    * of the RQMC estimator. To obtain the variance per run, it must be multiplied by n.
    */
   public double[] getVariances() {
      return variance;
   }

   /**
    * Returns the vector of log_base of variances from the last experiment.
    */
   public double[] getLogVariances() {
      return logVar;
   }

    /**
    * Performs an RQMC experiment with the given model, with this series of RQMC point sets.  
    * For each set in the series, computes m replicates of the RQMC estimator, 
    * the computes the average and the variance of these m replicates, 
    * and the logs of n and of the variance in the given base.
    */
   public void testVarianceRate (MonteCarloModelDouble model, int m) {
		int n;
		Tally statReps = new Tally();
		Chrono timer = new Chrono();
		numReplicates = m;
		this.model = model;
	    if (displayExec) {
			System.out.println("\n ============================================= ");
	    	System.out.println("RQMC simulation for mean estimation:  ");
	    	System.out.println("Model: " + model.toString());
	    	System.out.println(" Number of indep copies m  = " + m);
	    	System.out.println(" Point sets: " + theSets[0].toString() + "\n");
			System.out.println("    n     CPU time         mean      log(var) ");	    	
	    }
		for (int s = 0; s < numSets; s++) { // For each cardinality n
			n = theSets[s].getNumPoints();
			size[s] = n;
			logn[s] = Math.log(n) / logOfBase;
			// System.out.println(" n = " + n + ", log n = " + logn[s] + "\n"); // ****
			// System.out.println("  " + n + "     " + timer.format());
			RQMCExperiment.simulReplicatesRQMC (model, theSets[s], m, statReps);
			mean[s] = statReps.average();
			variance[s] = statReps.variance();
		    logVar[s] = Math.log(variance[s]) / logOfBase;
		    if (displayExec) {
			   System.out.println("  " + n + "     " + timer.format() + 
			              "   " + PrintfFormat.f(10, 5, mean[s]) + 
			              "   " + PrintfFormat.f(7, 2, logVar[s]));
		    }
		}	   
        cpuTime = timer.format();	 
   }

   /**
    * Similar to testVarianceRate, but with control variates, all centered at 0.
    */
   public void testVarianceRateCV (MonteCarloModelCV model, int m) {
		int numCV = model.getNumberCV();
		ListOfTalliesWithCV<Tally> statWithCV = ListOfTalliesWithCV.createWithTally(1, numCV);
		statWithCV.setExpectedValue (0, 0.0);  // The CV is centered to 0.
		int n;
		Chrono timer = new Chrono();
		numReplicates = m;
		this.model = model;
	    if (displayExec) {
	    	System.out.println("\n ============================================= ");
	    	System.out.println("RQMC simulation for mean estimation with control variates:  ");
	    	System.out.println("Model: " + model.toString());
	    	System.out.println(" Number of indep copies m  = " + m);
	    	System.out.println(" Point sets: " + theSets[0].toString() + "\n");
			System.out.println("    n     CPU time         mean     log(var) ");	    	
	    }
		for (int s = 0; s < numSets; s++) { // For each cardinality n
			n = theSets[s].getNumPoints();
			size[s] = n;
			logn[s] = Math.log(n) / logOfBase;
			// System.out.println(" n = " + n + ", log n = " + logn[s] + "\n"); // ****
			// System.out.println("  " + n + "     " + timer.format());
			RQMCExperiment.simulReplicatesRQMCCV (model, theSets[s], m, statWithCV);
			statWithCV.estimateBeta();    // This is where the var. and covar. are computed!
			mean[s] = statWithCV.averageWithCV(0);
			variance[s] = statWithCV.covarianceWithCV(0, 0);
		    logVar[s] = Math.log(variance[s]) / logOfBase;
		    if (displayExec) {
			   System.out.println("  " + n + "     " + timer.format() + 
			              "   " + PrintfFormat.f(10, 5, mean[s]) + 
					      "   " + PrintfFormat.f(7, 2, logVar[s]));
		    }
		}	   
        cpuTime = timer.format();	   
   }

   /**
    * Performs a linear regression of log(variance) vs log(n), and returns the 
    * coefficients (constant and slope) in two-dimensional vector.
    * The first numSkip values in the array are skipped (not used) to make the regression.
    * This is useful if we want to focus the regression on larger values of n. 
    */
   public double[] regressionLogVariance (int numSkip) {
		double[] x2 = new double[numSets-numSkip], y2 = new double[numSets-numSkip];
		for (int i = 0; i < numSets-numSkip; ++i) {
			x2[i] = logn[i+numSkip];
			y2[i] = logVar[i+numSkip];
		}
		return LeastSquares.calcCoefficients(x2, y2, 1);
	}
   
   /**
    * Takes the regression coefficients of log(variance) in #regCoeff and returns a two-line string 
    * that reports on these coefficients.  
    * @param  regCoeff  the regression coefficients.
    * @return  Report as a string.
    */
    public String formatRegression (double[] regCoeff) {
		StringBuffer sb = new StringBuffer("");
		// double[] regCoeff = regressionLogVariance (numSkipRegression);
		sb.append("  Slope of log(var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		return sb.toString();
	}
	
   /**
    * Produces and returns a report on the last experiment.
    * @param numSkip  The first numSkip values of n are skipped for the regression
    * @param details  If true, gives values (mean, log variance,...) for each n.
    * @return  Report as a string.
    */
	public String reportVarianceRate (int numSkip, boolean details) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for mean estimation: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" RQMC point sets: " + theSets[0].toString() + "\n\n");
		sb.append("RQMC variance \n");
		if (details) sb.append(dataLogForPlot());
		sb.append (formatRegression (regressionLogVariance (numSkip)));
		// sb.append("  Slope of log(var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		// sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		sb.append("  Total CPU Time = " + cpuTime + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
	}
	
	/**
	 * Takes the data from the most recent experiment and returns it in a @ref PgfDataTable.
	 * This will typically be used to plot the data.
	 * 
	 * @param tableName  Name (short identifier) of the table.
	 * @return Report as a string.
	 */
	public PgfDataTable toPgfDataTable(String tableName, String tableLabel) {
        double[][] data = new double[numSets][5];
		for (int s = 0; s < numSets; s++) { // For each cardinality n
			data[s][0] = size[s];
		    data[s][1] = mean[s];
		    data[s][2] = variance[s];
		    data[s][3] = logn[s];
		    data[s][4] = logVar[s];
		}
		return new PgfDataTable (tableName, tableLabel, tableFields, data);
	}

	public PgfDataTable toPgfDataTable(String tableLabel) {
		return toPgfDataTable (title, tableLabel);
	}

	

	/**
	 * Returns the data on the mean and variance for each n, in an appropriate format to produce a
	 * plot with the pgfplot package.
	 * 
	 * @return Report as a string.
	 */
	public String dataForPlot() {
		StringBuffer sb = new StringBuffer("");
		sb.append("    n      mean       variance \n");
		for (int s = 0; s < numSets; s++)  // For each cardinality n
			sb.append(" " + size[s] + " " + PrintfFormat.f(10, 5, mean[s]) + " "
			        + PrintfFormat.f(10, 5, variance[s]) + "\n");
		return sb.toString();
	}

	/**
	 * Similar to dataForPlot, but for the log(variance) in terms of log n.
	 * 
	 * @return Report as a string.
	 */
	public String dataLogForPlot() {
		StringBuffer sb = new StringBuffer("");
		sb.append("   log(n)      mean       log(variance) \n");
		for (int s = 0; s < numSets; s++)  // For each cardinality n
			sb.append(" " + logn[s] + " " + PrintfFormat.f(10, 5, mean[s]) + " "
			        + PrintfFormat.f(10, 5, logVar[s]) + "\n");
		return sb.toString();
	}
	
	/*
	 * Returns the data on the mean and variance for each n, in an appropriate format to produce a
	 * plot with the pgfplot package.   This is OBSOLETE.
	 * 
	 * @return Report as a string.

	public String XformatPgfCurve (String curveName) {
		StringBuffer sb = new StringBuffer("");
		sb.append("      \\addplot+[no marks] table[x=n,y=variance] {" + "\n");
		sb.append( dataForPlot() + " } \n");
		sb.append("      \\addlegendentry{" + curveName + "}\n"); 
		sb.append("      % \n");
		return sb.toString();
	}
	 */
	

	/**
	 * Performs an experiment (testVarianceRate) for each point set series in the given list,
	 * and returns a report as a string. 
	 * 
	 * @param model
	 * @param list
	 * @param m
	 * @return  a report on the experiment.
	 */
	public String testVarianceRateManyPointTypes (MonteCarloModelDouble model, 
			ArrayList<RQMCPointSet[]> list,
			int m, int numSkip, 
			boolean makePgfTable, boolean printReport, boolean details,
			ArrayList<PgfDataTable> listCurves) {
		StringBuffer sb = new StringBuffer("");
	    // if (makePgfTable)  
	    //	listCurves = new ArrayList<PgfDataTable>();
		for(RQMCPointSet[] ptSeries : list) {
			init (ptSeries, base);
			testVarianceRate (model, m);
			if (printReport)  System.out.println(reportVarianceRate (numSkip, details));			
			if (printReport) sb.append (reportVarianceRate (numSkip, details));			
            if (makePgfTable == true)  listCurves.add (toPgfDataTable 
            		(ptSeries[0].getLabel()));
		}
		return sb.toString();
	}

	/**
	 * Performs an experiment (testVarianceRate) for each point set series in the given list,
	 * and returns a report as a string. 
	 */
	public String toString () {
		return ("RQMC Experiment:" + title);   // theSets[0].toString();
	}
}