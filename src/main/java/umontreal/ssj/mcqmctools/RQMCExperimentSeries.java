/*
 * Class:        RQMCPointSetSeries
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
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.list.lincv.ListOfTalliesWithCV;
import umontreal.ssj.util.Chrono;
import umontreal.ssj.util.Num;
import umontreal.ssj.util.PrintfFormat;
import java.util.ArrayList;


/**
 * This class offers facilities to perform experiments to study the convergence
 * of the variance when estimating a mean (expectation) with a series of RQMC 
 * point sets usually of the same type, but different size @f$n@f.
 * The series of RQMC point sets of different sizes can be passed in an array 
 * to the constructor.   The method @f$testVarianceRate@f$ performs an experiment 
 * with a given model and the series of point sets.  One can recover the average,
 * variance, their logs in base 2, etc., in arrays, as well as the estimated 
 * linear regression of log(variance) as a function of log(n). 
 */

public class RQMCExperimentSeries {
	int numSets = 0;  // Number of point sets in the series.
    RQMCPointSet[] theSets = null;   
	double[] size = new double[numSets];    // values of n
	double[] mean = new double[numSets];    // averge performance for each point set 
	double[] variance = new double[numSets]; // variance for each point set
	double[] log2n = new double[numSets];   // log_2 n
	double[] log2Var = new double[numSets]; // log_2 of variance
	boolean displayExec = false;   // When true, prints a display of execution in real time
	int numReplicates;    // last value of m
	MonteCarloModelDouble model;
	int numSkipRegression = 0; // Number of values of n that are skipped for the regression
	String cpuTime;       // time for last experiment


   /**
    * Constructor with a give series of RQMC point sets.
    *  @param theSets      the RQMC point sets
    */
   public RQMCExperimentSeries (RQMCPointSet[] theSets) {
	   init(theSets);
   }

   /**
    * When set to true, a real-time display of execution results and CPU times 
    * will be printed on the default output.
    */
   public void setExecutionDisplay (boolean display) {
	      this.displayExec = display;
   }
   
   /**
    * Resets the array of RQMC point sets for this object, and initializes 
    * (re-creates) the arrays that will contain the results.
    */
   public void init(RQMCPointSet[] theSets) {
	     this.numSets = theSets.length;
	     this.theSets = theSets;
		 size = new double[numSets]; //  n for each point set
	     mean = new double[numSets]; //  average for each point set
	     variance = new double[numSets]; // variance for each point set
	     log2n = new double[numSets];    // log_2 n
	     log2Var = new double[numSets];  // log_2 of the variance
   }
   
   /**
    * Returns the point set number i associated to this object (starts at 0).
    *  @return the ith point set associated to this object
    */
   public RQMCPointSet getSet(int i) {
      return theSets[i];
   }

   /**
    * Performs an RQMC experiment with the given model, with this series of RQMC point sets.  
    * For each set in the series, computes the average, the variance, its log in base 2.
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
			System.out.println("    n     CPU time         mean       log2(var) ");	    	
	    }
		for (int s = 0; s < numSets; s++) { // For each cardinality n
			n = theSets[s].getNumPoints();
			size[s] = n;
			log2n[s] = Num.log2(n);
			// System.out.println(" n = " + n + ", Lg n = " + log2n[s] + "\n"); // ****
			// System.out.println("  " + n + "     " + timer.format());
			RQMCExperiment.simulReplicatesRQMC (model, theSets[s], m, statReps);
			mean[s] = statReps.average();
			variance[s] = statReps.variance();
		    log2Var[s] = Num.log2(variance[s]);
		    if (displayExec) {
			   System.out.println("  " + n + "     " + timer.format() + 
			              "   " + PrintfFormat.f(10, 5, mean[s]) + 
					      "   " + PrintfFormat.f(7, 2, log2Var[s]));
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
			System.out.println("    n     CPU time         mean       log2(var) ");	    	
	    }
		for (int s = 0; s < numSets; s++) { // For each cardinality n
			n = theSets[s].getNumPoints();
			size[s] = n;
			log2n[s] = Num.log2(n);
			// System.out.println(" n = " + n + ", Lg n = " + log2n[s] + "\n"); // ****
			// System.out.println("  " + n + "     " + timer.format());
			RQMCExperiment.simulReplicatesRQMCCV (model, theSets[s], m, statWithCV);
			statWithCV.estimateBeta();    // This is where the var. and covar. are computed!
			mean[s] = statWithCV.averageWithCV(0);
			variance[s] = statWithCV.covarianceWithCV(0, 0);
		    log2Var[s] = Num.log2(variance[s]);
		    if (displayExec) {
			   System.out.println("  " + n + "     " + timer.format() + 
			              "   " + PrintfFormat.f(10, 5, mean[s]) + 
					      "   " + PrintfFormat.f(7, 2, log2Var[s]));
		    }
		}	   
        cpuTime = timer.format();	   
   }

   /**
    * Sets the number of (small) values of n to skip for regression.
    */
   public void setNumSkipRegression (int numSkip) {
      numSkipRegression = numSkip;
   }

   /**
    * Sets the number of (small) values of n to skip for regression.
    */
   public void setDisplayExec (boolean displayExec) {
      this.displayExec = displayExec;
   }

   /**
    * Returns the vector of means from last experiment.
    */
   public double[] getMeans() {
      return mean;
   }

   /**
    * Returns the vector of variances from last experiment.
    */
   public double[] getVariances() {
      return variance;
   }

   /**
    * Returns the vector of log_2 of variances from last experiment.
    */
   public double[] getLog2Variances() {
      return log2Var;
   }

   /**
    * Returns the vector of values of n.
    */
   public double[] getValuesn() {
      return size;
   }

   /**
    * Returns the vector of log_2(n).
    */
   public double[] getLog2n() {
      return log2n;
   }

   /**
    * Performs a linear regression of log_2(variance) vs log_2(n), and returns the 
    * coefficients (constant and slope) in two-dimensional vector.
    * The first numSkip values in the array are skipped (not used) to make the regression.
    * This is useful if we want to focus the regression on larger values of n. 
    */
   public double[] regressionLogVariance (int numSkip) {
		double[] x2 = new double[numSets-numSkip], y2 = new double[numSets-numSkip];
		for (int i = 0; i < numSets-numSkip; ++i) {
			x2[i] = log2n[i+numSkip];
			y2[i] = log2Var[i+numSkip];
		}
		return LeastSquares.calcCoefficients(x2, y2, 1);
	}
   
   /**
    * Produces and returns a report on the last experiment.
    * @param numSkip  The first numSkip values of n are skipped for the regression
    * @param details  If true, gives values (mean, log variance,...) for each n.
    * @return  Report as a string.
    */
	public String report (boolean details) {
		StringBuffer sb = new StringBuffer("");
		sb.append("\n ============================================= \n");
		sb.append("RQMC simulation for mean estimation: \n ");
		sb.append("Model: " + model.toString() + "\n");
		sb.append(" Number of indep copies m  = " + numReplicates + "\n");
		sb.append(" Point sets: " + this.toString() + "\n\n");
		sb.append("RQMC variance \n");
		if (details) {
			sb.append("    n      mean       log2(var) \n");
			for (int s = 0; s < numSets; s++) { // For each cardinality n
				sb.append(" " + size[s] + " " + PrintfFormat.f(10, 5, mean[s]) +
				          " " + PrintfFormat.f(7, 2, log2Var[s]) + "\n");
			}
		}
		double[] regCoeff = regressionLogVariance (numSkipRegression);
		sb.append("  Slope of log2(var) = " + PrintfFormat.f(8, 5, regCoeff[1]) + "\n");
		sb.append("    constant term      = " + PrintfFormat.f(8, 5, regCoeff[0]) + "\n\n");
		sb.append("  Total CPU Time = " + cpuTime + "\n");
		sb.append("-----------------------------------------------------\n");		
		return sb.toString();
	}
	
	/**
	 * Performs an experiment (testVarianceRate) for each point set series in the given list,
	 * and returns a report as a string. 
	 * 
	 * @param model
	 * @param list
	 * @param m
	 */
	public String TestRQMCManyPointTypes (MonteCarloModelDouble model, 
			ArrayList<RQMCPointSet[]> list, int m, boolean details) {
		StringBuffer sb = new StringBuffer("");
		numReplicates = m;	
		for(RQMCPointSet[] ptSeries : list) {
			init (ptSeries);
			testVarianceRate (model, m);
			sb.append (report (details));			
		}
		return sb.toString();
	}
	
	/**
	 * Performs an experiment (testVarianceRate) for each point set series in the given list,
	 * and returns a report as a string. 
	 */
	public String toString () {
		return theSets[0].toString();
	}
}