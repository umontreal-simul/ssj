package ift6561examples;

import java.io.*;
//import java.util.Scanner;
import umontreal.ssj.rng.*;
import umontreal.ssj.stat.*;
//import umontreal.ssj.probdist.*;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;


/**
 * @author Pierre L'Ecuyer
 * 
 */
public class San13CMCCV extends San13Prob {

    double estimate; // Cond. prod. of exceeding x.
	int[] indexV = new int[8];  // The 8 control variates.
	double[] meansCV = new double[8];  // CV expectations.
	TallyStore[] statsCV = new TallyStore [8];  // Stats on control variates.
  	TallyStore statsProb = new TallyStore ("Original CMC estimator");

    public San13CMCCV (double x, String fileName) throws IOException {
		super(x, fileName);
		// The CVs: 0, 1, 2, 3, 7, 10, 11, 12
		indexV[0] = 0; indexV[1] = 1; indexV[2] = 2; indexV[3] = 3;
		indexV[4] = 7; indexV[5] = 10; indexV[6] = 11; indexV[7] = 12;
		for (int k = 0; k < 8; k++)
		{
			statsCV[k] = new TallyStore();
			meansCV[k] = dist[indexV[k]].getMean();
		}
	}

	public void simulateRuns (int n, RandomStream stream) {
		for (int k = 0; k < 8; k++)
			statsCV[k].init();
		statsProb.init();
		for (int i = 0; i < n; i++)
		{
			simulate(stream);
			// stream.resetNextSubstream();
		}
	}


	public void simulate(RandomStream stream)
	{
		int k;
		// Generate: 0, 1, 2, 3, 7, 10, 11, 12
		for (int j = 0; j < 8; j++)
		{
			k = indexV[j];
			V[k] = dist[k].inverseF(stream.nextDouble());
			statsCV[j].add(V[k]);
			if (V[k] < 0.0)
				V[k] = 0.0;
		}
		double prod = 1.0;
		// Incomplete path lengths
		paths[0] = V[1] + V[10];
		paths[1] = V[0] + V[2] + V[10];
		if (paths[0] > paths[1])
			paths[1] = paths[0];
		prod *= dist[5].cdf(x - paths[1]);
		paths[2] = V[0] + V[10];
		prod *= dist[4].cdf(x - paths[2]);
		paths[3] = V[0] + V[3] + V[7] + V[10];
		prod *= dist[9].cdf(x - paths[3]);
		paths[4] = V[0] + V[3] + V[7] + V[12];
		prod *= dist[8].cdf(x - paths[4]);
		paths[5] = V[0] + V[3] + V[11] + V[12];
		prod *= dist[6].cdf(x - paths[5]);
		estimate = 1.0 - prod;
		statsProb.add(estimate);
	}


	public String toString() {
		String s = "SAN network with 9 nodes and 13 links, from Elmaghraby (1977)\n"
				+ "Estimate prob longest path > x = " + x + ", using CMC + control variates.\n";
		return s;
	}

   	static Algebra alg = new Algebra();

   	/**
   	* Applies a vector of control variables to an estimator.
    	* The tally x contains the observations used to estimate the output
    	* average while the array of tallies c contains the
    	* observations of the control variates.
    	* The third array, ec, contains the known expectations
    	* of the control variates.
    	* The length c and ec should be equal while x, and c should
    	* contain the same number of observations.
    	*
    	* This method estimates the optimal vector of constants for
    	* the control variates, and prints information on the vector,
    	* and the variance with CV.     
    	* @param x the tally containing the observations for the output value.
    	* @param c the tally containing the observations of the control variates.
    	* @param ec the vector of expectations.
    	*/
	
	public static void applyCV (TallyStore x, TallyStore[] c, double[] ec) {
      // Construct and fill the matrix Cov[C]
      DoubleMatrix2D matC = new DenseDoubleMatrix2D (c.length, c.length);
      for (int i = 0; i < c.length; i++)
         matC.setQuick (i, i, c[i].variance ());
      for (int i = 0; i < c.length - 1; i++)
         for (int j = i + 1; j < c.length; j++) {
            double cov = c[i].covariance (c[j]);
            matC.setQuick (i, j, cov);
            matC.setQuick (j, i, cov);
         }
      // Construct and fill the vector Cov[C, X]
      DoubleMatrix2D matCX = new DenseDoubleMatrix2D (c.length, 1);
      for (int i = 0; i < c.length; i++)
         matCX.setQuick (i, 0, x.covariance (c[i]));

      DoubleMatrix2D mbeta;
      try {
         // Find Beta vector solving Cov[C, X] = Cov[C]*Beta
         mbeta = alg.solve (matC, matCX);
      }
      catch (IllegalArgumentException iae) {
         // This can happen, e.g., if the variance of a CV is (incorrectly) 0.
         System.out.println ("Cannot apply CV");
         System.out.println ();
         return;
      }
      
      // Compute average Xc = X - Beta^t * (C - E[C])
      double avgWithCV = x.average ();
      for (int i = 0; i < c.length; i++)
         avgWithCV -= mbeta.getQuick (i, 0) * (c[i].average () - ec[i]);
      // Compute variance Var[Xc] = Var[X] + Beta^t*Var[C]*Beta - 2Beta*Cov[C, X]
      double varWithCV = x.variance ();
      // viewDice transposes the matrix mbeta (which contains a single column),
      // and zMult performs the matrix multiplication.
      // The null second argument instructs Colt to create a new matrix for the result.
      // The result of the operation is a 1x1 matrix from which we extract the single element;
      // this is the second term of the controlled variance.
      varWithCV += mbeta.viewDice ().zMult (matC, null).zMult (mbeta, null).getQuick (0, 0);
      // A similar technique is used to compute the third term of the controlled variance.
      varWithCV -= 2*mbeta.viewDice ().zMult (matCX, null).getQuick (0, 0);

      // Print the results
      System.out.print ("Beta vector with CV                        : (");
      for (int i = 0; i < c.length; i++)
         System.out.printf ("%s%.3g", i > 0 ? ", " : "", mbeta.getQuick (i, 0));
      System.out.println (")");
      System.out.printf ("Average without CV                         : %8.5g%n", x.average ());
      System.out.printf ("Average with CV                            : %8.5g%n", avgWithCV);
      System.out.printf ("Variance without CV                        : %8.5g%n", x.variance ());
      System.out.printf ("Variance with CV                           : %8.5g%n", varWithCV);
      System.out.printf ("Variance reduction factor                  : %8.5g%n", x.variance () / varWithCV);

	  int n = x.numberObs();
	  double delta = Math.sqrt(x.variance() / n);
	  double LB_sansCV = x.average() - 1.96 * delta;
      double UB_sansCV = x.average() + 1.96 * delta;
	  delta = Math.sqrt(varWithCV / n);
	  double LB_avecCV = avgWithCV - 1.96 * delta;
      double UB_avecCV = avgWithCV + 1.96 * delta;
	  System.out.printf("IC 95 pourcent sans CV: (%8.5g, %8.5g)%n", LB_sansCV, UB_sansCV);
	  System.out.printf("IC 95 pourcent avec CV: (%8.5g, %8.5g)%n", LB_avecCV, UB_avecCV);
	  // System.out.println("IC 95% avec CV: (" + LB_avecCV + ", " + UB_avecCV + ")");     
	  System.out.println();
   }

	public static void main(String[] args) throws IOException {
		int n = 1000 * 100;
		San13CMCCV san = new San13CMCCV (90.0, "san13a.dat");
		san.simulateRuns(n, new LFSR113());
		System.out.println (san.statsProb.reportAndCIStudent(0.95, 6));
		applyCV (san.statsProb, san.statsCV, san.meansCV);

	}
}
