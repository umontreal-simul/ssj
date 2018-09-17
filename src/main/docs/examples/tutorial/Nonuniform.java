package tutorial;

import umontreal.ssj.rng.*;
import java.io.*;
import umontreal.ssj.charts.HistogramChart;
import umontreal.ssj.probdist.*;
import umontreal.ssj.randvar.*;
import umontreal.ssj.stat.*;

public class Nonuniform {
   // The parameter values are hardcoded here to simplify the program.
   double lambda = 5.0;   double p = 0.2;
   double alpha = 2.0;    double beta = 1.0;
   double mu = 5.0;       double sigma = 0.5;

   RandomStream stream = new LFSR113();
   RandomVariateGenInt genN = new RandomVariateGenInt
	      (stream, new PoissonDist (lambda));       // For N
   RandomVariateGen genY = new GammaAcceptanceRejectionGen
	      (stream, new GammaDist (alpha, beta));    // For Y_j
   RandomVariateGen genW = new RandomVariateGen
	      (stream, new LognormalDist (mu, sigma));  // For W_j

   // Generates and returns X.
   public double simulate () {
      int N;  int M;  int j;  double X = 0.0;
      N = genN.nextInt();
      M = GeometricDist.inverseF (p, stream.nextDouble());  // Uses static method
      for (j = 0; j < N; j++) X += genY.nextDouble();
      for (j = 0; j < M; j++) X += genW.nextDouble();
      return X;
   }

   // Performs n indep. runs and collects statistics in statX.
   public void simulateRuns (int n, TallyStore statX) {
	   for (int i=0; i<n; i++) statX.add (simulate ());
    }

   public static void main (String[] args) throws IOException {
	  int n = 100000;
	  TallyStore statX = new TallyStore (n); // To store the n observations of X.
      (new Nonuniform ()).simulateRuns (n, statX);  // Simulate X n times.
      System.out.println (statX.report (0.95, 1));

      // Compute and print the empirical quantiles.
      statX.quickSort();
      double[] data = statX.getArray();  // The sorted observations.
      System.out.printf (" 0.10 quantile: %9.1f%n", data[(int)(0.10 * n)]);
      System.out.printf (" 0.50 quantile: %9.1f%n", data[(int)(0.50 * n)]);
      System.out.printf (" 0.90 quantile: %9.1f%n", data[(int)(0.90 * n)]);
      System.out.printf (" 0.99 quantile: %9.1f%n", data[(int)(0.99 * n)]);
      
      // Make a histogram of the empirical distribution of X.
      HistogramChart hist = new HistogramChart("Histogram of distribution of $X$",
				"Values of $X$", "Frequency", statX.getArray(), n);
  	  double[] bounds = { 0, 4000, 0, 25000 }; // Range of x and y to be displayed.
	  hist.setManualRange(bounds);
	  (hist.getSeriesCollection()).setBins(0, 40, 0, 4000); // 40 bins over [0, 4000].
	  hist.view(800, 500);  // View on screen.
	  
	  // Make a Latex file that contains the histogram.
	  String histLatex = hist.toLatex(12.0, 8.0);  // Width and height of plot in cm.
	  Writer file = new FileWriter("src/main/docs/examples/tutorial/NonuniformHist.tex"); 
	  file.write(histLatex);
	  file.close();
   }
}
