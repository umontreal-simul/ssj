package tutorial;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.randvar.*;
import umontreal.ssj.stat.*;

public class Nonuniform {
   // The parameter values are hardwired here to simplify the program.
   double lambda = 5.0;   double p = 0.2;
   double alpha = 2.0;    double beta = 1.0;
   double mu = 5.0;       double sigma = 1.0;

   RandomStream stream = new LFSR113();
   RandomVariateGenInt genN = new RandomVariateGenInt
	      (stream, new PoissonDist (lambda));       // For N
   RandomVariateGen genY = new GammaAcceptanceRejectionGen
	      (stream, new GammaDist (alpha, beta));    // For Y_j
   RandomVariateGen genW = new RandomVariateGen
	      (stream, new LognormalDist (mu, sigma));  // For W_j

   // Generates and returns X.
   public double generateX () {
      int N;  int M;  int j;  double X = 0.0;
      N = genN.nextInt();
      M = GeometricDist.inverseF (p, stream.nextDouble());  // Uses static method
      for (j = 0; j < N; j++) X += genY.nextDouble();
      for (j = 0; j < M; j++) X += genW.nextDouble();
      return X;
   }

   // Performs n indep. runs and collects statistics in statX.
   public void simulateRuns (int n) {
      TallyStore statX = new TallyStore (n);
      for (int i=0; i<n; i++) statX.add (generateX ());
      System.out.println (statX.report ());
      statX.quickSort();
      double[] data = statX.getArray();
      System.out.printf ("0.10 quantile: %9.3f%n", data[(int)(0.10 * n)]);
      System.out.printf ("0.50 quantile: %9.3f%n", data[(int)(0.50 * n)]);
      System.out.printf ("0.90 quantile: %9.3f%n", data[(int)(0.90 * n)]);
      System.out.printf ("0.99 quantile: %9.3f%n", data[(int)(0.99 * n)]);
   }

   public static void main (String[] args) {
      (new Nonuniform ()).simulateRuns (10000);
   }
}
