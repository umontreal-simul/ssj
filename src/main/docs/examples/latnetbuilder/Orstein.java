import org.apache.commons.math3.linear.*;

import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.Chrono;

public class Orstein {
   double sigma;
   double theta;
   int s;
   RealVector vals;
   RealMatrix A;
   
   public Orstein (double sigma, double theta, int s) {
      this.sigma = sigma;
      this.theta = theta;
      this.s = s;
      
		Array2DRowRealMatrix mat = new Array2DRowRealMatrix(s, s);
		
		for (int i=0; i<s; i++) {
			for (int j=0; j<s; j++) {
				mat.setEntry(i, j, sigma*sigma / (2 * theta) * (- 1 + Math.exp(2 * theta * (Math.min(i, j) +1) / s )) * Math.exp(- theta * (i+j+2) / s) );
			}
		}
		
		EigenDecomposition decomp = new EigenDecomposition(mat);
		RealMatrix D = decomp.getD();
		RealMatrix P = decomp.getV();
		for (int i=0; i<s; i++) {
			D.setEntry(i, i, Math.sqrt(D.getEntry(i, i)));
		}
		this.A = P.multiply(D);
   }

   // Generates the process S.
   public void generatePath (RandomStream stream) {
	   ArrayRealVector v = new ArrayRealVector(s);
       for (int j = 0; j < s; j++) {
    	   double x = stream.nextDouble();
    	   v.setEntry(j, NormalDist.inverseF01 (x));
       }
       vals = A.operate(v);
   }

   public double getErrorOnEnergy () {
       double energy = 0.0;
       for (int j = 0; j < s; j++) {
    	  // energy += Math.pow(vals.getEntry(j), 2);
    	  energy += Math.abs(vals.getEntry(j));
       }
       // double thEnergy = sigma*sigma/(2 * theta) * (1 - (1- Math.exp(-2*theta)) / (s * (Math.exp(2*theta/s) - 1)) );
       // return energy/s - thEnergy;
       return energy/s;
   }

   // Performs n indep. runs using stream and collects statistics in statValue.
   public void simulateRuns (int n, RandomStream stream, Tally statValue) {
      statValue.init();
      for (int i=0; i<n; i++) {
         generatePath (stream);
         statValue.add (getErrorOnEnergy ());
         stream.resetNextSubstream();
      }
   }

   public static void main (String[] args) {
      int s = 5;
      double sigma = 1;
      double theta = 32;
      Orstein process = new Orstein (sigma, theta, s);
      Tally statValue = new Tally ("Stats on value of Sobol test func");

      Chrono timer = new Chrono();
      int n = 1000000;
      process.simulateRuns (n, new MRG32k3a(), statValue);
      statValue.setConfidenceIntervalStudent();
      System.out.println (statValue.report (0.95, 3));
      System.out.println ("Total CPU time:      " + timer.format() + "\n");
   }

}
