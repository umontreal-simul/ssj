package tutorial;

import java.io.IOException;
import umontreal.ssj.rng.*;
import umontreal.ssj.hups.*;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.Chrono;

// An extension of AsianGBM that uses RQMC point sets.
public class AsianGBMQMC extends AsianGBM {

   public AsianGBMQMC (double r, double sigma, double strike,
                    double s0, int s, double[] zeta) {
       super (r, sigma, strike, s0, s, zeta);
   }

   // Makes m independent randomizations of the RQMC point set prqmc.
   // For each of them, performs one simulation run for each point
   // of prqmc, and adds the average over these points to the collector statRQMC.
   public void simulateRunsRQMC (int m, RQMCPointSet prqmc, Tally statRQMC) {
      Tally statValue  = new Tally ("stat on value of Asian option");
      int n = prqmc.getNumPoints();
      PointSetIterator stream = prqmc.iterator ();
      for (int j=0; j<m; j++) {
    	  prqmc.randomize();
          stream.resetStartStream();
          simulateRuns (n, stream, statValue);
          statRQMC.add (statValue.average());
      }
   }


   public static void main (String[] args)  throws IOException {
      int d = 12;
      double[] zeta = new double[d+1];
      for (int j=0; j<=d; j++)
         zeta[j] = (double)j / (double)d;
      AsianGBMQMC process = new AsianGBMQMC (0.05, 0.5, 100.0, 100.0, d, zeta);
      Tally statMC  = new Tally ("value of Asian option");
      Tally statRQMC = new Tally ("RQMC averages for Asian option under GBM");
      Chrono timer = new Chrono();

      // We first perform a Monte Carlo experiment, to compare with RQMC.
      int n = 100000;
      System.out.println ("Ordinary MC:\n");
      process.simulateRuns (n, new MRG32k3a(), statMC);
      statMC.setConfidenceIntervalStudent();
      System.out.println (statMC.report (0.95, 3));
      System.out.println ("Total CPU time: " + timer.format());
      double varMC = statMC.variance();
      double cpuMC = timer.getSeconds() / n;  // CPU seconds per run.
      System.out.println ("------------------------\n");

      // Then we make a RQMC experiment, and compare the work-normalized variances.
      timer.init();
      DigitalNet p = new SobolSequence (16, 31, d); // n = 2^{16} points in d dim.
      PointSetRandomization rand = new LMScrambleShift (new MRG32k3a());
      RQMCPointSet prqmc = new RQMCPointSet (p, rand);
      n = p.getNumPoints();           // Number of RQMC points.
      int m = 20;                     // Number of RQMC randomizations.
      process.simulateRunsRQMC (m, prqmc, statRQMC);
      System.out.println ("QMC with Sobol point set with " + n +
          " points and affine matrix scramble:\n");
      statRQMC.setConfidenceIntervalStudent();
      System.out.println (statRQMC.report (0.95, 3));
      System.out.println ("Total CPU time: " + timer.format() + "\n");
      double varQMC = p.getNumPoints() * statRQMC.variance();
      double cpuQMC = timer.getSeconds() / (m * n);
      System.out.printf ("Variance ratio:   %9.4g%n", varMC/varQMC);
      System.out.printf ("Efficiency ratio: %9.4g%n",
           (varMC * cpuMC) / (varQMC * cpuQMC));
   }
}
