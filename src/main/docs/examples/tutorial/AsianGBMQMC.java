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

   // Makes m independent randomizations of the digital net p using stream
   // noise. For each of them, performs one simulation run for each point
   // of p, and adds the average over these points to the collector statQMC.
   public void simulateRQMC (int m, RQMCPointSet prqmc, Tally statRQMC) {
      Tally statValue  = new Tally ("stat on value of Asian option");
      PointSetIterator stream = prqmc.iterator ();
      for (int j=0; j<m; j++) {
    	  prqmc.randomize();
          stream.resetStartStream();
          simulateRuns (prqmc.getNumPoints(), stream, statValue);
          statRQMC.add (statValue.average());
      }
   }


   public static void main (String[] args)  throws IOException {
      int s = 12;
      double[] zeta = new double[s+1];
      for (int j=0; j<=s; j++)
         zeta[j] = (double)j / (double)s;
      AsianGBMQMC process = new AsianGBMQMC (0.05, 0.5, 100.0, 100.0, s, zeta);
      Tally statValue  = new Tally ("value of Asian option");
      Tally statRQMC = new Tally ("RQMC averages for Asian option under GBM");

      Chrono timer = new Chrono();
      int n = 100000;
      System.out.println ("Ordinary MC:\n");
      process.simulateRuns (n, new MRG32k3a(), statValue);
      statValue.setConfidenceIntervalStudent();
      System.out.println (statValue.report (0.95, 3));
      System.out.println ("Total CPU time: " + timer.format());
      double varMC = statValue.variance();
      double cpuMC = timer.getSeconds() / n;  // CPU seconds per run.
      System.out.println ("------------------------\n");

      timer.init();
      DigitalNet p = new SobolSequence (16, 31, s); // 2^{16} points.
      PointSetRandomization rand = new LMScrambleShift (new MRG32k3a());
      RQMCPointSet prqmc = new RQMCPointSet (p, rand);
      n = p.getNumPoints();           // Number of RQMC points.
      int m = 20;                     // Number of RQMC randomizations.
      process.simulateRQMC (m, prqmc, statRQMC);
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
