package ift6561examples;

import java.io.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.mcqmctools.*;
import umontreal.ssj.stat.*;

/**
 * This class is for the specific stochastic activity network in class San13.
 * The goal is to estimate the probability that the length of the longest path
 * exceed a given constant x.
 */

public class San13Prob extends San13 {

   double x; // Threshold for longest path: we estimate P[T > x].

   // The constructor reads link length distributions in a file.
   public San13Prob(double x, String fileName) throws IOException {
      super(fileName);
      this.x = x;
   }

   // Returns the indicator that T > x.
   public double getPerformance() {
      if (maxPath > x)
         return 1.0;
      else
         return 0.0;
   }

   public void setx(double x) {
      this.x = x;
   }

   public String toString() {
      String s = "SAN network with 9 nodes and 13 links, from Elmaghraby (1977)\n" + "Estimate P[T > x} for x = " + x
            + "\n";
      return s;
   }

   public static void main(String[] args) throws IOException {
      int n = 100000;
      double x = 90.0;
      San13Prob san = new San13Prob(x, "src/main/docs/examples/ift6561examples/san13a.dat");
      Tally statT = new Tally("SAN13Prob example");

      // Run a simulation experiment with n runs.
      System.out.println(MonteCarloExperiment.simulateRunsDefaultReportStudent(san, n, new LFSR113(), statT, 0.95, 5));
   }
}
