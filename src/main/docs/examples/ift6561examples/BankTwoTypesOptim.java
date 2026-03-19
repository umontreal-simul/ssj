package ift6561examples;

import java.io.*;

import umontreal.ssj.stat.Tally;

/*
 * For the BankTwoTypes model, this program estimates the sensitivity (or derivative)
 * of three performance measures with respect to the parameters mu_a and s,
 * using finite differences with common random numbers (CRNs).
 * This program is similar to CallCenterCRN, but we are not doing IRNs.
 */
public class BankTwoTypesOptim extends BankTwoTypesCRN {
   
   // Constructor: reads data in file and creates generators.
   public BankTwoTypesOptim (String fileName) throws IOException {
      super(fileName);
   }
  
   // Estimate the expected cost as a function of s, with CRNs.
   public void simulateGridCRNs (int n, int nums, double s0, double delta) {
      double meanA, meanB, cost;
      Tally statCost = new Tally("Cost per day");
      System.out.println (" s        meanA             meanB              total cost \n");
      for (int j = 0; j < nums; j++) {
         s = s0 + j * delta;
         myStreams.resetStartStream();
         statWaitsA.init();
         statWaitsB.init();
         statCost.init();
         for (int i = 0; i < n; i++) {
            myStreams.resetNextSubstream();
            simulateOneDay();      
            statCost.add(statWaitsInDayA.sum() + 5.0 * statWaitsInDayB.sum());
         }
         meanA = statWaitsA.average() * nExpectedA;
         meanB = statWaitsB.average() * nExpectedB;
         cost = meanA + 5.0 * meanB;
         //  cost = statCost.average();
         System.out.println (s + "  " + meanA + "  " + meanB + "  " + cost);
         // System.out.println (statCost.reportAndCIStudent(0.95, 1));
      }
   }
   
   // Compares the expected cost for three values of s, with CRNs.
   public void simulateDiffCRNs (int n, double s1, double s2, double s3) {
      Tally statDiffCost12 = new Tally("Diff in cost per day for s1 vs s2");
      Tally statDiffCost23 = new Tally("Diff in cost per day for s2 vs s3");
      double cost1, cost2, cost3;
      for (int i = 0; i < n; i++) {
         s = s1;
         myStreams.resetNextSubstream();
         simulateOneDay();
         cost1 = statWaitsInDayA.sum() + 5.0 * statWaitsInDayB.sum();
         s = s2;
         myStreams.resetStartSubstream();
         simulateOneDay();
         cost2 = statWaitsInDayA.sum() + 5.0 * statWaitsInDayB.sum();
         statDiffCost12.add (cost2 - cost1);        
         s = s3;
         myStreams.resetStartSubstream();
         simulateOneDay();
         cost3 = statWaitsInDayA.sum() + 5.0 * statWaitsInDayB.sum();
         statDiffCost23.add (cost3 - cost2); 
         }
      System.out.println (statDiffCost12.reportAndCIStudent(0.95, 1));
      System.out.println (statDiffCost23.reportAndCIStudent(0.95, 1));
   }

   static public void main(String[] args) throws IOException {
      BankTwoTypesOptim cc = new BankTwoTypesOptim(
         args.length == 1 ? args[0] : "src/main/docs/examples/ift6561examples/BankTwoTypes.dat");
      int n = 10000;        // Number of simulation runs.
      int nums = 8;         // Number of values of s. 
      double s0 = 60.0;     // Smallest values of s.
      double delta = 30.0;  // Spacings between values of s.
      cc.simulateGridCRNs (n, nums, s0, delta);

      n = 1000000;   // Number of simulation runs.
      nums = 7;      // Number of values of s. 
      s0 = 120.0;    // Smallest values of s.
      delta = 10.0;  // Spacings between values of s.
      cc.simulateGridCRNs (n, nums, s0, delta);
      
      n = 1000000;
      cc.simulateDiffCRNs (n, 140.0, 150.0, 160.0);
      cc.simulateDiffCRNs (n, 135.0, 140.0, 145.0);
      }
}
