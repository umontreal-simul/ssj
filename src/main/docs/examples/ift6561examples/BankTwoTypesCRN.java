package ift6561examples;

import umontreal.ssj.randvar.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.stat.*;
import java.io.*;

/*
 * For the BankTwoTypes model, this program estimates the sensitivity (or derivative)
 * of three performance measures with respect to the parameters mu_a and s,
 * using finite differences with common random numbers (CRNs).
 * This program is similar to CallCenterCRN, but we are not doing IRNs.
 */
public class BankTwoTypesCRN extends BankTwoTypes {

   Tally statDiffAD  = new Tally("Diff of AD matches");
   Tally statDiffWaitsA = new Tally("Diff of waiting times for type A");
   Tally statDiffWaitsB = new Tally("Diff of waiting times for type B");
   // The following will be two lognormal generators that use the same stream,
   // but with different means for the lognormal (mu_a and mu_a + delta).
   LognormalGen genServA1, genServA2; 
   
   // Constructor: reads data in file and creates generators.
   public BankTwoTypesCRN (String fileName) throws IOException {
      super(fileName);
   }


   // Estimate the sensitivity w.r.t. mu_a, with CRNs.
   public void simulateDiffCRNMua (int n, double delta) {
      double meanAD1, meanAD2, meanA1, meanA2, meanB1, meanB2;
      double mua = genServA.getDistribution().getMean();
      double vara = genServA.getDistribution().getVariance();
      genServA1 = genServA;
      genServA2 = new LognormalGen (genServA.getStream(), 
            new LognormalDistFromMoments (mua + delta, vara));
      statDiffAD.init();
      statDiffWaitsA.init();
      statDiffWaitsB.init();
         
      for (int i = 0; i < n; i++) {
         genServA = genServA1;
         myStreams.resetNextSubstream();
         simulateOneDay();
         meanAD1 = (double) numMatchesAD;
         meanA1 = (double) statWaitsInDayA.sum() / nExpectedA;
         meanB1 = (double) statWaitsInDayB.sum() / nExpectedB;
         
         genServA = genServA2;
         myStreams.resetStartSubstream();
         simulateOneDay();
         meanAD2 = (double) numMatchesAD;
         meanA2 = (double) statWaitsInDayA.sum() / nExpectedA;
         meanB2 = (double) statWaitsInDayB.sum() / nExpectedB;
 
         statDiffAD.add(meanAD2 - meanAD1);
         statDiffWaitsA.add(meanA2 - meanA1);
         statDiffWaitsB.add(meanB2 - meanB1);
      }
   }

  
   // Estimate the sensitivity w.r.t. the threshold s, with CRNs.
   public void simulateDiffCRNs (int n, double delta) {
      double meanAD1, meanAD2, meanA1, meanA2, meanB1, meanB2;
      double s1, s2;
      s1 = s;  s2 = s + delta;
      statDiffAD.init();
      statDiffWaitsA.init();
      statDiffWaitsB.init();
      
      for (int i = 0; i < n; i++) {
         s = s1;
         myStreams.resetNextSubstream();
         simulateOneDay();
         meanAD1 = (double) numMatchesAD;
         meanA1 = (double) statWaitsInDayA.sum() / nExpectedA;
         meanB1 = (double) statWaitsInDayB.sum() / nExpectedB;
         
         s = s2;
         myStreams.resetStartSubstream();
         simulateOneDay();
         meanAD2 = (double) numMatchesAD;
         meanA2 = (double) statWaitsInDayA.sum() / nExpectedA;
         meanB2 = (double) statWaitsInDayB.sum() / nExpectedB;
 
         statDiffAD.add(meanAD2 - meanAD1);
         statDiffWaitsA.add(meanA2 - meanA1);
         statDiffWaitsB.add(meanB2 - meanB1);
      }
   }

  

   static public void main(String[] args) throws IOException {
      
      int n = 10000;  // Number of simulation runs.
      double delta;
      BankTwoTypesCRN cc = new BankTwoTypesCRN(
         args.length == 1 ? args[0] : "src/main/docs/examples/ift6561examples/BankTwoTypes.dat");

      delta = 5.0;
      cc.simulateDiffCRNMua (n, delta);
      System.out.println ("==================================================\n");
      System.out.println ("Sensitivity with respect to mu_a, with delta = " + delta + "\n");
      System.out.println (cc.statDiffAD.reportAndCIStudent(0.95, 4));
      System.out.println (cc.statDiffWaitsA.reportAndCIStudent(0.95, 5));
      System.out.println (cc.statDiffWaitsB.reportAndCIStudent(0.95, 5));

      delta = 60.0;
      cc.simulateDiffCRNs (n, delta);
      System.out.println ("==================================================\n");
      System.out.println ("Sensitivity with respect to s, with delta = " + delta + "\n");
      System.out.println (cc.statDiffAD.reportAndCIStudent(0.95, 4));
      System.out.println (cc.statDiffWaitsA.reportAndCIStudent(0.95, 4));
      System.out.println (cc.statDiffWaitsB.reportAndCIStudent(0.95, 4));
      System.out.println ("==================================================\n");
      }
}
