package ift6561examples;

import umontreal.ssj.util.*;
import umontreal.ssj.stat.*;
import java.io.*;

public class CallCenterCRN extends CallCenter {
   Tally statDiffIndep = new Tally("stats on difference with IRNs");
   Tally statDiffCRN = new Tally("stats on difference with CRNs");
   int[] numAgents1, numAgents2;

   public CallCenterCRN(String fileName) throws IOException {
      super(fileName);
      numAgents1 = new int[numPeriods];
      numAgents2 = new int[numPeriods];
      for (int j = 0; j < numPeriods; j++)
         numAgents1[j] = numAgents2[j] = numAgents[j];
   }

   // Set the number of agents in each period j to the values in num.
   public void setNumAgents(int[] num) {
      for (int j = 0; j < numPeriods; j++)
         numAgents[j] = num[j];
   }

   public void simulateDiffCRN(int n) {
      double value1, value2;
      statDiffIndep.init();
      statDiffCRN.init();
      for (int i = 0; i < n; i++) {
         setNumAgents(numAgents1);
         streamB.resetNextSubstream();
         streamArr.resetNextSubstream();
         streamPatience.resetNextSubstream();
         (genServ.getStream()).resetNextSubstream();
         simulateOneDay();
         value1 = (double) nGoodQoS / nCallsExpected;
         setNumAgents(numAgents2);
         streamB.resetStartSubstream();
         streamArr.resetStartSubstream();
         streamPatience.resetStartSubstream();
         (genServ.getStream()).resetStartSubstream();
         simulateOneDay();
         value2 = (double) nGoodQoS / nCallsExpected;
         statDiffCRN.add(value2 - value1);
         simulateOneDay();
         value2 = (double) nGoodQoS / nCallsExpected;
         statDiffIndep.add(value2 - value1);
      }
   }

   static public void main(String[] args) throws IOException {
      int n = 1000; // Number of replications.
      CallCenterCRN cc = new CallCenterCRN(
            args.length == 1 ? args[0] : "src/main/docs/examples/ift6561examples/CallCenter.dat");
      cc.numAgents2[5]++;
      cc.numAgents2[6]++;
      cc.simulateDiffCRN(n);
      System.out.println(cc.statDiffIndep.reportAndCIStudent(0.9, 5) + cc.statDiffCRN.reportAndCIStudent(0.9, 5));
      double varianceIndep = cc.statDiffIndep.variance();
      double varianceCRN = cc.statDiffCRN.variance();
      System.out.println("Variance ratio:  " + PrintfFormat.format(10, 1, 3, varianceIndep / varianceCRN));
   }
}
