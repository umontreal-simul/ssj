package ift6561examples;

import java.io.*;
// import java.util.Scanner;
import umontreal.ssj.rng.*;
// import umontreal.ssj.probdist.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.mcqmctools.*;

/**
 * This class simulates the same static reliability network as in
 * NetworkReliablity13, but does it using CMC. The first two links Y[0] and Y[1]
 * are not generated. This is for Example intro-san-reliability-CMC in the book.
 */

public class NetworkReliability13CMC1 extends NetworkReliablity13 {

   double probConnected;

   // The constructor.
   public NetworkReliability13CMC1(double rj) {
      super(rj);
   }

   public void simulate(RandomStream stream) {
      for (int k = 0; k < 13; k++) {
         Y[k] = (stream.nextDouble() < r[k]);
      }
      boolean I1, I2;
      probConnected = 0.0;
      // Connection indicators for nodes 1 and 2.
      I1 = ((Y[2] & Y[5] & Y[10]) | (Y[4] & Y[10]) | (Y[3] & Y[7] & Y[9] & Y[10]) | (Y[3] & Y[7] & Y[8] & Y[12])
            | (Y[3] & Y[6] & Y[11] & Y[12]));
      I2 = (Y[5] & Y[10]);
      if (I1)
         probConnected = r[0];
      if (I2)
         probConnected += r[1] * (1.0 - probConnected);
   }

   public double getPerformance() {
      return (1.0 - probConnected);
   }

   public String toString() {
      String s = "SAN network with 9 nodes and 13 links, from Elmaghraby (1977)\n"
            + "Estimate prob that nodes 0 and 8 are disconnected";
      return s;
   }

   public static void main(String[] args) throws IOException {
      int n = 1 * 1000;
      double rj = 0.999;
      NetworkReliability13CMC1 net = new NetworkReliability13CMC1(rj);
      Tally statR = new Tally("Network reliability");
      net.r[0] = net.r[1] = 0.95; // For the CMC example in the book.
      System.out.println(MonteCarloExperiment.simulateRunsDefaultReportStudent(net, n, new LFSR113(), statR, 0.95, 4));
      System.out.println(statR.report(0.95, 8));
   }
}
