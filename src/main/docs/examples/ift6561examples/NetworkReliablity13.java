package ift6561examples;

import java.io.*;
// import java.util.Scanner;
import umontreal.ssj.rng.*;
// import umontreal.ssj.probdist.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.mcqmctools.*;

/**
 * This class simulates a static reliability model using the same network as for
 * the specific San13 example with 9 nodes and 13 links, taken from Elmaghraby
 * (1977). The goal here is to estimate the probability that nodes 1 and 9 are
 * disconnected. This program is very specific to this example and uses a very
 * naive way to compute the shortest path, by enumerating all six paths! Here
 * the links are numbered from 0 to 12.
 */

public class NetworkReliablity13 implements MonteCarloModelDouble {

   double[] r = new double[13];
   boolean[] Y = new boolean[13];
   // We consider the 6 paths that can lead to the sink.
   boolean connected;

   // The constructor.
   public NetworkReliablity13(double rj) {
      for (int k = 0; k < 13; k++)
         r[k] = rj;
   }

   public void simulate(RandomStream stream) {
      for (int k = 0; k < 13; k++) {
         Y[k] = (stream.nextDouble() < r[k]);
      }
      connected = false;
      // Check all directed connection paths
      if ((Y[1] & Y[5] & Y[10]) | (Y[0] & Y[2] & Y[5] & Y[10]) | (Y[0] & Y[4] & Y[10])
            | (Y[0] & Y[3] & Y[7] & Y[9] & Y[10]) | (Y[0] & Y[3] & Y[7] & Y[8] & Y[12])
            | (Y[0] & Y[3] & Y[6] & Y[11] & Y[12]))
         connected = true;
   }

   public double getPerformance() {
      if (connected)
         return 0.0;
      else
         return 1.0;
   }

   public String toString() {
      String s = "SAN network with 9 nodes and 13 links, from Elmaghraby (1977)\n"
            + "Estimate prob that nodes 1 and 9 are disconnected";
      return s;
   }
   
   public String getTag() {
      return "SAN-Reliability13";
   }

   public static void main(String[] args) throws IOException {
      int n = 1 * 1000;
      double rj = 0.999;
      NetworkReliablity13 net = new NetworkReliablity13(rj);
      Tally statC = new Tally("Network reliability");
      net.r[0] = net.r[1] = 0.95; // This is for the CMC vs MC example.
      MonteCarloExperiment.simulateRunsDefaultReportStudent(net, n, new LFSR113(), statC, 0.95, 4);
      System.out.println(statC.report(0.95, 8));
   }
}
