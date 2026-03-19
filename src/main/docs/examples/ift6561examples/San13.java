package ift6561examples;

import java.io.*;
import java.util.Scanner;
import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.probdist.ContinuousDistribution;
import umontreal.ssj.probdist.DistributionFactory;
import umontreal.ssj.rng.*;

/**
 * This class simulates a specific stochastic activity network with 9 nodes and
 * 13 links, taken from Elmaghraby (1977) and used again in L'Ecuyer and Lemieux
 * (2000), "Variance Reduction via Lattice Rules". It simulates and returns the
 * length T of the longest path in the network. This program is very specific to
 * this example and uses a very naive way to compute the longest path, by
 * enumerating all six paths! Here the arcs are numbered from 0 to 12.
 **/

public class San13 implements MonteCarloModelDouble {

   int dim = 13; // model dimension (number of uniforms needed for a simulation).
   double[] V = new double[dim];
   ContinuousDistribution[] dist = new ContinuousDistribution[dim];
   // We consider the 6 paths that can lead to the sink.
   double[] paths = new double[6]; // Path lengths.
   double maxPath; // Length of the current longest path.

   // The constructor reads link length distributions in a file.
   public San13(String fileName) throws IOException {
      readDistributions(fileName);
   }

   public void readDistributions(String fileName) throws IOException {
      // Reads data and construct arrays.
      BufferedReader input = new BufferedReader(new FileReader(fileName));
      Scanner scan = new Scanner(input);
      for (int k = 0; k < 13; k++) {
         dist[k] = DistributionFactory.getContinuousDistribution(scan.nextLine());
      }
      scan.close();
   }

   public int getDimension() {
      return dim;
   }

   // Returns the length of longest path when the random lengths
   // are replaced by their means.
   public double deterministicLengths() {
      for (int k = 0; k < 13; k++) {
         V[k] = dist[k].getMean();
         if (V[k] < 0.0)
            V[k] = 0.0;
      }
      return computePathsAndT();
   }

   public void simulate(RandomStream stream) {
      for (int k = 0; k < 13; k++) {
         V[k] = dist[k].inverseF(stream.nextDouble());
         if (V[k] < 0.0)
            V[k] = 0.0;
      }
      computePathsAndT();
   }

   // Compute the lengths of all paths and returns the longest length T
   public double computePathsAndT() {
      // Path lengths
      paths[0] = V[1] + V[5] + V[10];
      paths[1] = V[0] + V[2] + V[5] + V[10];
      paths[2] = V[0] + V[4] + V[10];
      paths[3] = V[0] + V[3] + V[7] + V[9] + V[10];
      paths[4] = V[0] + V[3] + V[7] + V[8] + V[12];
      paths[5] = V[0] + V[3] + V[6] + V[11] + V[12];
      maxPath = paths[0];
      for (int p = 1; p < 6; p++)
         if (paths[p] > maxPath)
            maxPath = paths[p];
      return maxPath;
   }

   // Returns the length T of longest path.
   public double getPerformance() {
      return maxPath;
   }

   public String toString() {
      String s = "SAN network with 9 nodes and 13 links, from Elmaghraby (1977)\n"
            + "Simulates and returns the length T of the longest path \n";
      return s;
   }
   
   public String getTag() {
      return "SAN13";
   }

}
