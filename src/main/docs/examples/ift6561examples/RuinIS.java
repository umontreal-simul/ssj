package ift6561examples;
// Rare event simulation with IS: ruin probability with exponential
// claim sizes (rate beta) and exponential interarrival times
// (rate lambda).  

import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.util.*;

public class RuinIS {

   static final double lambda     = 1.0;    // Arrival rate of claims.
   static final double beta       = 0.5;    // 1 / expected claim size.
   static final double r0         = 200.0;  // Initial reserve.
   double theta;

   RandomVariateGen genArrivals;   // For claim arrivals.
   RandomVariateGen genSizes;      // For claim sizes.
   Tally statIs = new Tally ("Ruin probability with IS");

   public double simulRuin (double c) {
      double sum = 0.0;
      while (sum < r0)
         sum += genSizes.nextDouble() - c * genArrivals.nextDouble();
      return Math.exp (- theta * sum);
   }

   public RuinIS (double c, int n) {
      // Computes IS parameters and makes n simulation runs with IS,
      // with input rate c for the premiums.
      theta = (c - 2.0) / (2.0 * c);
      double lambdaIs = lambda + c * theta; 
      double betaIs = beta - theta;
      genArrivals = new RandomVariateGen
         (new MRG32k3a(), new ExponentialDist (lambdaIs));
      genSizes = new RandomVariateGen
         (new MRG32k3a(), new ExponentialDist (betaIs));
      for (int i=0; i < n; i++) {
         statIs.add (simulRuin (c));
         }
      System.out.println (" lambda = " + 
         PrintfFormat.format (8, 3, 1, lambda));
      System.out.println (" beta   = " + 
         PrintfFormat.format (8, 3, 1, beta));
      System.out.println (" c      = " + 
         PrintfFormat.format (8, 3, 1, c));
      System.out.println (" R(0)   = " + 
         PrintfFormat.format (8, 3, 1, r0));
      System.out.println (" n      = " + n);
      System.out.println ();
      System.out.println (statIs.formatCIStudent (0.90));
      System.out.println (" Variance with IS = " + 
         PrintfFormat.format (10, 2, 2, statIs.variance()));
      double p = statIs.average();
      System.out.println (" Sample size for 10% error with MC = " + 
         PrintfFormat.format (10, 2, 2, 100.0 * (1.0-p) / p));
      System.out.println (" Sample size for 10% error with IS = " + 
         PrintfFormat.format (10, 2, 2, 100.0 * statIs.variance() / (p*p)));
      System.out.println ();
      System.out.println ("---------------------------------------------");
   }

   public static void main (String[] args) { 
      new RuinIS (3.0, 10000);
      new RuinIS (5.0, 10000);
      new RuinIS (10.0, 10000);
   }
}
