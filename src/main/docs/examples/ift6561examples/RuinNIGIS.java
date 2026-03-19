package ift6561examples;
// Rare event simulation with IS: ruin probability with 
// NIG (alpha, beta, mu, delta) claim sizes
// and exponential interarrival times (rate lambda).  

import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.functions.*;
import umontreal.ssj.util.*;

public class RuinNIGIS {

   static final double lambda = 1.0; // Arrival rate of claims.
   static final double alpha2 = 5.0; // alpha square.
   static final double alpha = Math.sqrt(alpha2);
   static final double beta = 2.0;   // NIG parameters.
   static final double gamma = Math.sqrt(alpha2 - beta * beta);
   static final double mu = 1.0;
   static final double delta = 1.0;
   double c;
   double theta;

   RandomVariateGen genArrivals; // For claim arrivals.
   RandomVariateGen genSizes;    // For claim sizes.
   Tally statIs = new Tally("Ruin probability with IS");
   Tally statClaims = new Tally("Stats on claim sizes");

   public RuinNIGIS(double c) {
      // Computes IS parameters and makes n simulation runs with IS,
      // with input rate c for the premiums.
      this.c = c;
      theta = RootFinder.brentDekker (0.001, alpha-beta, new RootMGF(), 1.0E-6);   
      
      double lambdaIs = lambda + c * theta;
      genArrivals = new RandomVariateGen(new MRG32k3a(), new ExponentialDist(lambdaIs));
      double betaIs = beta + theta;
      double gammaIs = Math.sqrt (alpha2 - betaIs * betaIs);   
      genSizes = new NormalInverseGaussianIGGen(
          new InverseGaussianGen (new MRG32k3a(), delta/gammaIs, delta*delta), 
          new NormalGen (new MRG32k3a(), 0.0, 1.0), betaIs, mu);
     
      System.out.println(" lambda = " + PrintfFormat.format(8, 5, 1, lambda));
      System.out.println(" alpha2 = " + PrintfFormat.format(8, 5, 1, alpha2));
      System.out.println(" beta   = " + PrintfFormat.format(8, 5, 1, beta));
      System.out.println(" gamma  = " + PrintfFormat.format(8, 5, 1, gamma));
      System.out.println(" mu     = " + PrintfFormat.format(8, 5, 1, mu));
      System.out.println(" delta  = " + PrintfFormat.format(8, 5, 1, delta));
      System.out.println(" theta  = " + PrintfFormat.format(8, 5, 1, theta));
      System.out.println(" lambdaIs = " + PrintfFormat.format(8, 5, 1, lambdaIs));
      System.out.println(" betaIs = " + PrintfFormat.format(8, 5, 1, betaIs));
      System.out.println(" gammaIs= " + PrintfFormat.format(8, 5, 1, gammaIs));
      System.out.println(" c      = " + PrintfFormat.format(8, 5, 1, c));
      System.out.println();
      
      double mean = mu + delta * beta / gamma;
      double var = delta * alpha2 / (gamma * gamma * gamma);
      System.out.println(" mean claim = " + PrintfFormat.format(10, 5, 1, mean));
      System.out.println(" var claim  = " + PrintfFormat.format(10, 5, 1, var));
      mean = mu + delta * betaIs / gammaIs;
      var = delta * alpha2 / (gammaIs * gammaIs * gammaIs);
      System.out.println(" mean claim IS = " + PrintfFormat.format(10, 5, 1, mean));
      System.out.println(" var claim  IS = " + PrintfFormat.format(10, 5, 1, var));
      
      statClaims.init();
      for (int i = 0; i < 10000; i++) {
         statClaims.add(genSizes.nextDouble());
      }
      System.out.println (statClaims.reportAndCIStudent(0.95));
   }

   
   // Defines the function for which we find a root to solve for M(theta)=1. 
   public class RootMGF implements MathFunction {
      // double c;     // Function parameter c.
      
      // public RootMGF (double c) { this.c = c; }
      
      public double evaluate (double theta) {
         return mu * theta + delta * gamma - Math.sqrt(alpha2 - (beta+theta)*(beta+theta)) 
               - Math.log (1.0 + c*theta);
      }
   }

   public double simulate (double r0) {
      double sum = 0.0;
      while (sum < r0)
         sum += genSizes.nextDouble() - c * genArrivals.nextDouble();
      return Math.exp(-theta * sum);
   }

   public void simulateRuns (double r0, int n) {
      statIs.init();
      for (int i = 0; i < n; i++) {
         statIs.add(simulate (r0));
      }
      System.out.println(" Simulation results: ");
      System.out.println(" R(0)   = " + PrintfFormat.format(8, 3, 1, r0));
      System.out.println(" n      = " + n);
      System.out.println(statIs.formatCIStudent(0.90));
      System.out.println(" Variance with IS = " + PrintfFormat.format(10, 2, 2, statIs.variance()));
      double p = statIs.average();
      System.out
            .println(" Sample size for 10% error with MC = " + PrintfFormat.format(10, 2, 2, 100.0 * (1.0 - p) / p));
      System.out.println(" Sample size for 10% error with IS = "
            + PrintfFormat.format(10, 2, 2, 100.0 * statIs.variance() / (p * p)));
      System.out.println();
      System.out.println("---------------------------------------------");
   }

    
   public static void main(String[] args) {
      double r0 = 200.0;
      new RuinNIGIS(3.1).simulateRuns(r0, 1000);
      new RuinNIGIS(3.3).simulateRuns(r0, 1000);
      new RuinNIGIS(3.6).simulateRuns(r0, 1000);
      new RuinNIGIS(4.0).simulateRuns(r0, 1000);
      new RuinNIGIS(5.0).simulateRuns(r0, 1000);
      
      r0 = 20.0;
      new RuinNIGIS(3.1).simulateRuns(r0, 1000);
      new RuinNIGIS(3.3).simulateRuns(r0, 1000);
      new RuinNIGIS(3.6).simulateRuns(r0, 1000);
      new RuinNIGIS(4.0).simulateRuns(r0, 1000);
      new RuinNIGIS(5.0).simulateRuns(r0, 1000);
   }
}
