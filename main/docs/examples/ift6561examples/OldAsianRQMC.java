package ift6561examples;

import umontreal.ssj.rng.*;
import umontreal.ssj.hups.*;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stochprocess.*;

public class OldAsianRQMC extends AsianOption {
	
   public OldAsianRQMC (double r, int d, double[] obsTimes, double strike)  {
       super (r, d, obsTimes, strike);
   }

   // Makes m independent randomizations of the digital net p using stream
   // noise. For each of them, performs one simulation run for each point
   // of p, and adds the average over these points to the collector statQMC.
   
   public void simulateQMC (int m, PointSet p, PointSetRandomization rand,
		   RandomStream noise, Tally statQMC) {
	   Tally statValue  = new Tally ("stat on value of Asian option");
	   PointSetIterator stream = p.iterator ();
	   for (int j=0; j<m; j++) {	
		   rand.randomize(p);
		   stream.resetStartStream();
		   simulateRuns (p.getNumPoints(), stream, statValue);
		   statQMC.add (statValue.average());
		   }
	   }
   

   public static void main (String[] args) {
	   int d = 16;
	   double[] obsTimes = new double[d + 1];
	   obsTimes[0] = 0.0;
	   for (int j = 1; j <= d; j++)
		   obsTimes[j] = (double) j / (double) d;		
	   int n = 1000;
	   RandomStream stream = new MRG32k3a();
	   double mu = -0.1436; 			
	   double sigma = 0.12136; 		
	   double theta = mu;
	   double r = 0.1; 				
	   double s0 = 100.0; 				
	   double v = 0.3;	
	   int m = 32;                     // Number of QMC randomizations.
	   OldAsianRQMC asian = new OldAsianRQMC(0.1, d, obsTimes, 101.0);


	  StochasticProcess process = new GeometricVarianceGammaProcess(s0, r,
		   new VarianceGammaProcessDiff( 0,  theta,  sigma, v, 
		   new GammaProcessSymmetricalBridge(0,r,v,stream), 
		   new GammaProcessSymmetricalBridge(0,r,v,stream)));

	  asian.setProcess(process);
	  Tally statValue = new Tally("Stats on value of Asian option");		
	  Tally statQMC = new Tally ("QMC averages for Asian option");

      System.out.println ("Ordinary MC:\n");
      asian.simulateRuns (n, stream, statValue);
      statValue.setConfidenceIntervalStudent();
      System.out.println (statValue.report (0.95, 3));
      double varMC = statValue.variance();
      System.out.println ("------------------------\n");

     // SobolSequence p = new SobolSequence (14, 31, 32); // 2^{14} points.
      //LMScrambleShift rand = new LMScrambleShift(new MRG32k3a());
      RandomShift rand = new RandomShift(new MRG32k3a());


     // KorobovLattice p = new KorobovLattice (16381, 5693, 32);
     // BakerTransformedPointSet pb = new BakerTransformedPointSet (p);
     
      int a[] = {1, 6229, 2691, 3349, 5893, 7643, 7921, 7055, 4829, 5177, 5459, 4863, 4901, 2833, 2385, 3729, 981, 957, 4047, 1013, 1635, 2327, 7879, 2805, 2353, 1081, 3999, 879, 5337, 7725, 4889, 5103};
      Rank1Lattice p = new Rank1Lattice((int) Math.pow(2, 14), a, 32);
		
      n = p.getNumPoints();
		
		//Point set with y = 2/(2+j)
		
      asian.simulateQMC (m, p, rand, stream, statQMC);
      
      System.out.println ("QMC with Sobol point set with " + n +
          " points and affine matrix scramble:\n");
      statQMC.setConfidenceIntervalStudent();
      System.out.println (statQMC.report (0.95, 3));
      double varQMC = p.getNumPoints() * statQMC.variance();
      System.out.printf ("Variance ratio:   %9.4g%n", varMC/varQMC);
   }
}
