package latnetbuilder;

import umontreal.ssj.rng.*;
import umontreal.ssj.hups.*;
import umontreal.ssj.latnetbuilder.DigitalNetSearch;
import umontreal.ssj.latnetbuilder.PolynomialLatticeSearch;
import umontreal.ssj.stat.Tally;

import java.lang.Math;

public class SobolTestFuncRQMC extends SobolTestFunc {

   public SobolTestFuncRQMC (double c, int s) {
       super (c, s);
   }

   // Makes m independent randomizations of the digital net p using stream
   // noise. For each of them, performs one simulation run for each point
   // of p, and adds the average over these points to the collector statQMC.
   public void simulateQMC (int m, DigitalNet p,
                            RandomStream noise, Tally statQMC) {
      Tally statValue  = new Tally ("stat on value of Asian option");
      PointSetIterator stream = p.iterator ();
      for (int j=0; j<m; j++) {
         p.leftMatrixScramble (noise);
         p.addRandomShift (0, p.getDimension(), noise);
          stream.resetStartStream();
          simulateRuns (p.getNumPoints(), stream, statValue);
          statQMC.add (statValue.average());
      }
   }


   public static void main (String[] args) {
      for (int s = 5; s<= 65; s +=10) {
    	  for (double c = 0.1; c <=1.01; c += 0.1) {
    		  int k = 12;
    		  System.out.println("-----------------------------------");
    		  System.out.println("s: " + s);
    		  System.out.println("c: " + c);
    		  System.out.println("k: " + k);
    		  SobolTestFuncRQMC process = new SobolTestFuncRQMC (c, s);
    		  Tally statValue  = new Tally ("value of Sobol Test function");
    		  Tally statQMC = new Tally ("QMC averages for Sobol Test function");
    		  Tally statQMC2 = new Tally ("QMC averages for Sobol Test function");
    		  Tally statQMC3 = new Tally ("QMC averages for Sobol Test function");
    		  
    		  int n = 10000;
    		  System.out.println ("Ordinary MC:\n");
    		  process.simulateRuns (n, new MRG32k3a(), statValue);
    		  statValue.setConfidenceIntervalStudent();
    		  System.out.println (statValue.report (0.95, 3));
    		  double varMC = statValue.variance();
    		  System.out.println ("------------------------\n");
    		  
    		  ////////////////////////////////////////////////////////////////
    		  DigitalNet p = new SobolSequence (k, 31, s); // 2^{k} points.   
    		  
    		  n = p.getNumPoints();
    		  int m = 1000;                     // Number of QMC randomizations.
    		  process.simulateQMC (m, p, new MRG32k3a(), statQMC);
    		  System.out.println ("QMC with standard Sobol point set with " + n +
    				  " points and affine matrix scramble:\n");
    		  statQMC.setConfidenceIntervalStudent();
    		  System.out.println (statQMC.report (0.95, 3));
    		  double varQMC = p.getNumPoints() * statQMC.variance();
    		  System.out.printf ("Variance ratio:   %9.4g%n", varMC/varQMC);
    		  
    		  ////////////////////////////////////////////////////////////////
    		  System.out.println("Searching for best polynomial lattice rule...");
    		  DigitalNetSearch search = new PolynomialLatticeSearch("lattice");
    		  search.setPathToLatNetBuilder("/home/anaconda3/envs/latnetbuilder/bin/latnetbuilder");
    		  search.setDimension(s);
    		  search.setSizeParameter("2^" + k);
    		  search.setFigureOfMerit("CU:P2");
    		  search.setNormType("2");
    		  double w = c * c;
    		  search.addWeight("product:" + w);
    		  search.setExplorationMethod("fast-CBC");
    		  DigitalNetBase2 q = search.search();
    		  System.out.println("Found best polynomial lattice rule...");
    		  
    		  n = q.getNumPoints();
    		  process.simulateQMC (m, q, new MRG32k3a(), statQMC2);
    		  System.out.println ("QMC with custom polynomial point set with " + n +
    				  " points and affine matrix scramble:\n");
    		  statQMC2.setConfidenceIntervalStudent();
    		  System.out.println (statQMC2.report (0.95, 3));
    		  varQMC = q.getNumPoints() * statQMC2.variance();
    		  System.out.printf ("Variance ratio:   %9.4g%n", varMC/varQMC);
    		  
    		  ////////////////////////////////////////////////////////////////
			  System.out.println("Searching for best sobol net...");
			  search = new DigitalNetSearch("sobol");
			  search.setPathToLatNetBuilder("/home/pmarion/Documents/Stage_3A/latsoft/bin/latnetbuilder");
			  search.setDimension(s);
			  search.setSizeParameter("2^" + k);
			  search.setFigureOfMerit("t-value");
			  search.addWeight("order-dependent:0:" + w + "," + Math.pow(w, 2) + "," + Math.pow(w, 3));
			  search.setNormType("inf");
			  search.setExplorationMethod("random:1000");
			  q = search.search();
			  System.out.println("Found best sobol net...");
			  
			  n = q.getNumPoints();
			  process.simulateQMC (m, q, new MRG32k3a(), statQMC3);
			  System.out.println ("QMC with random sobol net with " + n +
					  " points and affine matrix scramble: \n");
			  statQMC3.setConfidenceIntervalStudent();
			  System.out.println (statQMC3.report (0.95, 3));
			  varQMC = q.getNumPoints() * statQMC3.variance();
			  System.out.printf ("Variance ratio:   %9.4g%n", varMC/varQMC);
    	  }
      
    	 }
   }
}
