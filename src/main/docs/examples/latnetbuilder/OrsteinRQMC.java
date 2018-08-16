import umontreal.ssj.rng.*;
import umontreal.ssj.hups.*;
import umontreal.ssj.latnetbuilder.DigitalNetSearch;
import umontreal.ssj.latnetbuilder.PolynomialLatticeSearch;
import umontreal.ssj.stat.Tally;

public class OrsteinRQMC extends Orstein {

	public OrsteinRQMC (double sigma, double theta, int s) {
       super (sigma, theta, s);
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
      int s = 10;
      double sigma = 1;
      double theta = 8;
      OrsteinRQMC process = new OrsteinRQMC (sigma, theta, s);
      
      Tally statValue  = new Tally ("value of Asian option");
      int n = 100000;
      System.out.println ("Ordinary MC:\n");
      process.simulateRuns (n, new MRG32k3a(), statValue);
      statValue.setConfidenceIntervalStudent();
      System.out.println (statValue.report (0.95, 3));
      double varMC = statValue.variance();
      System.out.println ("Total Variance: " + varMC);
      
      int m = 1000; 

      for (int k=5; k<20; k++) {
    	  System.out.println ("------------------------\n");
    	  System.out.println("k: " + k);
    	  Tally statQMC = new Tally ("QMC averages for Asian option");
    	  DigitalNet p = new SobolSequence (k, 31, s); // 2^{k} points. 
          
          n = p.getNumPoints();
          process.simulateQMC (m, p, new MRG32k3a(), statQMC);
          System.out.println ("QMC with standard Sobol point set with " + n +
              " points and affine matrix scramble:\n");
          statQMC.setConfidenceIntervalStudent();
          System.out.println (statQMC.report (0.95, 3));
          double varQMC = p.getNumPoints() * statQMC.variance();
          System.out.printf ("Variance ratio:   %9.4g%n", varMC/varQMC);
    	  
    	  ///////////////////////////////////////////////
          Tally statQMC2 = new Tally ("QMC averages for Asian option");
    	  System.out.println("Searching for best polynomial lattice rule...");
    	  DigitalNetSearch search = new PolynomialLatticeSearch("lattice");
    	  earch.setPathToLatNetBuilder("/home/anaconda3/envs/latnetbuilder/bin/latnetbuilder");
    	  search.setDimension(s);
    	  search.setSizeParameter("2^" + k);
    	  search.setFigureOfMerit("CU:P2");
    	  search.setNormType("2");
    	  search.addWeight("product:0.3:1,0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.45");
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
      }
   }
}
