package ift6561examples;
import umontreal.ssj.stochprocess.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.hups.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.util.Chrono;

public class TestAsianVGRQMC2017 {

	 
	   public static void main (String[] args) {
	        
	       int d = 8;
		   double T1 = 1.0 / d;
		   double T = 1.0;
		   
		   double mu = -0.1; 			
		   double sigma = 0.15; 		
		   double theta = mu;
		   double r = 0.05; 				
		   double nu = 0.2;	
		   double K = 52.0;
		   double s0 = 50.0;

		   int n = 100000;  // 2^14  for Monte Carlo
		   int m = 100;         // Number of RQMC replications.

		   AsianOption asian = new AsianOption (r, d, T1, T, K);
		   RandomStream noise = new MRG32k3a();
		   Chrono timer = new Chrono();
		   Tally statValue = new Tally("Stats on value of Asian option for MC");		
		   Tally statRQMC = new Tally ("RQMC averages for Asian option");


/*  Note that VarianceGammaProcessAlternate is not yet in SSJ. 
 *  We need a general class that sample any process with subordinator (with random clock)
 *  in an alternate way.  Can be implemented as a wrapper.
 * */		   
		  StochasticProcess spBGSS = new GeometricVarianceGammaProcess(s0, r,
			  new VarianceGammaProcessAlternate(0.0,
			     new BrownianMotion(0.0, mu, sigma, noise), 
			     new GammaProcess(0.0, 1.0, nu, noise)));

		  StochasticProcess spBGBS = new GeometricVarianceGammaProcess(s0, r,
			  new VarianceGammaProcessAlternate(0.0,
				   new BrownianMotionBridge(0.0, mu, sigma, noise), 
				   new GammaProcessSymmetricalBridge(0.0, 1.0, nu, noise)));

		  double dummy = 1.0;  // Dummy parameters used for creating the gamma processes.
		  StochasticProcess spDGBS = new GeometricVarianceGammaProcess(s0, r,
			  new VarianceGammaProcessDiff(0.0,  theta,  sigma, nu, 
				   new GammaProcessSymmetricalBridge(0.0, dummy, dummy, noise), 
				   new GammaProcessSymmetricalBridge(0.0, dummy, dummy, noise)));
		  
	      RandomShift randShift = new RandomShift(new MRG32k3a());
	      LMScrambleShift randLMS = new LMScrambleShift(new MRG32k3a());

	      SobolSequence pSobol = new SobolSequence (14, 31, 2*d); // 2^{14} points.

	      KorobovLattice pKor = new KorobovLattice (16381, 5693, 2*d);
	      BakerTransformedPointSet pKorBaker = new BakerTransformedPointSet (pKor);
	     
		  // Lattice point set found with gamma_j = 2/(2+j).
	      int a1[] = {1, 6229, 2691, 3349, 5893, 7643, 7921, 7055, 4829, 5177, 5459, 4863, 
	                  4901, 2833, 2385, 3729, 981, 957, 4047, 1013, 1635, 2327, 7879, 2805, 
	                  2353, 1081, 3999, 879, 5337, 7725, 4889, 5103};
	      Rank1Lattice pLat1 = new Rank1Lattice(16 * 1024, a1, 2*d);
	      BakerTransformedPointSet pLat1Baker = new BakerTransformedPointSet (pLat1);

	      System.out.println ("Pricing an Asian Option under a VG process. ");
	      // System.out.println ("The exact mean and MC variance are approximately: ");
	      // System.out.println ("Mean = 5.725 and variance = 29.89. \n");	      

	      // Monte Carlo experiments first, with BGSS.
	      n = 1000000;   // for MC.
	      System.out.println ("Ordinary MC:\n");	      
		  asian.setProcess(spBGSS);
	      MonteCarloExperiment.simulateRunsDefaultReport(asian, n, noise, statValue, timer); 
	      double varMC = statValue.variance();
	      double secondsMC = timer.getSeconds() / n;

		  // asian.setProcess(spBGBS);
	      // MonteCarloExperiment.simulateRunsDefaultReport(asian, n, noise, statValue, timer); 

		  // asian.setProcess(spDGBS);
	      // MonteCarloExperiment.simulateRunsDefaultReport(asian, n, noise, statValue, timer); 

	      System.out.println ("-----------------------------------------------------\n");
	      
/*	      // Brownian gamma sequential sampling
	      System.out.println ("\n *****  Brownian gamma sequential sampling (BGSS)   *****\n");
		  asian.setProcess(spBGSS);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pSobol, randShift, 
	    		  noise,  statRQMC, varMC, secondsMC);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pSobol, randLMS, 
	    		  noise,  statRQMC, varMC, secondsMC);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pKor, randShift, 
	    		  noise,  statRQMC, varMC, secondsMC);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pKorBaker, randShift, 
	    		  noise,  statRQMC, varMC, secondsMC);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pLat1, randShift, 
	    		  noise,  statRQMC, varMC, secondsMC);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pLat1Baker, randShift, 
	    		  noise,  statRQMC, varMC, secondsMC);

          // Brownian gamma bridge sampling
	      System.out.println ("\n *****  Brownian gamma bridge sampling (BGBS)   *****\n");
		  asian.setProcess(spBGBS);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pSobol, randShift, 
	    		  noise,  statRQMC, varMC, secondsMC);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pSobol, randLMS, 
	    		  noise,  statRQMC, varMC, secondsMC);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pKor, randShift, 
	    		  noise,  statRQMC, varMC, secondsMC);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pKorBaker, randShift, 
	    		  noise,  statRQMC, varMC, secondsMC);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pLat1, randShift, 
	    		  noise,  statRQMC, varMC, secondsMC);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pLat1Baker, randShift, 
	    		  noise,  statRQMC, varMC, secondsMC);

          // Double gamma bridge sampling
	      System.out.println ("\n *****  Double gamma bridge sampling (DGBS)   *****\n");
		  asian.setProcess(spDGBS);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pSobol, randShift, 
	    		  noise,  statRQMC, varMC, secondsMC);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pSobol, randLMS, 
	    		  noise,  statRQMC, varMC, secondsMC);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pKor, randShift, 
	    		  noise,  statRQMC, varMC, secondsMC);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pKorBaker, randShift, 
	    		  noise,  statRQMC, varMC, secondsMC);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pLat1, randShift, 
	    		  noise,  statRQMC, varMC, secondsMC);
	      RQMCExperiment.simulReplicatesRQMCDefaultReportCompare (asian, m, pLat1Baker, randShift, 
	    		  noise,  statRQMC, varMC, secondsMC);
	      
*/	      
          // Estimating the derivative w.r.t. nu.
	      System.out.println ("\n *****   Derivative estimation   *****\n");

		  Tally statDiff = new Tally("Stats on difference for MC");		
		  statDiff.setConfidenceIntervalStudent();
		  statRQMC.setConfidenceIntervalStudent();
		  System.out.println(pSobol.toString());
		  System.out.println(randLMS.toString());
		  AsianOption asian2 = new AsianOption (r, d, T1, T, K);

		  double delta = 0.01;
		  StochasticProcess spDGBSdelta = new GeometricVarianceGammaProcess(s0, r,
				  new VarianceGammaProcessDiff(0.0,  theta,  sigma, nu + delta, 
					   new GammaProcessSymmetricalBridge(0.0, dummy, dummy, noise), 
					   new GammaProcessSymmetricalBridge(0.0, dummy, dummy, noise)));
		  asian.setProcess(spDGBS);
		  asian2.setProcess(spDGBSdelta);
	      RQMCExperiment.simulFDReplicatesRQMC (asian, asian2, delta, m, pSobol, randLMS, 
	    		  noise,  statRQMC);
	      System.out.println ("delta = " + delta);
		  System.out.println(statRQMC.report(0.95, 6));
	      System.out.println ("Variance per run: " + statRQMC.variance() * pSobol.getNumPoints() + "\n");
		  
		  delta = 0.001;
		  spDGBSdelta = new GeometricVarianceGammaProcess(s0, r,
				  new VarianceGammaProcessDiff(0.0,  theta,  sigma, nu + delta, 
					   new GammaProcessSymmetricalBridge(0.0, dummy, dummy, noise), 
					   new GammaProcessSymmetricalBridge(0.0, dummy, dummy, noise)));
		  asian.setProcess(spDGBS);
		  asian2.setProcess(spDGBSdelta);
	      RQMCExperiment.simulFDReplicatesRQMC (asian, asian2, delta, m, pSobol, randLMS, 
	    		  noise,  statRQMC);
	      System.out.println ("delta = " + delta);
		  System.out.println(statRQMC.report(0.95, 6));
	      System.out.println ("Variance per run: " + statRQMC.variance() * pSobol.getNumPoints() + "\n");
		  
		  delta = 0.0001;
		  spDGBSdelta = new GeometricVarianceGammaProcess(s0, r,
				  new VarianceGammaProcessDiff(0.0,  theta,  sigma, nu + delta, 
					   new GammaProcessSymmetricalBridge(0.0, dummy, dummy, noise), 
					   new GammaProcessSymmetricalBridge(0.0, dummy, dummy, noise)));
		  asian.setProcess(spDGBS);
		  asian2.setProcess(spDGBSdelta);
	      RQMCExperiment.simulFDReplicatesRQMC (asian, asian2, delta, m, pSobol, randLMS, 
	    		  noise,  statRQMC);
	      System.out.println ("delta = " + delta);
		  System.out.println(statRQMC.report(0.95, 6));
	      System.out.println ("Variance per run: " + statRQMC.variance() * pSobol.getNumPoints() + "\n");
		  
		    // To estimate a derivative via a finite difference.
		  MonteCarloExperiment.simulFDReplicatesCRN (asian, asian2, delta, n, noise, statDiff);
	      System.out.println ("Ordinary MC with CRNs");
		  System.out.println(statDiff.report(0.95, 6));
	      System.out.println ("Variance per run: " + statDiff.variance() + "\n");
				
		  MonteCarloExperiment.simulFDReplicatesIRN (asian, asian2, delta, n, noise, statDiff);
	      System.out.println ("Ordinary MC with IRNs");
		  System.out.println(statDiff.report(0.95, 6));
	      System.out.println ("Variance per run: " + statDiff.variance() + "\n");
			
	      System.out.println ("***   THE END   ***\n");	   
	   }
	}
