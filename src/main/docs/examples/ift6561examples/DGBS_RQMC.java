package ift6561examples;

// Exercice 2(a)
// package travail_final;

import umontreal.ssj.hups.BakerTransformedPointSet;
import umontreal.ssj.hups.DigitalNet;
import umontreal.ssj.hups.KorobovLattice;
import umontreal.ssj.hups.LMScrambleShift;
import umontreal.ssj.hups.PointSetIterator;
import umontreal.ssj.hups.RandomShift;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stochprocess.GammaProcess;


public class DGBS_RQMC {
	double strike;    // Strike price.
	   int d;            // Number of observation times.
	   double r ; // taux d inetret
	   double sigma;// volatilit� du MB
	   double theta;// moyenne du MB
	   double s0;// Log of the GBM process: logS[t] = log (S[t]).
	   double nu;// Variance du processus Gamma
	   double discount;
	   double T;// horizon de temps
	  
	   double[]  z;
	   public DGBS_RQMC(double nu,double theta,double r, double sigma, double strike, double s0, int d, double T,double[]  z) {
	       this.strike = strike;
	       this.d = d;
	       this.r=r;
	       this.s0=s0;
	       this.nu=nu;
	       this.sigma=sigma;
	       this.theta=theta;
	       this.T=T;
	       this.s0=s0;
	       discount = Math.exp (-r * T);
	       this.z=z;
	     }   
		   
	     public double getPayoff (double[] S ) {	       double average =0.0;  // Average of the GBM process.

		       for (int i  = 1; i <= d; i++) average=average +S[i];
	            average =average /d;  
		       if (average > strike) {return discount * (average - strike);}
		       else return 0.0; 
		   }




	     public void    simulateNrunsDGBS  (int n , RandomStream Stream, Tally  statDGBS,double [] z)  {
	    	    statDGBS.init();
	    		  
	    		  
	    		  int[] p=new int [d+1] ;
	    		  for (int l=0; l<=d; l++) p[l]=l;
	    		  double w=Math.log(1-theta*nu-sigma*sigma*0.5*nu)/nu;
	    		  
	    		  double mup=(Math.sqrt(theta*theta+(2*sigma*sigma/nu))+theta)/2;
	    		   double mun=(Math.sqrt(theta*theta+(2*sigma*sigma/nu))-theta)/2;
	    		   double volp=mup*mup*nu;
	    		   double voln=mun*mun*nu; 
	    		   GammaProcess vpos =new GammaProcess(0,mup,volp,Stream) ;
	    		   GammaProcess gneg =new GammaProcess(0,mun,voln,Stream) ;

	    	      for (int i = 0; i <n; i++) {
	    	   	   double[] S1=new double [d+1] ;
	    	   	   double[] S2=new double [d+1] ;
	    	   	   double[] S=new double [d+1] ;
	    	   	   vpos.setObservationTimes(z, d);
	    	   	   vpos.generatePath(Stream);
	    	   	   vpos.getSubpath(S1, p);
	    	   	   gneg.setObservationTimes(z, d);
	    	   	   gneg.generatePath(Stream);
	    	   	   gneg.getSubpath(S2, p);	
	    			   for (int l=1; l<=d; l++)    S[l]=s0*Math.exp(r*z[l]+S1[l]-S2[l]+w*z[l]);
	    			   double payoff=getPayoff(S);	 	 
	    			   statDGBS.add (payoff); 
	    	          Stream.resetNextSubstream();  }	
	    	         }

	     public void simulateQMC1 (int m, DigitalNet p,
	    	 		RandomStream noise, Tally statQMC, Tally statDGBS,double [] z)  {
	    		
	    		
	    		   statQMC.init();
	    		   PointSetIterator stream = p.iterator ();
	    		   RandomShift l= new RandomShift(noise);
	    		   
	    	 		for (int i=0; i< m; i++) {
	    	 		// d�calade al�atoire 
	    	 		l.randomize(p);
	    	 	    stream.resetStartStream();
	    	 		simulateNrunsDGBS (p.getNumPoints(),stream,  statDGBS, z);
	    	 		statQMC.add (statDGBS.average()); 
	    	 		
	    	 		}
	    	 		}
	    	public void simulateQMC2 (int m, DigitalNet p,
	    	 		RandomStream noise, Tally statQMC2,Tally  statDGBS,double [] z)  {
	    		    System.out.println(0.0);
	    		    statQMC2.init();
	    		    // Cette m�thode fait un left matrix scramble 
	    		    //suivi d'un d�calage digital al�atoire
	    		    LMScrambleShift  l=new LMScrambleShift(noise) ;
	    		    PointSetIterator stream = p.iterator ();
	    	 		for (int i=0; i< m; i++) {
	    	 		l.randomize(p);
	    	 		stream.resetStartStream();
	    	 		simulateNrunsDGBS (p.getNumPoints(),stream,  statDGBS, z);
	    	 		statQMC2.add (statDGBS.average());
	    	 		}	
	    	}
	    	//RQMC avec korobov Lattice et d�calage al�atoire 
	    	public void simulateQMC3 (int m, KorobovLattice p,
	    		RandomStream noise, Tally statQMC3,Tally  statDGBS,double [] z)  {
	    	    // d�calage al�atoire seulement
	    	    statQMC3.init();
	    	    RandomShift l= new RandomShift(noise);
	    	    PointSetIterator stream = p.iterator ();
	    		for (int i=0; i< m; i++) {
	    		// d�calade al�atoire 
	    		l.randomize(p);
	    		stream.resetStartStream();
	    		simulateNrunsDGBS (p.getNumPoints(),stream,  statDGBS, z);
	    		statQMC3.add (statDGBS.average());
	    			}
	    	}
	    	//RQMC avec korobov Lattice et d�calage al�atoire et la transformation du patissier
	    	public void simulateQMC4 (int m, KorobovLattice p,
	    		RandomStream noise, Tally statQMC4,Tally  statDGBS,double [] z)  {
                  
	    	    statQMC4.init();
	    	    PointSetIterator stream = p.iterator ();
	    	    RandomShift l= new RandomShift(noise);

	    	    
	    		for (int i=0; i< m; i++) {
	    		// d�calade al�atoire 
	    	    l.randomize(p);
	    		//p_tr.addRandomShift (0, p_tr.getDimension(), noise);
	    	    BakerTransformedPointSet pb = new BakerTransformedPointSet (p);

	    		stream.resetStartStream();
	    		simulateNrunsDGBS (pb.getNumPoints(),stream,  statDGBS, z);
	    		statQMC4.add (statDGBS.average());
	    			
	    			}}



	    	public static void main(String[] args) {
	    		int n=1000;
	    	    int d =16;
	    	    double[] z = new double[d+1];
	    	    for (int t=0; t<=d; t++)
	    	       z[t] = (double)t / (double)d; 
	    	    double strike=101;   
	    		double r=0.1 ; 
	    		double sigma=0.12136;
	    		double theta=-0.1436;
	    		double s0=100;
	    		double nu=0.3;
	    		double T=1;
	    		int m=32;// nombre de r�p�titions de RQMC
	    		int s=32;//la dimension
	    		DigitalNet p = new SobolSequence (14,31,32);//2^14 points
	    		KorobovLattice k=new  KorobovLattice(16381,5693,32);

	    	    Tally statDGBS = new Tally ("Stats of the method statDGBS");
	    	   
	    	    Tally statQMC1=new Tally ("Stats of the method with sobol net and  digital shift ");
	    	   
	    	    DGBS_RQMC Acrn = new  DGBS_RQMC(  nu, theta, r,  sigma,  strike, s0,  d,  T,z);
	    	    Acrn.simulateNrunsDGBS  (n ,  new MRG32k3a(),   statDGBS, z) ;
	    		double v1=statDGBS.variance();
	    	    //1//....D�commenter pour voir les r�sultats de RQMC avec Sobol net et d�calage al�atoire....//

	    		//statDGBS.init();
	    		//Acrn.simulateQMC2 ( m, p,
	    		    	   //new MRG32k3a(),statQMC1,statDGBS,z) ; 
	    		    //System.out.println(("R�sultats avec la M�thode Sobol Net avec d�calage al�atoire"));
	    		   // System.out.println("Pour la m�thode DGBS");
	    		   // statQMC1.setConfidenceIntervalStudent();
	    		    //System.out.println(statQMC1.report());
	    		    //System.out.println ("Rapport de r�duction de la variance pour la m�thode DGBS: "  +v1/( p.getNumPoints()*statQMC1.variance()));
	    		
	    	    //2//....D�commenter pour voir les r�sultats de RQMC avec Sobol net et d�calage al�atoire et left matrix scramble....//

	    		//statDGBS.init();
	    		//Acrn.simulateQMC2 ( m, p,
	    		    	//   new MRG32k3a(),statQMC1,statDGBS,z) ; 
	    		    //System.out.println(("R�sultats avec la M�thode Sobol Net avec d�calage al�atoire et left matrix sramble"));
	    		   // System.out.println("Pour la m�thode BGBS");
	    		   // statQMC1.setConfidenceIntervalStudent();
	    		    //System.out.println(statQMC1.report());
	    		    //System.out.println ("Rapport de r�duction de la variance pour la m�thode DGBS: "  +v1/( p.getNumPoints()*statQMC1.variance()));

	    	    //3//....D�commenter pour voir les r�sultats de RQMC avec RQMC avec korobov+random shift....//

	    	    
	    	    //statDGBS.init();
	    	    //Acrn.simulateQMC3 (m, k,
	    	    	//   new MRG32k3a(),statQMC1,statDGBS,z) ; 
	    	    //System.out.println(("R�sultats avec korobov lattice et d�calage al�atoire"));
	    	    //System.out.println("Pour la m�thode BGBS");
	    	    //statQMC1.setConfidenceIntervalStudent();
	    	    //System.out.println(statQMC1.report());
	    	    //System.out.println ("Rapport de r�duction de la variance pour la m�thode BGBS: "  +v1/( k.getNumPoints()*statQMC1.variance()));
	    	    

	    	    
	    	    
	    	  // 4 //....D�commenter pour voir les r�sultats de RQMC avec korobov+random shift+trans du pattissier....//
	    	    
	    	    
	    		//Acrn.simulateQMC4 (m, k,
	    		//   	   new MRG32k3a(),statQMC4,statDGBS4,z) ; 
	    		//    System.out.println(("R�sultats avec korobov lattice, la transformation du p�tissier et d�calage al�atoire"));
	    		//    System.out.println("Pour la m�thode DGBS");
	    		//   statQMC4.setConfidenceIntervalStudent();
	    		//   System.out.println(statQMC4.report());
	    		//    System.out.println ("Rapport de r�duction de la variance pour la m�thode BGBS: "  +v1/( k.getNumPoints()*statQMC4.variance()));
	    		    
	    		    

	    	   
	    	   
	    	              }
	    	}













