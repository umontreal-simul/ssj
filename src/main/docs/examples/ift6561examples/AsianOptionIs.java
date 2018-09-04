package ift6561examples;
// Exercice 2 (c)


import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.randvar.GammaAcceptanceRejectionGen;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;

public class AsianOptionIs {
	
		   double strike;    // Strike price.
		   int d;            // Number of observation times.
		   double r ; // taux d inetret
		   double sigma;// volatilit� du MB
		   double theta;// moyenne du MB
		   double s0;
		   double nu;// Variance du processus Gamma
		   double T;// horizon de temps
		   double theta_para;
		   double discount;
		   RandomStream stream=new MRG32k3a();
		  

		  
         
          
		   public AsianOptionIs(double nu,double theta,double r, double sigma, double strike, double s0, int d, double T,double theta_para) {
		   this.strike = strike;
		   this.d=d;
		   this.sigma=sigma;
		   this.theta=theta;
		   this.T=T;
		   this.s0=s0;
           discount = Math.exp (-r * T);
           this.nu=nu;
           this.r=r;
           this.theta_para=theta_para;
}         
		   
		   
		   public void simulateIs ( RandomStream stream, Tally  statIS,Tally  statMC,double [] z,int d,
				   GammaAcceptanceRejectionGen gpos,GammaAcceptanceRejectionGen gneg, GammaAcceptanceRejectionGen gposIS,GammaAcceptanceRejectionGen gnegIS) {
			   double mup=(Math.sqrt(theta*theta+(2*sigma*sigma/nu))+theta)/2;
 			   double mun=(Math.sqrt(theta*theta+(2*sigma*sigma/nu))-theta)/2;
 			  double mupg=(1/nu);
 			  double volpg=(1/(mup*nu));
 			  double mung=(1/nu);
			  double volng=(1/(mun*nu));
 			   
 			  
			   // On commence par estime en utilisant monte Carlo Standard
			   double w=Math.log(1-theta*nu-sigma*sigma*0.5*nu)/nu;
			   //System.out.println(w);
			   double s=Math.log(strike/s0);
               //gneg.setObservationTimes(z, d);// z contient l'instant 0 et 1
			   //gneg.generatePath();
			   double y1=gneg.nextDouble();// On prend G+(1)
			   //double y3=GammaAcceptanceRejectionGen.nextDouble(stream,mun*mun/voln,mun/voln);
			   //gpos.setObservationTimes(z, d);// z contient l'instant 0 et 1
			   //gpos.generatePath();
			   double y2=gpos.nextDouble();// On prend G-(1)
                if (s<=r+w+y2-y1) statMC.add(discount*(s0*Math.exp(r+w+y2-y1)-strike));
			   else statMC.add(0.0);
			   // Estimation par importance sampling
			   //gnegIS.setObservationTimes(z, d);// z contient l'instant 0 et 1
			   //gnegIS.generatePath();
		        double y1IS=gnegIS.nextDouble();// On prend GIS+(1)
			   //gposIS.setObservationTimes(z, d);// z contient l'instant 0 et 1
			   //gposIS.generatePath();
	            double y2IS=gposIS.nextDouble();// On prend GIS-(1)
			   // Rapport de vraisemblance 
	            double L=Math.exp(theta_para*(y1IS-y2IS))*Math.pow((volpg/(volpg-theta_para)),mupg)*Math.pow((volng/(volng+theta_para)),mung);
                //System.out.println(volpg*volng*0.1029);

	           if (s<=r+w+y2IS-y1IS) statIS.add(L*discount*(s0*Math.exp(r+w+y2IS-y1IS)-strike));
               else statIS.add(0.0);
	           
	          
	 
			   }
		       public void simulateISNRUNS(int N, RandomStream stream, Tally  statIS,Tally  statMC,double [] z,int d){
                statIS.init();
                statMC.init();
               double mup=(Math.sqrt(theta*theta+(2*sigma*sigma/nu))+theta)/2;
 			   double mun=(Math.sqrt(theta*theta+(2*sigma*sigma/nu))-theta)/2;
 			  double mupg=(1/nu);
 			  double volpg=(1/(mup*nu));
 			   double mung=(1/nu);
			  double volng=(1/(mun*nu));
 			   GammaAcceptanceRejectionGen     gpos=new GammaAcceptanceRejectionGen (stream, new GammaDist(1/nu,1/(mup*nu)));
 			   GammaAcceptanceRejectionGen    gneg=new  GammaAcceptanceRejectionGen(stream, new GammaDist(1/nu,1/(mun*nu)));
 			   GammaAcceptanceRejectionGen    gposIS=new  GammaAcceptanceRejectionGen(stream,new GammaDist(mupg,volpg-theta_para)); 
 			   GammaAcceptanceRejectionGen   gnegIS=new  GammaAcceptanceRejectionGen(stream,new GammaDist(mung,volng+theta_para));
 			 // MathFunction f=mun*x+volp;
 			 // UnivariateRealFunction f = new SinFunction();
         		for (int l=1; l<=N; l++)simulateIs (  stream,  statIS, statMC, z,d,gpos,gneg,gposIS,gnegIS);
                System.out.println("Pour une valeur du strike �gale � " + strike );

         		statIS.setConfidenceIntervalStudent();
               System.out.println(statIS.report(0.95,3));
                statMC.setConfidenceIntervalStudent();
                System.out.println(statMC.report(0.95,3));
                System.out.println(" Le rapport de r�duction de la variance avec la m�thode IS est "+statMC.variance()/statIS.variance());
		   
		   }
         		public static void main(String[] args) {
         			int N=100000000;
         			
         		    int d =1;
         		    double[] z = new double[d+1];

         		   for (int t=1; t<=d; t++)
         		       z[t] = (double)t / (double)d; 
         		    
         			double r=0.1 ; 
         			double sigma=0.12136;
         			double theta=-0.1436;
         			double s0=100;
         			double nu=0.3;
         			double T=1;
         		    Tally statMC = new Tally ("Stats of the the estimation with standard MC ");
         		    Tally statIS = new Tally ("Stats of the the estimation with IS ");
         		    //k=180,theta= 25.5652
         		    //k=140,THETA=17.3194
         		    AsianOptionIs Ais = new AsianOptionIs(nu, theta, r,  sigma, 140,  s0,  d,  T, 17.3194);
         		    Ais.simulateISNRUNS (N,new MRG32k3a(), statIS, statMC, z,d);
     
		   
         		 }
}
         		    
