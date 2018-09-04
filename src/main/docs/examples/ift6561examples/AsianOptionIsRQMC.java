package ift6561examples;

import umontreal.ssj.hups.DigitalNet;
import umontreal.ssj.hups.PointSetIterator;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;

	
public class AsianOptionIsRQMC extends AsianOptionIs2 {
	
	
	public AsianOptionIsRQMC(double nu, double theta, double r, double sigma,
		double strike, double s0, int d, double T, double theta_para) {
		super(nu, theta, r, sigma, strike, s0, d, T, theta_para);
		
	}

    public void simulateQMC(int m, DigitalNet p,
 		RandomStream noise, Tally statQMC,Tally statIS)  {
 	    // d�calage al�atoire seulement
	    statQMC.init();
	    PointSetIterator stream = p.iterator ();
	    int d =1;
		    double[] z = new double[d+1];

		   for (int t=1; t<=d; t++)
		       z[t] = (double)t / (double)d; 
 		for (int i=0; i< m; i++) {
 	   
 	    p.leftMatrixScramble(noise);
 		p.addRandomShift (0, p.getDimension(), noise);
 	    stream.resetStartStream();
 		simulateISNRUNS(p.getNumPoints(),  stream,  statIS, z, d);
 		statQMC.add (statIS.average());
 		 }	}	


public static void main(String[] args) {
	int d=1;
	 DigitalNet p = new SobolSequence (14,31,32);
	 int m=20;
	 double r=0.1 ; 
		double sigma=0.12136;
		double theta=-0.1436;
		double s0=100;
		double nu=0.3;
		double T=1;
	    Tally statIS = new Tally ("Stats of the the estimation with standard IS only ");
	    Tally statRQMC = new Tally ("Stats of the the estimation with IS and RQMC ");
	    AsianOptionIsRQMC Ais = new AsianOptionIsRQMC(nu, theta, r,  sigma, 140,  s0,  d,  T, 17.3194);
        Ais.simulateQMC(m, p, new MRG32k3a(), statRQMC, statIS);
        System.out.println("pour K="+140);
        statRQMC.setConfidenceIntervalStudent();
        System.out.println(statRQMC.report());
        System.out.println("Rapport addictionnel de r�duction de la variance ="+ 0.113*0.113/(statRQMC.variance()*p.getNumPoints()));
        System.out.println("pour K="+180);
        AsianOptionIsRQMC Ais2 = new AsianOptionIsRQMC(nu, theta, r,  sigma, 180,  s0,  d,  T, 25.5652);
        Ais2.simulateQMC(m, p, new MRG32k3a(), statRQMC, statIS);
        statRQMC.setConfidenceIntervalStudent();
        System.out.println(statRQMC.report());
        System.out.println("Rapport addictionnel de r�duction de la variance ="+ (1.7E-4*1.7E-4)/(p.getNumPoints()*statRQMC.variance()));

}	








}	