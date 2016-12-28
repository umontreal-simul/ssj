package probdistmulti.norta;

import umontreal.ssj.probdist.*;
import umontreal.ssj.probdistmulti.norta.*;

public class ExampleNortaInitDisc
{
    public static void main (String[] args) {
	final double rX = 0.43;              // Target rank correlation rX
	final double tr = 1.0 - 1.0e-6;	     // Quantile upper limit

	// Define the two marginal distributions
	DiscreteDistributionInt dist1 = new NegativeBinomialDist(15.68, 0.3861);
	DiscreteDistributionInt dist2 = new NegativeBinomialDist(60.21, 0.6211);

	NI1 ni1Obj = new NI1(rX, dist1, dist2, tr, 1.0e-4);
	System.out.println("Result with method NI1:  rho_Z = "
			   + String.format("%.14g", ni1Obj.computeCorr()));
	NI2a ni2aObj = new NI2a(rX, dist1, dist2, tr, 0.005, 1.0e-4);
	System.out.println("Result with method NI2b: rho_Z = "
			   + String.format("%.14g", ni2aObj.computeCorr()));
	NI2b ni2bObj = new NI2b(rX, dist1, dist2, tr, 5, 1.0e-4);
	System.out.println("Result with method NI2a: rho_Z = "
			   + String.format("%.14g", ni2bObj.computeCorr()));
	NI3 ni3Obj = new NI3(rX, dist1, dist2, tr, 1.0e-4);
	System.out.println("Result with method NI3:  rho_Z = "
			   + String.format("%.14g", ni3Obj.computeCorr()));
    }
}
