package ift6561examples;
import umontreal.ssj.stochprocess.*;
import umontreal.ssj.rng.*;

public class VarianceGammaProcessAlternate extends VarianceGammaProcess {


	public VarianceGammaProcessAlternate (double s0, double theta, double sigma,
	                              double nu, RandomStream stream) {
	        super (s0, theta, sigma, nu, stream);
	     }


    public VarianceGammaProcessAlternate (double s0, BrownianMotion BM,
	                              GammaProcess Gamma) {
		   super (s0, BM, Gamma);
         }

	   /**
	    * Generates the path in an alternate way, and returns the path. 
	    * Equivalent to calling the old  `generatePath(double[] uniform01)`.
	    */
       // @overload
	   public double[] generatePath() {
		   int d = getNbObservationTimes();
		   double[] points = new double[2*d]; 
		   RandomStream stream = getStream();
		   for (int j=0; j<2*d; j++)
			   points[j] = stream.nextDouble(); 
		   return generatePath (points);
		}
   }
