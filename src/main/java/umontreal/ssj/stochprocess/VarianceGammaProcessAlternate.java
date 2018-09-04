package umontreal.ssj.stochprocess;

import umontreal.ssj.rng.*;

/**
 * This is a @ref VarianceGammaProcess for which the successive random numbers
 * are used in a different order to generate the sample path.
 * The first one is used for the first generated value of the gamma process, 
 * the second one for the first generated value of the Brownian process, 
 * the third one for the second generated value of the gamma process,
 * the fourth one for the second value of the Brownian process, and so on.
 * Only the order in which the uniform random numbers are used in the method `generatePath` differs.
 * These numbers are generated at the beginning and then reordered.
 * This can make a difference when we use RQMC methods.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
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
	    * Generates the sample path by using the uniform random numbers in an alternate way, 
	    * and returns the path of the VG process. 
	    */
	   public double[] generatePath() {
		   int d = getNbObservationTimes();
		   double[] points = new double[2*d]; 
		   RandomStream stream = getStream();
		   for (int j=0; j<2*d; j++)
			   points[j] = stream.nextDouble(); 
		   return generatePath (points);
		}
   }
