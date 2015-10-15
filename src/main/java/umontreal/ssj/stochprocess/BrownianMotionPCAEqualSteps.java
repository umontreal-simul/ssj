/*
 * Class:        BrownianMotionPCAEqualSteps
 * Description:  
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since

 * SSJ is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, either version 3 of the License, or
 * any later version.

 * SSJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * A copy of the GNU General Public License is available at
   <a href="http://www.gnu.org/licenses">GPL licence site</a>.
 */
package umontreal.ssj.stochprocess;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.randvar.*;

/**
 * Same as BrownianMotionPCA, but uses a trick to speed up the calculation
 * when the time steps are equidistant.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class BrownianMotionPCAEqualSteps extends BrownianMotion {

    double dt;
    protected double[][]  A;     // sigmaCov = AA' (PCA decomposition).
    protected double[]    z;     // vector of standard normals.
    protected boolean     isDecompPCA;
    protected double[]    sortedEigenvalues;

    /**
     * Constructs a new `BrownianMotionPCAEqualSteps`.
     */
    public BrownianMotionPCAEqualSteps (double x0, double mu, double sigma,
                              RandomStream stream) {
        super (x0, mu, sigma, stream);
        isDecompPCA = false;
    }

    /**
     * Constructs a new `BrownianMotionPCAEqualSteps`.
     */
    public BrownianMotionPCAEqualSteps (double x0, double mu, double sigma,
                              NormalGen gen) {
        super (x0, mu, sigma, gen);
        isDecompPCA = false;
    }



    public double nextObservation() {
        throw new UnsupportedOperationException 
	    ("nextObservation() not defined for PCA.");
    }



    public double[] generatePath() {
       if(!isDecompPCA) {init();}  // if the decomposition is not done, do it...
       for (int j = 0; j < d; j++)
           z[j] = gen.nextDouble ();
       for (int j = 0; j < d; j++) {
           double sum = 0.0;
           for (int k = 0; k < d; k++)
               sum += A[j][k] * z[k];
           path[j+1] = x0 + mu * t[j+1] + sum;
       }
       observationIndex   = d;
       observationCounter = d;
       return path;
    }

    public double[] generatePath(double[] QMCpointsBM) {
       if(!isDecompPCA) {init();}  // if the decomposition is not done, do it...
       for (int j = 0; j < d; j++)
           z[j] = NormalDist.inverseF01(QMCpointsBM[j]);
       for (int j = 0; j < d; j++) {
           double sum = 0.0;
           for (int k = 0; k < d; k++)
               sum += A[j][k] * z[k];
           path[j+1] = x0 + mu * t[j+1] + sum;
       }
       observationIndex   = d;
       observationCounter = d;
       return path;
    }


    public void setObservationTimes(double[] t, int d){
	    super.setObservationTimes(t,d);
	    this.dt = t[1] - t[0];
	    for(int i = 1; i < d; i++)
		if( Math.abs((t[i+1] - t[i])/dt - 1.0) > 1e-7 )
		    throw new IllegalArgumentException("Not equidistant times");
    }


    public void setObservationTimes(double dt, int d){
	this.dt = dt;
	super.setObservationTimes(dt,d);
    }


    protected void init() {
	super.init();
	if(observationTimesSet){
	    final double twoOverSqrt2dP1 = 2.0/Math.sqrt(2.0*d+1.0);
	    final double piOver2dP1 = Math.PI/(2.0*d+1.0);

	    z = new double[d];
	    // A contains the eigenvectors (as columns), times sqrt(eigenvalues).
	    A = new double[d][d];
	    sortedEigenvalues = new double[d];
	    for(int ic = 1; ic <= d; ic++){
		final double tempSin = Math.sin( (2*ic-1)*piOver2dP1/2.0 );
		sortedEigenvalues[ic-1] = dt/4.0/tempSin/tempSin * sigma*sigma;
		for(int ir = 1; ir <= d; ir++){
		    A[ir-1][ic-1] = twoOverSqrt2dP1 * Math.sin( (2*ic-1)*piOver2dP1*ir );
		    A[ir-1][ic-1] *= Math.sqrt( sortedEigenvalues[ic-1] );
		}
	    }
	    double[][] AA = new double[d][d];
	    for(int ic = 0; ic < d; ic++){
		for(int ir = 0; ir < d; ir++){
		    double sum = 0.0;
		    for(int k = 0; k < d; k++) sum += A[ir][k]*A[ic][k];
		    AA[ir][ic] = sum;
		}
	    }
	    isDecompPCA = true;
	}
    }
       public double[] getSortedEigenvalues(){
       return sortedEigenvalues;
   }
}