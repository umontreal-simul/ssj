/*
 * Class:        VarianceGammaProcessDiffPCA
 * Description:  
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since        2008
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package umontreal.ssj.stochprocess;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.randvar.*;

/**
 * Same as  @ref VarianceGammaProcessDiff, but the two inner
 * @ref GammaProcess ’es are of PCA type. Also, `generatePath(double[]
 * uniforms01)` distributes the uniform random variates to the
 * @ref GammaProcessPCA ’s according to their eigenvalues, i.e. the
 * @ref GammaProcessPCA with the higher eigenvalue gets the next uniform
 * random number. If one should decide to create a
 * @ref VarianceGammaProcessDiffPCA by giving two  @ref GammaProcessPCA ’s to
 * an objet of the class  @ref VarianceGammaProcessDiff, the uniform random
 * numbers would not be given this way to the  @ref GammaProcessPCA ’s; this
 * might give less variance reduction when used with QMC.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class VarianceGammaProcessDiffPCA extends VarianceGammaProcessDiff {
    int[] indexEigenUp;
    int[] indexEigenDw;

   /**
    * Constructs a new  @ref VarianceGammaProcessDiffPCA with parameters
    * @f$\theta= \mathtt{theta}@f$, @f$\sigma= \mathtt{sigma}@f$,
    * @f$\nu= \mathtt{nu}@f$ and initial value @f$S(t_0) =
    * \mathtt{s0}@f$. There is only one
    * @ref umontreal.ssj.rng.RandomStream here which is used for the two
    * inner  @ref GammaProcessPCA ’s. The other parameters are set as in
    * @ref VarianceGammaProcessDiff.
    */
   public VarianceGammaProcessDiffPCA (double s0, double theta,
                                       double sigma, double nu,
                                       RandomStream stream) {
     super(s0, theta, sigma, nu, 
	  new GammaProcessPCA (0.0, 1.0, 1.0, stream),
	  new GammaProcessPCA (0.0, 1.0, 1.0, stream));
    // Params mu, nu of the 2 gamma processes are redefined in init()
    // which will be called after a call to 'setObservTimes'
}

   /**
    * Constructs a new  @ref VarianceGammaProcessDiffPCA with parameters
    * @f$\theta= \mathtt{theta}@f$, @f$\sigma= \mathtt{sigma}@f$,
    * @f$\nu= \mathtt{nu}@f$ and initial value @f$S(t_0) =
    * \mathtt{s0}@f$. As in  @ref VarianceGammaProcessDiff, the
    * @ref umontreal.ssj.rng.RandomStream of `gneg` is replaced by the one
    * of `gpos` to avoid any confusion.
    */
   public VarianceGammaProcessDiffPCA (double s0, double theta,
                                       double sigma, double nu,
                                       GammaProcessPCA gpos,
                                       GammaProcessPCA gneg) {
    super(s0, theta, sigma, nu, gpos, gneg); // from VarianceGammaProcessDiff
    // Params mu, nu of the 2 gamma processes are redefined in init()
    // which will be called after a call to 'setObservTimes'
}

   /**
    * This method is not implemented is this class since the path cannot
    * be generated sequentially.
    */
   public double nextObservation() {
        throw new UnsupportedOperationException 
        ("Impossible with PCA, use generatePath() instead.");
    }
 

   public double[] generatePath() {
        double[] u = new double[2*d];
        for(int i =0; i < 2*d; i++)
            u[i] = getStream().nextDouble();
        return generatePath(u);
    }

   public double[] generatePath(double[] uniform01)  {
        int dd = uniform01.length;
        int d = dd / 2;

        if(dd % 2 != 0){
            throw new IllegalArgumentException (
                     "The Array uniform01 must have a even length");
        }

        double[] QMCpointsUP = new double[d];
        double[] QMCpointsDW = new double[d];

        for(int i = 0; i < d; i++){
             QMCpointsUP[i] = uniform01[ indexEigenUp[i] ];
             QMCpointsDW[i] = uniform01[ indexEigenDw[i] ];
        }
        gpos.resetStartProcess();
        gneg.resetStartProcess();

        double[] pathUP = gpos.generatePath(QMCpointsUP);
        double[] pathDOWN = gneg.generatePath(QMCpointsDW);

        for (int i=0; i<d; i++) {
           path[i+1] = x0 + pathUP[i+1] - pathDOWN[i+1];
        }
        observationIndex   = d;
        observationCounter = d;
        return path;
    }


   protected void init() {
        super.init ();  // from VarianceGammaProcessDiff
	if( observationTimesSet){
        // Two lines below (casts) should be reinstated after fix inheritance PCA/PCABridge.
	    double[] eigenValUp = ((GammaProcessPCA)gpos).getBMPCA().getSortedEigenvalues();
	    double[] eigenValDw = ((GammaProcessPCA)gneg).getBMPCA().getSortedEigenvalues();
	    indexEigenUp = new int[d];
	    indexEigenDw = new int[d];

	    int iUp = 0;
	    int iDw = 0;
	    for(int iQMC = 0; iQMC < 2*d; iQMC++){
		if(iUp == d) {indexEigenDw[iDw] = iQMC; iDw++;continue;}
        if(iDw == d) {indexEigenUp[iUp] = iQMC; iUp++;continue;}
        if( eigenValUp[iUp] >= eigenValDw[iDw] ){
		    indexEigenUp[iUp] = iQMC; 
		    iUp++;
		}
		else{
		    indexEigenDw[iDw] = iQMC; 
		    iDw++;
		}
	    }
	}

    }

}