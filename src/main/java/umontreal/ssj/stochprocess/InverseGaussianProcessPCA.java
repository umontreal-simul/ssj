/*
 * Class:        InverseGaussianProcessPCA
 * Description:  
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since
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
 * Approximates a principal component analysis (PCA) decomposition of the
 * `InverseGaussianProcess`. The PCA decomposition of a
 * @ref BrownianMotionPCA with a covariance matrix identical to the one of
 * our `InverseGaussianProcess` is used to generate the path of our
 * `InverseGaussianProcess` @cite fLEC08a&thinsp;. Such a path is a perfectly
 * random path and it is hoped that it will provide reduction in the
 * simulation variance when using quasi-Monte Carlo.
 *
 * The method `nextObservation()` cannot be used with PCA decompositions
 * since the whole path must be generated at once.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class InverseGaussianProcessPCA extends InverseGaussianProcess {

    protected BrownianMotionPCA bmPCA;

   /**
    * Constructs a new `InverseGaussianProcessPCA`. The initial value `s0`
    * will be overridden by @f$t[0]@f$ when the observation times are set.
    */
   public InverseGaussianProcessPCA (double s0, double delta, double gamma,
                                     RandomStream stream) {
        super(s0, delta, gamma, stream);
        bmPCA = new BrownianMotionPCA(0.,0.,delta,stream);
        numberOfRandomStreams = 1;
    }


   public double[] generatePath () {
        double[] uniformIncrement = new double[d];
        double[] BMpath = bmPCA.generatePath();

        for(int i = 0; i < d; i++)
        {
            double dt    = bmPCA.mudt[i]; //bmTime[i + 1] - bmTime[i];
            double sigma = bmPCA.sigmasqrdt[i] ;//Math.sqrt(dt) * bmSigma;
            uniformIncrement[i] =
            NormalDistQuick.cdf01( (BMpath[i+1] - BMpath[i] - bmPCA.mu * dt)/sigma );
        }
        path[0] = x0;
        for(int i = 0; i < d; i++)
            path[i+1] = path[i] +
                InverseGaussianDist.inverseF(imu[i], ilam[i], uniformIncrement[i]);

        observationIndex   = d;
        observationCounter = d;
        return path;
    }

/**
 * Instead of using the internal stream to generate the path, uses an array
 * of uniforms @f$U[0,1)@f$. The length of the array should be equal to the
 * length of the number of periods in the observation times. This method is
 * useful for  @ref NormalInverseGaussianProcess.
 */
public double[] generatePath (double[] uniforms01) {
        double[] uniformIncrement = new double[d];
        double[] BMpath = bmPCA.generatePath(uniforms01);

        for(int i = 0; i < d; i++)
        {
            double dt    = bmPCA.mudt[i]; //bmTime[i + 1] - bmTime[i];
            double sigma = bmPCA.sigmasqrdt[i] ;//Math.sqrt(dt) * bmSigma;
            uniformIncrement[i] =
            NormalDistQuick.cdf01( (BMpath[i+1] - BMpath[i] - bmPCA.mu * dt)/sigma );
        }
        path[0] = x0;
        for(int i = 0; i < d; i++)
            path[i+1] = path[i] +
                InverseGaussianDist.inverseF(imu[i], ilam[i], uniformIncrement[i]);

        observationIndex   = d;
        observationCounter = d;
        return path;
    }

   /**
    * Not implementable for PCA.
    */
   public double nextObservation() {
        throw new UnsupportedOperationException("Not implementable for PCA.");
    }

   /**
    * Sets the observation times of both the
    * @ref InverseGaussianProcessPCA and the inner <br>
    * @ref BrownianMotionPCA.
    */
   public void setObservationTimes (double t[], int d) {
        super.setObservationTimes(t,d);
        bmPCA.setObservationTimes(t,d);
    }


   public RandomStream getStream() {
        if( stream != bmPCA.getStream() )
            throw new IllegalStateException("Two different streams or more are present");
        return stream;
    }


    public void setStream (RandomStream stream) {
        super.setStream(stream);
        bmPCA.setStream(stream);
    }

/**
 * Sets the brownian motion PCA. The observation times will be overriden when
 * the method `observationTimes()` is called on the
 * @ref InverseGaussianProcessPCA.
 */
public void setBrownianMotionPCA (BrownianMotionPCA bmPCA) {
        this.bmPCA = bmPCA;
    }

   /**
    * Returns the  @ref BrownianMotionPCA.
    */
   public BrownianMotion getBrownianMotionPCA() {
        return bmPCA; 
    }

}