/* 
 * The conditional bridge distribution is known @cite fWEB03a,
 * but is not integrable, therefore there can be no bridge method
 * using simply inversion, with a single 
 * \externalclass{umontreal.ssj.rng}{RandomStream}.
 */

/*
 * Class:        InverseGaussianProcessBridge
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
 * Samples the path by bridge sampling: first finding the process value at
 * the final time and then the middle time, etc. The method
 * `nextObservation()` returns the path value in that non-sequential order.
 * This class uses two  @ref umontreal.ssj.rng.RandomStream ’s to generate a
 * path @cite fWEB03a&thinsp;.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class InverseGaussianProcessBridge extends InverseGaussianProcessMSH {

    // Careful: mu and lambda are completely different from parent class.
    protected double[] imu2;
    protected double[] imuLambdaZ;
    protected double[] imuOver2LambdaZ;

    protected int[] wIndexList;
    protected int   bridgeCounter = -1; // Before 1st observ

   /**
    * Constructs a new `InverseGaussianProcessBridge`. The initial value
    * `s0` will be overridden by @f$t[0]@f$ when the observation times are
    * set.
    */
   public InverseGaussianProcessBridge (double s0, double delta, 
                                        double gamma, RandomStream stream,
                                        RandomStream otherStream) {
        super(s0, delta, gamma, stream, otherStream);
        numberOfRandomStreams = 2;
    }

   /**
    * Generates the path. The two inner
    * @ref umontreal.ssj.rng.RandomStream ’s are sampled alternatively.
    */
   public double[] generatePath() {
        bridgeCounter    = -1;
        observationIndex =  0;
        path[0]          = x0;
        for (int j = 0; j < d; j++)   nextObservation();
        return path;
    }

   /**
    * Instead of using the internal streams to generate the path, it uses
    * two arrays of uniforms @f$U[0,1)@f$. The length of the arrays
    * `unifNorm` and `unifOther` should be equal to the number of time
    * steps, excluding @f$t_0@f$.
    */
   public double[] generatePath (double[] unifNorm, double[] unifOther) {
        double s = x0;
        RandomStream cacheStreamNormal = stream;
        RandomStream cacheStreamOther  = otherStream;
        stream = new NonRandomStream(unifNorm);
        normalGen.setStream(stream);
        otherStream = new NonRandomStream(unifOther);

        bridgeCounter    = -1;
        observationIndex =  0;
        path[0]          = x0;
        for (int j = 0; j < d; j++)   nextObservation();

        stream = cacheStreamNormal;
        normalGen.setStream(stream);
        otherStream = cacheStreamOther;
        return path;
    }

   /**
    * Returns the next observation in the bridge order, not the sequential
    * order.
    */
   public double nextObservation() {
        double s;
        if (bridgeCounter == -1) 
        {
            double temp = delta * (t[d] - t[0]);
            s = x0 + InverseGaussianMSHGen.nextDouble(otherStream,
                                                normalGen, temp/gamma, temp*temp);
            bridgeCounter    = 0;
            observationIndex = d;
        }
        else
        {
            int j = bridgeCounter * 3;
            int oldIndexL = wIndexList[j];
            int newIndex  = wIndexList[j + 1];
            int oldIndexR = wIndexList[j + 2];

            // Use the fact that \chi^2_1 is equivalent to a normal squared.
            double q = normalGen.nextDouble();
            q *= q;

            double z        = path[oldIndexR];
            double mu       = imu[newIndex];
            double muLambda = imuLambdaZ[newIndex]/z;
            double mu2      = imu2[newIndex];
            double muOver2Lambda = imuOver2LambdaZ[newIndex]*z;

            double root = mu + muOver2Lambda*mu*q -
            muOver2Lambda*Math.sqrt(4.*muLambda*q + mu2*q*q);
            double probabilityRoot1 = mu*(1.+root)/(1.+mu)/(root+mu);
            // Check if reject first root for 2nd one: root2=mu^2/root1.
            if (otherStream.nextDouble() > probabilityRoot1)
                root = mu2/root;
            s = path[oldIndexL] + (path[oldIndexR] - path[oldIndexL]) /(1.0 + root);
            bridgeCounter++;
            observationIndex = newIndex;
        }
        observationCounter = bridgeCounter + 1;
        path[observationIndex] = s;
        return s;
    }


   public void resetStartProcess () {
        observationIndex   = 0;
        observationCounter = 0;
        bridgeCounter = -1;
    }


    protected void init () {
        // imu[] etc. in super.init() will be overriden here.
        // Necessary nonetheless.
        super.init();

        double tauX;
        double tauY;
        double mu;
        double lambdaZ;
        if (observationTimesSet) {
            wIndexList  = new int[3*d];

            int[] ptIndex = new int[d+1];
            int   indexCounter = 0;
            int   newIndex, oldIndexL, oldIndexR;

            ptIndex[0] = 0;
            ptIndex[1] = d;

            // Careful: mu and lambda are completely different from parent class.
            imu            = new double[d+1];
            ilam           = new double[d+1];
            imu2           = new double[d+1];
            imuLambdaZ     = new double[d+1];
            imuOver2LambdaZ= new double[d+1];

            for (int powOfTwo = 1; powOfTwo <= d/2; powOfTwo *= 2) {
                /* Make room in the indexing array "ptIndex" */
                for (int j = powOfTwo; j >= 1; j--) ptIndex[2*j] = ptIndex[j];

                /* Insert new indices and Calculate constants */
                for (int j = 1; j <= powOfTwo; j++) {
                    oldIndexL = 2*j - 2;
                    oldIndexR = 2*j;
                    newIndex  = (int) (0.5*(ptIndex[oldIndexL] + ptIndex[oldIndexR]));

                    tauX      = t[newIndex] - t[ptIndex[oldIndexL]];
                    tauY      = t[ptIndex[oldIndexR]] - t[newIndex];
                    mu        = tauY/tauX;
                    lambdaZ   = delta*delta*tauY*tauY;
                    imu[newIndex]    = mu;
                    ilam[newIndex]   = lambdaZ;
                    imu2[newIndex]   = mu*mu;
                    imuLambdaZ[newIndex]      = mu*lambdaZ;
                    imuOver2LambdaZ[newIndex] = mu/2./lambdaZ;

                    ptIndex[oldIndexL + 1]       = newIndex;
                    wIndexList[indexCounter]   = ptIndex[oldIndexL];
                    wIndexList[indexCounter+1] = newIndex;
                    wIndexList[indexCounter+2] = ptIndex[oldIndexR];

                    indexCounter += 3;
                }
            }

            /* Check if there are holes remaining and fill them */
            for (int k = 1; k < d; k++) {
                if (ptIndex[k-1] + 1 < ptIndex[k]) {
                // there is a hole between (k-1) and k.

                    tauX      = t[ptIndex[k-1]+1] - t[ptIndex[k-1]];
                    tauY      = t[ptIndex[k]]     - t[ptIndex[k-1]+1];
                    mu        = tauY/tauX;
                    lambdaZ   = delta*delta*tauY*tauY;
                    imu[ptIndex[k-1]+1]    = mu;
                    ilam[ptIndex[k-1]+1]   = lambdaZ;
                    imu2[ptIndex[k-1]+1]   = mu*mu;
                    imuLambdaZ[ptIndex[k-1]+1]      = mu*lambdaZ;
                    imuOver2LambdaZ[ptIndex[k-1]+1] = mu/2./lambdaZ;

                    wIndexList[indexCounter]   = ptIndex[k]-2;
                    wIndexList[indexCounter+1] = ptIndex[k]-1;
                    wIndexList[indexCounter+2] = ptIndex[k];
                    indexCounter += 3;
                }
            }
        }
    }

/**
 * Only returns a stream if both inner streams are the same.
 */
public RandomStream getStream() {
        if( stream != otherStream)
            throw new IllegalStateException("Two different streams or more are present");
        return stream;
    }

   /**
    * Sets the streams.
    */
   public void setStream (RandomStream stream, RandomStream otherStream) {
        this.stream = stream;
        this.otherStream = otherStream;
        normalGen.setStream(stream);
    }

   /**
    * Sets both inner streams to the same `stream`.
    */
   public void setStream (RandomStream stream) {
        setStream(stream, stream);
    }

}