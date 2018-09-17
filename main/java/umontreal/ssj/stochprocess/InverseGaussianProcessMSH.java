/*
 * Class:        InverseGaussianProcessMSH
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
 * Uses a faster generating method (MSH) @cite rMIC76a&thinsp; than the
 * simple inversion of the distribution function used by
 * @ref InverseGaussianProcess. It is about 60 times faster. However it
 * requires two  @ref umontreal.ssj.rng.RandomStream ’s instead of only one
 * for  @ref InverseGaussianProcess. The second stream is called
 * `otherStream` below and it is used to randomly choose between two roots at
 * each time step.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class InverseGaussianProcessMSH extends InverseGaussianProcess {

    // otherStream is used to decide between the two roots in method MSH.
    protected RandomStream otherStream;
    // Needed for the MSH method of generating inverse gaussian.
    protected NormalGen normalGen;

   /**
    * Constructs a new `InverseGaussianProcessMSH`. The initial value `s0`
    * will be overridden by @f$t[0]@f$ when the observation times are set.
    */
   public InverseGaussianProcessMSH (double s0, double delta, double gamma,
                                     RandomStream stream, 
                                     RandomStream otherStream) {
        super(); // dummy
        this.x0 = s0;
        setParams(delta, gamma);
        this.stream = stream;
        this.otherStream = otherStream;
        normalGen = new NormalGen(stream); 
        numberOfRandomStreams = 2;
    }

   /**
    * Generates the path. It is done by successively calling
    * `nextObservation()`, therefore the two
    * @ref umontreal.ssj.rng.RandomStream s are sampled alternatively.
    */
   public double[] generatePath() {
        double s = x0;
        for (int i = 0; i < d; i++) 
        {
            s += InverseGaussianMSHGen.nextDouble(otherStream,
                                                normalGen, imu[i], ilam[i]);
            path[i+1] = s;
        }
        observationIndex   = d;
        observationCounter = d;
        return path;
    }

   /**
    * Instead of using the internal streams to generate the path, uses two
    * arrays of uniforms @f$U[0,1)@f$. The length of the arrays should be
    * equal to the number of periods in the observation times. This method
    * is useful for  @ref NormalInverseGaussianProcess.
    */
   public double[] generatePath (double[] unifNorm, double[] unifOther) {
        double s = x0;
        // The class NonRandomStream is defined below.
        RandomStream nonRandOther = new NonRandomStream(unifOther);
        // this.stream should keep in memory the original stream of the Normal.
        normalGen.setStream(new NonRandomStream(unifNorm));
        for (int i = 0; i < d; i++) {
            s += InverseGaussianMSHGen.nextDouble(nonRandOther,
                                                normalGen, imu[i], ilam[i]);
            path[i+1] = s;
        }
        observationIndex   = d;
        observationCounter = d;
        normalGen.setStream(stream);  // reset to original stream
        return path;
    }

   /**
    * Not implemented, requires two  @ref umontreal.ssj.rng.RandomStream
    * ’s.
    */
   public double[] generatePath (double[] uniforms01) {
       throw new UnsupportedOperationException("Use generatePath with 2 streams");
    }


   public double nextObservation() {
        double s = path[observationIndex];
        s += InverseGaussianMSHGen.nextDouble(otherStream, normalGen,
                  imu[observationIndex], ilam[observationIndex]);
        observationIndex++;
        observationCounter = observationIndex;
        path[observationIndex] = s;
        return s;
    }

/**
 * Only returns a stream if both inner  @ref umontreal.ssj.rng.RandomStream
 * ’s are the same.
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
        super.setStream(stream);
        normalGen.setStream(stream);
        setOtherStream(otherStream);
    }

   /**
    * Sets both inner streams to `stream`.
    */
   public void setStream (RandomStream stream) {
        super.setStream(stream);
        normalGen.setStream(stream);
        setOtherStream(stream);
    }

   /**
    * Sets the `otherStream`, which is the stream used to choose between
    * the two roots in the MSH method.
    */
   public void setOtherStream (RandomStream otherStream) {
        this.otherStream = otherStream;
    }

   /**
    * Returns the `otherStream`, which is the stream used to choose
    * between the two quadratic roots from the MSH method.
    */
   public RandomStream getOtherStream() {
        return otherStream;
    }

   /**
    * Sets the normal generator. It also sets one of the two inner streams
    * to the stream of the normal generator.
    */
   public void setNormalGen (NormalGen normalGen) {
        this.normalGen = normalGen;
        stream = normalGen.getStream();
    }

   /**
    * Returns the normal generator.
    */
   public NormalGen getNormalGen() {
        return normalGen;
    }


/**
 *   NonRandomStream:     
 * Given a double array, this class will return those values
 * as if it where a random stream.
 * Careful: Will not hard copy the array given as input.
 * And not checking for end of array for the time being.
 * And not checking j>i.
 */
    protected class NonRandomStream implements RandomStream
    {
       double[] array;
       int position;

       public NonRandomStream(double[] array)
       {
	  this.array = array;
	  position = 0;
       }

       public NonRandomStream(double value)
       {
	  this.array = new double[]{value};
	  position = 0;
       }

       public double nextDouble()
       {
          return array[position++];
       }
    
       public void nextArrayOfDouble(double[] u, int start, int n)
       {
	  for(int i = 0; i < n; i++)
	     u[start+i] = array[position++];
       }

       public void nextArrayOfInt(int i, int j, int[] u, 
			         int start, int n)
       {
	  double diff = (double)(j - i);
	  for(int ii = 0; ii < n; ii++)
	      u[start+ii] = i + 
		(int)Math.round(diff * array[position++]);
       }
    
       public int nextInt(int i, int j)
       {
	  return (int)Math.round( (double)(j-i) * array[position]);
       }

    
       public void resetNextSubstream()
       {
       }

       public void resetStartStream()
       {
	  position = 0;
       }

       public void resetStartSubstream()
       {
       }

       public String toString()
       {
	 return new String("NonRandomStream of length " +
		      array.length);
       }
    }

}