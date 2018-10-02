/*
 * Class:        RandomStart
 * Description:  Randomizes a sequence with a random starting point
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
package umontreal.ssj.hups;
 import umontreal.ssj.rng.RandomStream;
 import java.lang.IllegalArgumentException;

/**
 * This class implements a  @ref umontreal.ssj.hups.PointSetRandomization
 * that randomizes a sequence simply by taking a random starting point. 
 * For now, this only applies to the Halton sequence, but it could eventually
 * be generalized to other types of sequences.  The point set
 * must be an instance of  @ref umontreal.ssj.hups.HaltonSequence or an
 * IllegalArgumentException is thrown.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class RandomStart implements PointSetRandomization {

   protected RandomStream stream;

   /**
    * Empty constructor.
    */
   public RandomStart() {
   }

   /**
    * Sets internal variable `stream` to the given `stream`.
    *  @param stream       stream to use in the randomization
    */
   public RandomStart (RandomStream stream) {
       this.stream = stream;
   }

   /**
    * This method calls  umontreal.ssj.hups.HaltonSequence.init(double[]).
    * If `p` is not a  @ref umontreal.ssj.hups.HaltonSequence, an
    * IllegalArgumentException is thrown.
    *  @param p            Point set to randomize
    */
   public void randomize (PointSet p) {
      if (p instanceof HaltonSequence) {
         double[] x0 = new double[p.getDimension()];
         stream.nextArrayOfDouble(x0, 0, x0.length);
         ((HaltonSequence)p).setStart (x0);
      } else {
         throw new IllegalArgumentException("RandomStart" +
                     " can only randomize a HaltonSequence");
      }
   }

   /**
    * Sets the internal  @ref umontreal.ssj.rng.RandomStream to `stream`.
    *  @param stream       stream to use in the randomization
    */
   public void setStream (RandomStream stream) {
      this.stream = stream;
   }

   /**
    * Returns the internal  @ref umontreal.ssj.rng.RandomStream.
    *  @return stream used in the randomization
    */
   public RandomStream getStream() {
      return stream;
   }

}