/*
 * Class:        RandomShift
 * Description:  Applies a random shift on a point set
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

/**
 * This class implements a  @ref umontreal.ssj.hups.PointSetRandomization.
 * The  @ref umontreal.ssj.rng.RandomStream is stored internally. The method
 * #randomize(PointSet) simply calls
 * {@link umontreal.ssj.hups.PointSet.addRandomShift(RandomStream)
 * addRandomShift(stream)}.
 *
 * This class can be used as a base class to implement a specific
 * randomization by overriding method  #randomize(PointSet).
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class RandomShift implements PointSetRandomization {
   protected RandomStream stream;

   /**
    * Empty constructor.
    */
   public RandomShift() {
   }

   /**
    * Sets the internal  @ref umontreal.ssj.rng.RandomStream to `stream`.
    *  @param stream       stream to use in the randomization
    */
   public RandomShift (RandomStream stream) {
       this.stream = stream;
   }

   /**
    * This method calls
    * {@link umontreal.ssj.hups.PointSet.addRandomShift(RandomStream)
    * addRandomShift(stream)}.
    *  @param p            Point set to randomize
    */
   public void randomize (PointSet p) {
      p.addRandomShift(stream);
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