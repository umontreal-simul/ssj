/*
 * Class:        EmptyRandomization
 * Description:  implements an empty PointSetRandomization
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
 // import umontreal.ssj.rng.MRG32k3a;
 import umontreal.ssj.rng.RandomStream;

/**
 * This class implements an empty
 * @ref umontreal.ssj.hups.PointSetRandomization. The method
 * #randomize(PointSet) does nothing. The internal stream is never used. This
 * class can be used in methods where a randomization is needed but you don’t
 * want one.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class EmptyRandomization implements PointSetRandomization {
   protected RandomStream stream;   //  = new MRG32k3a();

   /**
    * This method does nothing.
    *  @param p            Point set to randomize
    */
   public void randomize (PointSet p) {
      // Does nothing
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