/*
 * Interface:    PointSetRandomization
 * Description:  Used to randomize a PointSet
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001--2018  Pierre L'Ecuyer and Universite de Montreal
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
 * This interface is for a *randomization* that can be used to 
 * randomize a  @ref umontreal.ssj.hups.PointSet.
 * One can implement method  #randomize(PointSet) in any way. This method
 * must use an internal  @ref umontreal.ssj.rng.RandomStream. This stream can
 * be set in the constructor, but the methods  #getStream and
 * #setStream(RandomStream) must also be implemented.
 *
 * The method  #randomize(PointSet) can be implemented using combinations of
 * the randomization methods from the point set such as
 * umontreal.ssj.hups.PointSet.addRandomShift,
 * umontreal.ssj.hups.DigitalNet.leftMatrixScramble,
 * umontreal.ssj.hups.DigitalNet.stripedMatrixScramble, etc.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public interface PointSetRandomization {

   /**
    * This method must randomize `p`.
    *  @param p            Point set to randomize
    */
   public void randomize (PointSet p);

   /**
    * Sets the internal  @ref umontreal.ssj.rng.RandomStream to `stream`.
    *  @param stream       stream to use in the randomization
    */
   public void setStream (RandomStream stream);

   /**
    * Returns the internal  @ref umontreal.ssj.rng.RandomStream.
    *  @return stream used in the randomization
    */
   public RandomStream getStream();

}