/*
 * Interface:    PointSetRandomization
 * Description:  Used to randomize a PointSet
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
package umontreal.ssj.hups;
import umontreal.ssj.rng.RandomStream;

/**
 * This interface is used to randomize a  @ref umontreal.ssj.hups.PointSet.
 * One can implement method  #randomize(PointSet) in any way. This method
 * must use an internal  @ref umontreal.ssj.rng.RandomStream. This stream can
 * be set in the constructor, but the methods  #getStream and
 * #setStream(RandomStream) must be implemented.
 *
 * The method  #randomize(PointSet) must be implemented using combinations of
 * the randomization methods from the point set such as
 * umontreal.ssj.hups.PointSet.addRandomShift,
 * umontreal.ssj.hups.DigitalNet.leftMatrixScramble,
 * umontreal.ssj.hups.DigitalNet.stripedMatrixScramble, …
 *
 * If more than one  @ref PointSetRandomization is applied to the same point
 * set, the randomizations will concatenate if they are of different types,
 * but only the last of each type will remain.
 *
 * @remark **Pierre:** There should be examples to illustrate how this works
 * and how to use it.
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