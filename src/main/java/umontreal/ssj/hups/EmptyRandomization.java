/*
 * Class:        EmptyRandomization
 * Description:  implements an empty PointSetRandomization
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
 import umontreal.ssj.rng.MRG32k3a;
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
   protected RandomStream stream = new MRG32k3a();

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