/*
 * Class:        RandomStart
 * Description:  Randomizes a sequence with a random starting point
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
 import java.lang.IllegalArgumentException;

/**
 * This class implements a  @ref umontreal.ssj.hups.PointSetRandomization
 * that randomizes a sequence with a random starting point. The point set
 * must be an instance of  @ref umontreal.ssj.hups.HaltonSequence or an
 * IllegalArgumentException is thrown. For now, only the Halton sequence is
 * allowed, but there may be others later.
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