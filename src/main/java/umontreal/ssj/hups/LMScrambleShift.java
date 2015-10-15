/*
 * Class:        LMScrambleShift
 * Description:  performs a left matrix scramble and adds a random digital shift
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
 * that performs a left matrix scrambling and adds a random digital shift.
 * Point set must be a  @ref umontreal.ssj.hups.DigitalNet or an
 * IllegalArgumentException is thrown.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class LMScrambleShift extends RandomShift {

   /**
    * Empty constructor.
    */
   public LMScrambleShift() {
   }

   /**
    * Sets internal variable `stream` to the given `stream`.
    *  @param stream       stream to use in the randomization
    */
   public LMScrambleShift (RandomStream stream) {
       super(stream);
   }

   /**
    * This method calls
    * umontreal.ssj.hups.DigitalNet.leftMatrixScramble(RandomStream), then
    * umontreal.ssj.hups.DigitalNet.addRandomShift(RandomStream). If `p`
    * is not a  @ref umontreal.ssj.hups.DigitalNet, an
    * IllegalArgumentException is thrown.
    *  @param p            Point set to randomize
    */
   public void randomize (PointSet p) {
      if (p instanceof DigitalNet) {
         ((DigitalNet)p).leftMatrixScramble (stream);
         ((DigitalNet)p).addRandomShift (stream);
      } else {
         throw new IllegalArgumentException("LMScrambleShift"+
                                            " can only randomize a DigitalNet");
      }
   }

}