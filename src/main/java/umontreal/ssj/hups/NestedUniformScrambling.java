/*
 * Class:        NestedUniformScrambling
 * Description:  performs Owen's nested uniform scrambling
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2016  David Munger, Pierre L'Ecuyer and Universite de Montreal
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
 * This class implements a PointSetRandomization that performs Owen's nested
 * uniform scrambling \cite vOWE95a, \cite vOWE03a .
 * The point set must be a CachedPointSet of a DigitalNetBase2.
 *
 * To use this randomization, one should should call the randomize() function of
 * an instance of NestedUniformScrambling with, as its argument, an instance of
 * CachedPointSet with a reference to a DigitalNetBase2 instance.
 *
 * @warning Calling CachedPointSet.randomize() with an instance of
 * NestedUniformScrambling as its arguments will not work, because
 * CachedPointSet.randomize() calls randomize() on its reference point set.
 * However, this randomization modifies the cached values and not the original
 * point set.
 */
public class NestedUniformScrambling implements PointSetRandomization {

   private RandomStream stream;

   /**
    * Empty constructor.
    */
   public NestedUniformScrambling() {
      this(null);
   }

   /**
    * Sets internal variable `stream` to the given `stream`.
    *  @param stream       stream to use in the randomization
    */
   public NestedUniformScrambling (RandomStream stream) {
       this.stream = stream;
   }

   @Override public RandomStream getStream() {
      return stream;
   }

   @Override public void setStream(RandomStream stream) {
      this.stream = stream;
   }

   /**
    * This method calls scrambles the points in the cached point set.
    * If `p` is not a CachedPointSet of a DigitalNetBase2, an IllegalArgumentException
    * is thrown.
    *  @param p            Point set to randomize
    */
   public void randomize (PointSet p) {
      if (p instanceof CachedPointSet) {
         CachedPointSet cp = (CachedPointSet) p;
         if (cp.getRefPointSet() instanceof DigitalNetBase2) {
            ((DigitalNetBase2) cp.getRefPointSet()).nestedUniformScramble(stream, cp.getArray());
            return;
         }
      }
      throw new IllegalArgumentException("NestedUniformScrambling" +
         " can only randomize a CachedPointSet of a DigitalNetBase2");
   }
}
