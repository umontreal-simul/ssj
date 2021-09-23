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
 * This @ref PointSetRandomization class provides the nested uniform scrambling (NUS) randomization proposed 
 * by Owen (\cite vOWE95a, \cite vOWE03a) for digital nets.
 * Since the scrambled points are all stored explicitly, it can only be applied
 * to a @ref CachedPointSet that contains a @ref DigitalNetBase2.
 * The proper way to use it is to construct a @ref CachedPointSet `p` that contains 
 * the digital net, and call @ref NestedUniformScrambling.randomize(p) to randomize.
 * The actual implementation is in @ref DigitalNetBase2.nestedUniformScramble().
 *
 * Note that calling CachedPointSet.randomize() with an instance of
 * NestedUniformScrambling as its arguments will not work, because
 * CachedPointSet.randomize() calls randomize() on its reference point set (the digital net)
 * whereas NUS should modify the cached values instead.
 */
public class NestedUniformScrambling implements PointSetRandomization {

   private RandomStream stream;
   private int numBits;

   /**
    * Empty constructor.
    */
   public NestedUniformScrambling() {
      this(null);
   }

   /**
    * Same as @ref NestedUniformScrambling (stream, 0). 
    *  @param stream       stream to use for the randomization
    */
   public NestedUniformScrambling (RandomStream stream) {
      this(stream, 0);
   }

   /**
    * Create a @ref NestedUniformScrambling instance, using `stream` as the source of
    * randomness, and randomizing only the first `numBits` output bits of
    * the underlying @ref DigitalNetBase2.
    *
    *  @param stream       stream to use in the randomization
    *  @param numBits      number of output bits to scramble (it can be smaller
    *                      than, equal to or larger than the number of output
    *                      bits in the DigitalNetBase2).  If this parameter is zero, 
    *                      `outDigits` bits will be scrambled (up to 31 bits).
    */
   public NestedUniformScrambling (RandomStream stream, int numBits) {
       this.stream = stream;
       this.numBits = numBits;
   }

   public RandomStream getStream() {
      return stream;
   }

   public void setStream (RandomStream stream) {
      this.stream = stream;
   }

   /**
    * Set the number of bits to randomize to `numBits`
    */
   public void setNumBits(int numBits) {
      this.numBits = numBits;
   }

   /**
    * Scrambles the points of the @ref DigitalNetBase2 contained in the @ref CachedPointSet `p`
    * and caches the scrambled points in `p`.  
    * This `p` must be a @ref CachedPointSet of a @ref DigitalNetBase2.
    *  @param p            Point set to randomize
    */
   public void randomize (PointSet p) {
      if (p instanceof CachedPointSet) {
         CachedPointSet cp = (CachedPointSet) p;
         if (cp.getParentPointSet() instanceof DigitalNetBase2) {
            ((DigitalNetBase2) cp.getParentPointSet()).nestedUniformScramble(stream, cp.getArray(), numBits);
            return;
         }
      }
      throw new IllegalArgumentException("NestedUniformScrambling" +
         " can only randomize a CachedPointSet of a DigitalNetBase2");
   }
   
   /**
    * Returns a descriptor of this object.
    */
   public String toString () {
		return "Nested uniform scrambling";
	}
}
