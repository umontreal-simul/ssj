/*
 * Class:        RandomVariateGen
 * Description:  base class for all random variate generators over the reals 
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
package umontreal.ssj.randvar;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.probdist.Distribution;

/**
 * This is the base class for all random variate generators over the real
 * line. It specifies the signature of the  #nextDouble method, which is
 * normally called to generate a real-valued random variate whose
 * distribution has been previously selected. A random variate generator
 * object can be created simply by invoking the constructor of this class
 * with previously created  @ref umontreal.ssj.rng.RandomStream and
 * @ref umontreal.ssj.probdist.Distribution objects, or by invoking the
 * constructor of a subclass. By default, all random variates will be
 * generated via inversion by calling the
 * umontreal.ssj.probdist.Distribution.inverseF method for the distribution,
 * even though this can be inefficient in some cases. For some of the
 * distributions, there are subclasses with special and more efficient
 * methods to generate the random variates.
 *
 * For generating many random variates, creating an object and calling the
 * non-static method is more efficient when the generating algorithm involves
 * a significant setup. When no work is done at setup time, the static
 * methods are usually slightly faster.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_general
 */
public class RandomVariateGen {

   protected RandomStream stream;
   // the stream used for generating random variates

   protected Distribution dist;
   // the distribution used by this generator

 // This constructor is needed for subclasses with no associated distribution.
   protected RandomVariateGen() {}

/**
 * Creates a new random variate generator from the distribution `dist`, using
 * stream `s`.
 *  @param s            random stream used for generating uniforms
 *  @param dist         continuous distribution object of the generated
 *                      values
 */
public RandomVariateGen (RandomStream s, Distribution dist) {
      this.stream = s;
      this.dist   = dist;
   }

   /**
    * Generates a random number from the continuous distribution contained
    * in this object. By default, this method uses inversion by calling
    * the  umontreal.ssj.probdist.ContinuousDistribution.inverseF method
    * of the distribution object. Alternative generating methods are
    * provided in subclasses.
    *  @return the generated value
    */
   public double nextDouble() {
      return dist.inverseF (stream.nextDouble());
   }

   /**
    * Generates `n` random numbers from the continuous distribution
    * contained in this object. These numbers are stored in the array `v`,
    * starting from index `start`. By default, this method calls
    * #nextDouble() `n` times, but one can override it in subclasses for
    * better efficiency.
    *  @param v            array in which the variates will be stored
    *  @param start        starting index, in `v`, of the new variates
    *  @param n            number of variates to generate
    */
   public void nextArrayOfDouble (double[] v, int start, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n must be positive.");
      for (int i = 0; i < n; i++)
         v[start + i] = nextDouble();
   }

   /**
    * Returns the  @ref umontreal.ssj.rng.RandomStream used by this
    * generator.
    *  @return the stream associated to this object
    */
   public RandomStream getStream() { return stream; }

   /**
    * Sets the  @ref umontreal.ssj.rng.RandomStream used by this generator
    * to `stream`.
    */
   public void setStream (RandomStream stream) {
      this.stream = stream;
   }

   /**
    * Returns the  @ref umontreal.ssj.probdist.Distribution used by this
    * generator.
    *  @return the distribution associated to that object
    */
   public Distribution getDistribution() {
      return dist;
   }

   /**
    * Returns a `String` containing information about the current
    * generator.
    */
   public String toString () {
      if (dist != null)
         return getClass().getSimpleName() + " with  " + dist.toString();
      else
         return getClass().getSimpleName() ;
   }

}