/*
 * Class:        RandomVariateGenInt
 * Description:  base class for all generators of discrete random variates
over the integers
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
import umontreal.ssj.probdist.DiscreteDistributionInt;

/**
 * This is the base class for all generators of discrete random variates over
 * the set of integers. Similar to  @ref RandomVariateGen, except that the
 * generators produce integers, via the  #nextInt method, instead of real
 * numbers.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_general
 */
public class RandomVariateGenInt extends RandomVariateGen {

   protected RandomVariateGenInt() {}

/**
 * Creates a new random variate generator for the discrete distribution
 * `dist`, using stream `s`.
 *  @param s            random stream used for generating uniforms
 *  @param dist         discrete distribution object of the generated values
 */
public RandomVariateGenInt (RandomStream s, DiscreteDistributionInt dist) {
      this.stream = s;
      this.dist   = dist;
   }

   /**
    * Generates a random number (an integer) from the discrete
    * distribution contained in this object. By default, this method uses
    * inversion by calling the `inverseF` method of the distribution
    * object. Alternative generating methods are provided in subclasses.
    *  @return the generated value
    */
   public int nextInt() {
      return ((DiscreteDistributionInt) dist).inverseFInt (stream.nextDouble());
   }

   /**
    * Generates `n` random numbers from the discrete distribution
    * contained in this object. The results are stored into the array `v`,
    * starting from index `start`. By default, this method calls
    * #nextInt() `n` times, but one can reimplement it in subclasses for
    * better efficiency.
    *  @param v            array into which the variates will be stored
    *  @param start        starting index, in `v`, of the new variates
    *  @param n            number of variates being generated
    */
   public void nextArrayOfInt (int[] v, int start, int n) {
      if (n < 0)
         throw new IllegalArgumentException ("n must be positive.");
      for (int i = 0; i < n; i++)
         v[start + i] = nextInt();
   }

   /**
    * Returns the  @ref umontreal.ssj.probdist.DiscreteDistributionInt
    * used by this generator.
    *  @return the distribution associated to that object
    */
   public DiscreteDistributionInt getDistribution() {
      return (DiscreteDistributionInt) dist;
   }

}