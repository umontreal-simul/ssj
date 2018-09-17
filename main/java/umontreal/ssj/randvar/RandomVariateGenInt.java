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
    * Generates `n` random numbers from the discrete distribution
    * contained in this object, and returns them in a new array of size `n`.
    * By default, this method calls
    * #nextInt() `n` times, but one can override it in subclasses for
    * better efficiency.
    *  @param n            number of variates to generate
    *  @return  a new array with the generated numbers
    */
   public int[] nextArrayOfInt (int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n must be positive.");
      int[] v = new int[n];
      for (int i = 0; i < n; i++)
         v[i] = nextInt();
      return v;
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