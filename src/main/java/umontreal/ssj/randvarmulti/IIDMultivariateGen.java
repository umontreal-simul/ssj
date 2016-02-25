/*
 * Class:        IIDMultivariateGen
 * Description:  vector of independent identically distributed random variables
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       David Munger
 * @since        January 2011
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
package umontreal.ssj.randvarmulti;
import umontreal.ssj.randvar.RandomVariateGen;

/**
 * Extends  @ref RandomMultivariateGen for a vector of independent
 * identically distributed (i.i.d.) random variables.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class IIDMultivariateGen extends RandomMultivariateGen {

   /**
    * Constructs a generator for a <tt>d</tt>-dimensional vector of i.i.d.
    * variates with a common one-dimensional generator `gen1`.
    *  @param gen1         the one-dimensional generator
    *  @param d            dimension of the vector (number of i.i.d.
    *                      variates).
    */
   public IIDMultivariateGen (RandomVariateGen gen1, int d) {
      setGen1 (gen1);
      this.stream = gen1.getStream();
      dimension = d;
   }

   /**
    * Changes the dimension of the vector to `d`.
    */
   public void setDimension (int d) {
      dimension = d;
   }

   /**
    * Generates a vector of i.i.d. variates.
    */
   public void nextPoint (double[] p) {
      if (p.length != dimension)
         throw new IllegalArgumentException(String.format(
            "p's dimension (%d) does not mach dimension (%d)", p.length, dimension));

      for (int i = 0; i < dimension; i++)
         p[i] = gen1.nextDouble();
   }

   /**
    * Sets the common one-dimensional generator to `gen1`.
    */
   public void setGen1 (RandomVariateGen gen1) {
      if (gen1 == null)
         throw new NullPointerException ("gen1 is null");
      this.gen1 = gen1;
   }

   /**
    * Returns the common one-dimensional generator used in this class.
    */
   public RandomVariateGen getGen1() {
     return gen1;
   }

   /**
    * Returns a string representation of the generator.
    */
   public String toString() {
      return dimension + "-dimensional vector of i.i.d. " +
            gen1.getDistribution().toString();
   }

}