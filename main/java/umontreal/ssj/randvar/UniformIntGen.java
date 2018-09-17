/*
 * Class:        UniformIntGen
 * Description:  random variate generator for the uniform distribution over integers
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
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;

/**
 * This class implements a random variate generator for the *uniform*
 * distribution over integers, over the interval @f$[i,j]@f$. Its mass
 * function is
 * @anchor REF_randvar_UniformIntGen_eq_fmassuniformint
 * @f[
 *   p(x) = \frac{1}{j - i + 1} \qquad\mbox{ for } x = i, i + 1, …, j \tag{fmassuniformint}
 * @f]
 * and 0 elsewhere.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_discrete
 */
public class UniformIntGen extends RandomVariateGenInt {
   protected int left;     // the left limit of the interval
   protected int right;    // the right limit of the interval

   /**
    * Creates a uniform random variate generator over the integers in the
    * closed interval @f$[i, j]@f$, using stream `s`.
    */
   public UniformIntGen (RandomStream s, int i, int j) {
      super (s, new UniformIntDist (i, j));
      setParams (i, j);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public UniformIntGen (RandomStream s, UniformIntDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getI(), dist.getJ());
   }

   /**
    * Generates a new *uniform* random variate over the interval
    * @f$[i,j]@f$, using stream `s`, by inversion.
    */
   public static int nextInt (RandomStream s, int i, int j) {
      return UniformIntDist.inverseF (i, j, s.nextDouble());
   }

   /**
    * Returns the parameter @f$i@f$.
    */
   public int getI() {
      return left;
   }

   /**
    * Returns the parameter @f$j@f$.
    */
   public int getJ() {
      return right;
   }
   protected  void setParams (int i, int j) {
      if (j < i)
        throw new IllegalArgumentException ("j < i");
      left = i;
      right = j;
   }
}