/*
 * Class:        ConstantIntDist
 * Description:  constant integer distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Éric Buist
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
package umontreal.ssj.probdist;

/**
 * Represents a *constant* discrete distribution taking a single integer
 * value with probability 1. Its mass function is
 * @anchor REF_probdist_ConstantIntDist_eq_fconsint
 * @f[
 *   p(x) = \left\{\begin{array}{ll}
 *    1, 
 *    & 
 *    \qquad\mbox{for } x = c,
 *    \\ 
 *    0, 
 *    & 
 *    \qquad\mbox{elsewhere. } 
 *   \end{array}\right. \tag{fconsint}
 * @f]
 * Its distribution function is
 * @anchor REF_probdist_ConstantIntDist_eq_cdfconsint
 * @f[
 *   F(x) = \left\{\begin{array}{ll}
 *    0, 
 *    & 
 *    \qquad\mbox{ for } x < c
 *    \\ 
 *    1, 
 *    & 
 *    \qquad\mbox{ for } x \ge c. 
 *   \end{array}\right. \tag{cdfconsint}
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_discrete
 */
public class ConstantIntDist extends UniformIntDist {

   /**
    * Constructs a new constant distribution with probability 1 at `c`.
    */
   public ConstantIntDist (int c) {
      super (c, c);
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : c = " + i;
   }

}