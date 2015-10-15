/*
 * Class:        ConstantIntDist
 * Description:  constant integer distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Éric Buist
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