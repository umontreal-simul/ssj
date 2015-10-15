/*
 * Class:        ConstantGen
 * Description:  random variate generator for a constant distribution
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
package umontreal.ssj.randvar;

/**
 * This class implements a random variate generator that returns a constant
 * value. Its mass function is
 * @anchor REF_randvar_ConstantGen_eq_randcons
 * @f[
 *   p(x) = \left\{\begin{array}{ll}
 *    1, 
 *    & 
 *    \qquad\mbox{for } x = c,
 *    \\ 
 *    0, 
 *    & 
 *    \qquad\mbox{elsewhere. } 
 *   \end{array}\right. \tag{randcons}
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class ConstantGen extends RandomVariateGen {
   private double val;

   /**
    * Constructs a new constant generator returning the given value `val`.
    */
   public ConstantGen (double val) {
      this.val = val;
   }

   @Override
   public double nextDouble() {
      return val;
   }

   @Override
   public void nextArrayOfDouble (double[] v, int start, int n) {
      for (int i = 0; i < n; i++)
         v[start + i] = val;
   }
}