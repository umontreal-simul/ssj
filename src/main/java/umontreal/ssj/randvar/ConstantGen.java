/*
 * Class:        ConstantGen
 * Description:  random variate generator for a constant distribution
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