/*
 * Class:        LognormalSpecialGen
 * Description:  random variates from the lognormal distribution
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

/**
 * Implements methods for generating random variates from the *lognormal*
 * distribution using an arbitrary normal random variate generator. The
 * (non-static) `nextDouble` method calls the `nextDouble` method of the
 * normal generator and takes the exponential of the result.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class LognormalSpecialGen extends RandomVariateGen {

   NormalGen myGen;

   /**
    * Create a lognormal random variate generator using the normal
    * generator `g` and with the same parameters.
    */
   public LognormalSpecialGen (NormalGen g) {
      // Necessary to compile, but we do not want to use stream and dist
      super (g.stream, null);
      stream = null;
      myGen = g;
   }
 

   public double nextDouble() {
      return Math.exp (myGen.nextDouble());
   }
}