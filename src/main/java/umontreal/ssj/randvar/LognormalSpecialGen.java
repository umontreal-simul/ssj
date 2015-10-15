/*
 * Class:        LognormalSpecialGen
 * Description:  random variates from the lognormal distribution
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