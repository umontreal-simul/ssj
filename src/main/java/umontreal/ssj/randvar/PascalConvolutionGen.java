/*
 * Class:        PascalConvolutionGen
 * Description:  Pascal random variate generators using the convolution method
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
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;

/**
 * Implements *Pascal* random variate generators by the *convolution* method
 * (see @cite sLAW00a&thinsp;). The method generates @f$n@f$ geometric
 * variates with probability @f$p@f$ and adds them up.
 *
 * The algorithm is slow if @f$n@f$ is large.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_discrete
 */
public class PascalConvolutionGen extends PascalGen {

   /**
    * Creates a *Pascal* random variate generator with parameters @f$n@f$
    * and @f$p@f$, using stream `s`.
    */
   public PascalConvolutionGen (RandomStream s, int n, double p) {
      super (s, null);
      setParams (n, p);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public PascalConvolutionGen (RandomStream s, PascalDist dist) {
      super (s, dist);
   }
 
    
   public int nextInt() {
      int x = 0;
      for (int i = 0; i < n; i++)
         x += GeometricDist.inverseF (p, stream.nextDouble());
      return x;

   }

/**
 * Generates a new variate from the *Pascal* distribution, with parameters
 * @f$n = @f$&nbsp;`n` and @f$p = @f$&nbsp;`p`, using the stream `s`.
 */
public static int nextInt (RandomStream s, int n, double p) {
     return convolution (s, n, p);
   }
   private static int convolution (RandomStream stream, int n, double p) {
      int x = 0;
      for (int i = 0; i < n; i++)
         x += GeometricDist.inverseF (p, stream.nextDouble());
      return x;
   }
}