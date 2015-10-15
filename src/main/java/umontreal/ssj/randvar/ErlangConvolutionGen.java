/*
 * Class:        ErlangConvolutionGen
 * Description:  Erlang random variate generators using the convolution method
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
 * This class implements *Erlang* random variate generators using the
 * *convolution* method. This method uses inversion to generate @f$k@f$
 * exponential variates with parameter @f$\lambda@f$ and returns their sum.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class ErlangConvolutionGen extends ErlangGen {

   /**
    * Creates an Erlang random variate generator with parameters `k` and
    * @f$\lambda@f$ = `lambda`, using stream `s`.
    */
   public ErlangConvolutionGen (RandomStream s, int k, double lambda) {
      super (s, null);
      setParams (k, lambda);
   }

   /**
    * Creates an Erlang random variate generator with parameters `k` and
    * @f$\lambda= 1@f$, using stream `s`.
    */
   public ErlangConvolutionGen (RandomStream s, int k) {
      this (s, k, 1.0);
   }

   /**
    * Creates a new generator for the distribution `dist` and stream `s`.
    */
   public ErlangConvolutionGen (RandomStream s, ErlangDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getK(), dist.getLambda());
   }


   public double nextDouble() {
      return convolution (stream, k, lambda);
   }

   public static double nextDouble (RandomStream s, int k, double lambda) {
      return convolution (s, k, lambda);
   }

   private static double convolution (RandomStream s, int k, double lambda) {
      double x = 0.0;
      for (int i=0;  i<k;  i++)  
         x += ExponentialDist.inverseF (lambda, s.nextDouble());
      return x;
   }
}