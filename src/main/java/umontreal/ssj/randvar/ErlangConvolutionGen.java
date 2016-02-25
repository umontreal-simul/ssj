/*
 * Class:        ErlangConvolutionGen
 * Description:  Erlang random variate generators using the convolution method
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