/*
 * Class:        BinomialConvolutionGen
 * Description:  binomial random variate generators using the convolution method
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
import umontreal.ssj.probdist.*;
import umontreal.ssj.rng.*;

/**
 * Implements binomial random variate generators using the convolution
 * method. This method generates @f$n@f$ Bernouilli random variates with
 * parameter @f$p@f$ and adds them up. Its advantages are that it requires
 * little computer memory and no setup time. Its disadvantage is that it is
 * very slow for large @f$n@f$. It makes sense only when @f$n@f$ is small.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_discrete
 */
public class BinomialConvolutionGen extends BinomialGen {

   /**
    * Creates a *binomial* random variate generator with parameters
    * @f$n@f$ and @f$p@f$, using stream `s`.
    */
   public BinomialConvolutionGen (RandomStream s, int n, double p) {
      super (s, null);
      setParams (n, p);
   }

   /**
    * Creates a random variate generator for the *binomial* distribution
    * `dist` and stream `s`.
    */
   public BinomialConvolutionGen (RandomStream s, BinomialDist dist) {
      super (s, dist);
   }


   public int nextInt() { 
      int x = 0;
      for (int i = 0; i < n; i++) {
         double unif = stream.nextDouble();
         if (unif <= p)
            x++;
      }
      return x;
   }

/**
 * Generates a new integer from the binomial distribution with parameters
 * @f$n = @f$&nbsp;`n` and @f$p = @f$&nbsp;`p`, using the given stream `s`.
 */
public static int nextInt (RandomStream s, int n, double p) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      if (p < 0 || p > 1)
         throw new IllegalArgumentException ("p must be in [0,1]");
      return convolution (s, n, p);
   }
   private static int convolution (RandomStream stream, int n, double p) {
      int x = 0;
      for (int i = 0; i < n; i++) {
         double unif = stream.nextDouble();
         if (unif <= p)
            x++;
      }
      return x;
   }
}