/*
 * Class:        ChiSquareGen
 * Description:  random variate generators with the chi square distribution
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
 * This class implements random variate generators with the *chi square*
 * distribution with @f$n>0@f$ degrees of freedom. Its density function is
 * @anchor REF_randvar_ChiSquareGen_eq_Fchi2
 * @f[
 *   f(x) = \frac{e^{-x/2}x^{n/2-1}}{2^{n/2}\Gamma(n/2)} \qquad\mbox{ for } x > 0, \tag{Fchi2}
 * @f]
 * where @f$\Gamma(x)@f$ is the gamma function defined in (
 * {@link REF_randvar_GammaGen_eq_Gamma Gamma} ).
 *
 * The (non-static) `nextDouble` method simply calls `inverseF` on the
 * distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class ChiSquareGen extends RandomVariateGen {
   protected int n = -1;

   /**
    * Creates a *chi square* random variate generator with @f$n@f$ degrees
    * of freedom, using stream `s`.
    */
   public ChiSquareGen (RandomStream s, int n) {
      super (s, new ChiSquareDist(n));
      setParams (n);
   }

   /**
    * Create a new generator for the distribution `dist` and stream `s`.
    */
   public ChiSquareGen (RandomStream s, ChiSquareDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getN ());
   }

   /**
    * Generates a new variate from the chi square distribution with
    * @f$n@f$ degrees of freedom, using stream `s`.
    */
   public static double nextDouble (RandomStream s, int n) {
      return ChiSquareDist.inverseF (n, s.nextDouble());
   }

   /**
    * Returns the value of @f$n@f$ for this object.
    */
   public int getN() {
      return n;
   }
   protected void setParams (int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      this.n = n;
   }
}