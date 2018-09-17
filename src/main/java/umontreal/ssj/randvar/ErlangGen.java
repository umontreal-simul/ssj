/*
 * Class:        ErlangGen
 * Description:  random variate generators for the Erlang distribution
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
 * This class implements random variate generators for the *Erlang*
 * distribution with parameters @f$k > 0@f$ and @f$\lambda> 0@f$. This
 * Erlang random variable is the sum of @f$k@f$ exponentials with parameter
 * @f$\lambda@f$ and has mean @f$k/\lambda@f$.
 *
 * The (non-static) `nextDouble` method simply calls `inverseF` on the
 * distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class ErlangGen extends GammaGen {
   protected int    k = -1;

   /**
    * Creates an Erlang random variate generator with parameters `k` and
    * @f$\lambda@f$ = `lambda`, using stream `s`.
    */
   public ErlangGen (RandomStream s, int k, double lambda) {
      super (s, new ErlangDist(k, lambda));
      setParams (k, lambda);
   }

   /**
    * Creates an Erlang random variate generator with parameters `k` and
    * @f$\lambda= 1@f$, using stream `s`.
    */
   public ErlangGen (RandomStream s, int k) {
      this (s, k, 1.0);
   }

   /**
    * Creates a new generator for the distribution `dist` and stream `s`.
    */
   public ErlangGen (RandomStream s, ErlangDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getK(), dist.getLambda());
   }

   /**
    * Generates a new variate from the Erlang distribution with parameters
    * @f$k = @f$&nbsp;`k` and @f$\lambda= @f$&nbsp;`lambda`, using stream
    * `s`.
    */
   public static double nextDouble (RandomStream s, int k, double lambda) {
      return ErlangDist.inverseF (k, lambda, 15, s.nextDouble());
   }

   /**
    * Returns the parameter @f$k@f$ of this object.
    */
   public int getK() {
      return k;
   }

   /**
    * Sets the parameter @f$k@f$ and @f$\lambda@f$ of this object.
    */
   protected void setParams (int k, double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (k <= 0)
         throw new IllegalArgumentException ("k <= 0");
      this.lambda = lambda;
      this.k = k;
   }
}