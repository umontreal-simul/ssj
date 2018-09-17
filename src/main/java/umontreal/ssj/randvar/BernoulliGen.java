/*
 * Class:        BernoulliGen
 * Description:  random variate generators for the Bernoulli distribution
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
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
import umontreal.ssj.probdist.BernoulliDist;
import umontreal.ssj.rng.RandomStream;

/**
 * This class implements random variate generators for the *Bernoulli*
 * distribution (see class  @ref umontreal.ssj.probdist.BernoulliDist ).
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_discrete
 */
public class BernoulliGen extends RandomVariateGenInt {
   protected double p;

   /**
    * Creates a Bernoulli random variate generator with parameter @f$p@f$,
    * using stream `s`.
    */
   public BernoulliGen (RandomStream s, double p) {
      super (s, new BernoulliDist (p));
      setParams (p);
   }

   /**
    * Creates a random variate generator for the *Bernoulli* distribution
    * `dist` and the random stream `s`.
    */
   public BernoulliGen (RandomStream s, BernoulliDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getP());
   }

   /**
    * Generates a new integer from the *Bernoulli* distribution with
    * parameter @f$p = @f$&nbsp;`p`, using the given stream `s`.
    */
   public static int nextInt (RandomStream s, double p) {
      return BernoulliDist.inverseF (p, s.nextDouble());
   }

   /**
    * Returns the parameter @f$p@f$ of this object.
    */
   public double getP() {
      return p;
   }

   /**
    * Sets the parameter @f$p@f$ of this object.
    */
   protected void setParams (double p) {
      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in range [0, 1]");
      this.p = p;
   }
}