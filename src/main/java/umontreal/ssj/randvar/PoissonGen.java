/*
 * Class:        PoissonGen
 * Description:  random variate generators having the Poisson distribution
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
 * This class implements random variate generators having the *Poisson*
 * distribution. Its mass function is
 * @anchor REF_randvar_PoissonGen_eq_fmass_Poisson
 * @f[
 *   p(x) = \frac{e^{-\lambda} \lambda^x}{x!} \qquad\mbox{for } x=0,1,…, \tag{fmass-Poisson}
 * @f]
 * where @f$\lambda> 0@f$ is a real valued parameter equal to the mean.
 *
 * No local copy of the parameter @f$\lambda= @f$ `lambda` is maintained in
 * this class. The (non-static) `nextInt` method simply calls `inverseF` on
 * the distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_discrete
 */
public class PoissonGen extends RandomVariateGenInt {
   protected double lambda;

   /**
    * Creates a Poisson random variate generator with parameter
    * @f$\lambda= @f$ `lambda`, using stream `s`.
    */
   public PoissonGen (RandomStream s, double lambda) {
      super (s, new PoissonDist (lambda));
      setParams (lambda);
   }

   /**
    * Creates a new random variate generator using the Poisson
    * distribution `dist` and stream `s`.
    */
   public PoissonGen (RandomStream s, PoissonDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getLambda());
   }

   /**
    * A static method for generating a random variate from a *Poisson*
    * distribution with parameter @f$\lambda@f$ = `lambda`.
    */
   public static int nextInt (RandomStream s, double lambda) {
      return PoissonDist.inverseF (lambda, s.nextDouble());
   }

   /**
    * Returns the @f$\lambda@f$ associated with this object.
    */
   public double getLambda() {
      return lambda;
   }

   /**
    * Sets the parameter @f$\lambda= @f$ `lam` of this object.
    */
   protected void setParams (double lam) {
      if (lam <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      this.lambda = lam;
   }
}