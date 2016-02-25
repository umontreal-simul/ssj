/*
 * Class:        ExponentialGen
 * Description:  random variate generators for the exponential distribution
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
 * This class implements random variate generators for the *exponential*
 * distribution. The density is
 * @f[
 *   f(x) = \lambda e^{-\lambda x} \qquad\mbox{ for }x\ge0,
 * @f]
 * where @f$\lambda> 0@f$.
 *
 * The (non-static) `nextDouble` method simply calls `inverseF` on the
 * distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class ExponentialGen extends RandomVariateGen {
   protected double lambda;

   /**
    * Creates an exponential random variate generator with parameter
    * @f$\lambda@f$ = `lambda`, using stream `s`.
    */
   public ExponentialGen (RandomStream s, double lambda) {
      super (s, new ExponentialDist(lambda));
      setParams (lambda);
   }

   /**
    * Creates a new generator for the exponential distribution `dist` and
    * stream `s`.
    */
   public ExponentialGen (RandomStream s, ExponentialDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getLambda());
   }

   /**
    * Uses inversion to generate a new exponential variate with parameter
    * @f$\lambda= @f$&nbsp;`lambda`, using stream `s`.
    */
   public static double nextDouble (RandomStream s, double lambda) {
      return ExponentialDist.inverseF (lambda, s.nextDouble());
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