/*
 * Class:        JohnsonSystemG
 * Description:  Johnson system of distributions
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        july 2012
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
 * This class contains common parameters and methods for the random variate
 * generators associated with the *Johnson* system of distributions
 * @cite tJOH49a, @cite tJOH95a&thinsp;. See the definitions of
 * @ref umontreal.ssj.probdist.JohnsonSLDist,
 * @ref umontreal.ssj.probdist.JohnsonSBDist, and
 * @ref umontreal.ssj.probdist.JohnsonSUDist in package `probdist`.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
abstract class JohnsonSystemG extends RandomVariateGen {
   protected double gamma;
   protected double delta;
   protected double xi;
   protected double lambda;

   /**
    * Constructs a `JohnsonSystemG` object with shape parameters
    * @f$\gamma= \mathtt{gamma}@f$ and @f$\delta= \mathtt{delta}@f$,
    * location parameter @f$\xi= \mathtt{xi}@f$, and scale parameter
    * @f$\lambda= \mathtt{lambda}@f$.
    */
   protected JohnsonSystemG (RandomStream s, double gamma, double delta,
                             double xi, double lambda) {
      super (s, null);
      setParams (gamma, delta, xi, lambda);
   }

   /**
    * Constructs a `JohnsonSystemG` object with parameters obtained from
    * distribution `dist`.
    */
   protected JohnsonSystemG (RandomStream s, ContinuousDistribution dist) {
      super (s, dist);
   }

   /**
    * Returns the value of @f$\gamma@f$.
    */
   public double getGamma() {
      return gamma;
   }

   /**
    * Returns the value of @f$\delta@f$.
    */
   public double getDelta() {
      return delta;
   }

   /**
    * Returns the value of @f$\xi@f$.
    */
   public double getXi() {
      return xi;
   }

   /**
    * Returns the value of @f$\lambda@f$.
    */
   public double getLambda() {
      return lambda;
   }

   /**
    * Sets the value of the parameters @f$\gamma@f$, @f$\delta@f$,
    * @f$\xi@f$ and @f$\lambda@f$ for this object.
    */
   protected void setParams (double gamma, double delta, double xi,
                             double lambda) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (delta <= 0)
         throw new IllegalArgumentException ("delta <= 0");
      this.gamma = gamma;
      this.delta = delta;
      this.xi = xi;
      this.lambda = lambda;
   }

}