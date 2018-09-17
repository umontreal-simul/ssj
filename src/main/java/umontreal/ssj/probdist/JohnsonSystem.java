/*
 * Class:        JohnsonSystem
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
package umontreal.ssj.probdist;

/**
 * This class contains common parameters and methods for the *Johnson* system
 * of distributions @cite tJOH49a, @cite tJOH95a&thinsp; with shape
 * parameters @f$\gamma@f$ and @f$\delta> 0@f$, location parameter
 * @f$\xi@f$, and scale parameter @f$\lambda>0@f$. Denoting
 * @f$T=(X-\xi)/\lambda@f$, the variable @f$Z = \gamma+ \delta f(T)@f$ is
 * a standard normal variable, where @f$f(t)@f$ is one of the following
 * transformations:
 *
 * <center>
 *
 * <table class="SSJ-table SSJ-has-hlines">
 * <tr class="bt">
 *   <td class="c bl br">Family</td>
 *   <td class="c bl br">@f$f(t)@f$</td>
 * </tr><tr class="bt">
 *   <td class="c bl br">@f$S_L@f$</td>
 *   <td class="c bl br">@f$\ln(t)@f$</td>
 * </tr><tr>
 *   <td class="c bl br">@f$S_B@f$</td>
 *   <td class="c bl br">@f$\ln(t / (1-t))@f$</td>
 * </tr><tr>
 *   <td class="c bl br">@f$S_U@f$</td>
 *   <td class="c bl br">@f$\ln(t + \sqrt{1 + t^2})@f$</td>
 * </tr>
 * </table>
 *
 * </center>
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
abstract class JohnsonSystem extends ContinuousDistribution {
   protected double gamma;
   protected double delta;
   protected double xi;
   protected double lambda;

   /**
    * Constructs a `JohnsonSystem` object with shape parameters
    * @f$\gamma= \mathtt{gamma}@f$ and @f$\delta= \mathtt{delta}@f$,
    * location parameter @f$\xi= \mathtt{xi}@f$, and scale parameter
    * @f$\lambda= \mathtt{lambda}@f$.
    */
   protected JohnsonSystem (double gamma, double delta, double xi,
                            double lambda) {
      setParams0 (gamma, delta, xi, lambda);
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
    * @f$\xi@f$ and @f$\lambda@f$.
    */
   protected void setParams0(double gamma, double delta, double xi,
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

   /**
    * Return an array containing the parameters of the current
    * distribution. This array is put in regular order: [@f$\gamma@f$,
    * @f$\delta@f$, @f$\xi@f$, @f$\lambda@f$].
    */
   public double[] getParams () {
      double[] retour = {gamma, delta, xi, lambda};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : gamma = " + gamma + ", delta = " + delta + ", xi = " + xi + ", lambda = " + lambda;
   }

}