/*
 * Class:        JohnsonSystem
 * Description:  Johnson system of distributions
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        july 2012

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