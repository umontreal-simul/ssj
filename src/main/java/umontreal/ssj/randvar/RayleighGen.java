/*
 * Class:        RayleighGen
 * Description:  random variate generators for the Rayleigh distribution
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since

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
package umontreal.ssj.randvar;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;

/**
 * This class implements random variate generators for the *Rayleigh*
 * distribution. Its density is
 * @f[
 *   f(x) = \left\{ \begin{array}{ll}
 *    \frac{(x-a)}{\beta^2}  e^{-(x-a)^2/2\beta^2} \quad
 *    & 
 *    \mbox{ for } x \ge a 
 *    \\ 
 *    0 
 *    & 
 *    \mbox{ for } x < a, 
 *   \end{array} \right.
 * @f]
 * where @f$\beta> 0@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class RayleighGen extends RandomVariateGen {
   private double a;
   private double beta;

   /**
    * Creates a Rayleigh random variate generator with parameters @f$a
    * =@f$ `a` and @f$\beta=@f$ `beta`, using stream `s`.
    */
   public RayleighGen (RandomStream s, double a, double beta) {
      super (s, new RayleighDist(a, beta));
      setParams (a, beta);
   }

   /**
    * Creates a Rayleigh random variate generator with parameters @f$a =
    * 0@f$ and @f$\beta= @f$ `beta`, using stream `s`.
    */
   public RayleighGen (RandomStream s, double beta) {
      this (s, 0.0, beta);
   }

   /**
    * Creates a new generator for the Rayleigh distribution `dist` and
    * stream `s`.
    */
   public RayleighGen (RandomStream s, RayleighDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getA(), dist.getSigma());
   }

   /**
    * Uses inversion to generate a new variate from the Rayleigh
    * distribution with parameters @f$a = @f$&nbsp;`a` and @f$\beta=
    * @f$&nbsp;`beta`, using stream `s`.
    */
   public static double nextDouble (RandomStream s, double a, double beta) {
       return RayleighDist.inverseF (a, beta, s.nextDouble());
   }

   /**
    * Returns the parameter @f$a@f$.
    */
   public double getA() {
      return a;
   }

   /**
    * Returns the parameter @f$\beta@f$.
    */
   public double getSigma() {
      return beta;
   }

   /**
    * Sets the parameters @f$a = @f$&nbsp;`a` and
    * @f$\beta=@f$&nbsp;`beta` for this object.
    */
   public void setParams (double a, double beta) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      this.a  = a;
      this.beta = beta;
   }

}