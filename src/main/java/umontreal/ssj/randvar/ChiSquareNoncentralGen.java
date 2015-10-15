/*
 * Class:        ChiSquareNoncentralGen
 * Description:  random variate generators for the noncentral chi square distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
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
 * This class implements random variate generators for the *noncentral chi
 * square* distribution with @f$\nu> 0@f$ degrees of freedom and
 * noncentrality parameter @f$\lambda> 0@f$. See the definition in
 * @ref umontreal.ssj.probdist.ChiSquareNoncentralDist.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class ChiSquareNoncentralGen extends RandomVariateGen {
   protected double nu = -1.0;
   protected double lambda = -1.0;

   /**
    * Creates a *noncentral chi square* random variate generator with `nu`
    * @f$=\nu>0@f$ degrees of freedom and noncentrality parameter
    * `lambda` @f$= \lambda>0@f$, using stream `s`.
    */
   public ChiSquareNoncentralGen (RandomStream s, double nu, double lambda) {
      super (s, new ChiSquareNoncentralDist(nu, lambda));
      setParams (nu, lambda);
   }

   /**
    * Create a new generator for the distribution `dist` and stream `s`.
    */
   public ChiSquareNoncentralGen (RandomStream s,
                                  ChiSquareNoncentralDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getNu (), dist.getLambda());
   }

   /**
    * Generates a new variate from the noncentral chi square distribution
    * with `nu` = @f$\nu@f$ degrees of freedom and noncentrality
    * parameter `lambda` @f$=\lambda@f$, using stream `s`.
    */
   public static double nextDouble (RandomStream s,
                                    double nu, double lambda) {
      return ChiSquareNoncentralDist.inverseF (nu, lambda, s.nextDouble());
   }

   /**
    * Returns the value of @f$\nu@f$ of this object.
    */
   public double getNu() {
      return nu;
   }

   /**
    * Returns the value of @f$\lambda@f$ for this object.
    */
   public double getLambda() {
      return lambda;
   }

   /**
    * Sets the parameters @f$\nu=@f$ `nu` and @f$\lambda= @f$ `lambda`
    * of this object.
    */
   protected void setParams (double nu, double lambda) {
      if (nu <= 0.0)
         throw new IllegalArgumentException ("nu <= 0");
      if (lambda < 0.0)
         throw new IllegalArgumentException ("lambda < 0");
      this.nu = nu;
      this.lambda = lambda;
   }
}