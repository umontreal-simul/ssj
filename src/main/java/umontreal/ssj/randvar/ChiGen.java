/*
 * Class:        ChiGen
 * Description:  random variate generators for the chi distribution
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
 * This class implements random variate generators for the *chi*
 * distribution. It has @f$\nu>0@f$ degrees of freedom and its density
 * function is
 * (see @cite tJOH95a&thinsp;, page 417)
 * @anchor REF_randvar_ChiGen_eq_Fchi
 * @f[
 *   f (x) = \frac{e^{-x^2 /2} x^{\nu-1}}{2^{(\nu/2) - 1}\Gamma(\nu/2)} \qquad\mbox{for } x > 0, \tag{Fchi}
 * @f]
 * where @f$\Gamma(x)@f$ is the gamma function defined in (
 * {@link REF_randvar_GammaGen_eq_Gamma Gamma} ).
 *
 * The (non-static) `nextDouble` method simply calls `inverseF` on the
 * distribution (slow).
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class ChiGen extends RandomVariateGen {
   protected int nu = -1;

   /**
    * Creates a *chi* random variate generator with @f$\nu=@f$ `nu`
    * degrees of freedom, using stream `s`.
    */
   public ChiGen (RandomStream s, int nu) {
      super (s, new ChiDist(nu));
      setParams (nu);
   }

   /**
    * Create a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public ChiGen (RandomStream s, ChiDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getNu ());
   }

   /**
    * Generates a random variate from the chi distribution with @f$\nu=
    * @f$&nbsp;`nu` degrees of freedom, using stream `s`.
    */
   public static double nextDouble (RandomStream s, int nu) {
      if (nu <= 0)
         throw new IllegalArgumentException ("nu <= 0");
      return ChiDist.inverseF (nu, s.nextDouble());
   }

   /**
    * Returns the value of @f$\nu@f$ for this object.
    */
   public int getNu() {
      return nu;
   }
   protected void setParams (int nu) {
      if (nu <= 0)
         throw new IllegalArgumentException ("nu <= 0");
      this.nu = nu;
   }
}