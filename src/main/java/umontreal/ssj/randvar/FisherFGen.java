/*
 * Class:        FisherFGen
 * Description:  random variate generators for the Fisher F distribution
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
 * This class implements random variate generators for the *Fisher F*
 * distribution with @f$n@f$ and @f$m@f$ degrees of freedom, where @f$n@f$
 * and @f$m@f$ are positive integers. The density function of this
 * distribution is
 * @anchor REF_randvar_FisherFGen_eq_FisherF
 * @f[
 *   f(x) = \frac{\Gamma({(n + m)}/2)n^{n/2}m^{m/2}}{\Gamma(n/2)\Gamma(m/2)} \frac{x^{{(n - 2)}/2}}{(m + nx)^{{(n + m)}/2}}, \qquad\mbox{for } x > 0 \tag{FisherF}
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class FisherFGen extends RandomVariateGen {
   protected int    n = -1;
   protected int    m = -1;

   /**
    * Creates a *Fisher F* random variate generator with @f$n@f$ and
    * @f$m@f$ degrees of freedom, using stream `s`.
    */
   public FisherFGen (RandomStream s, int n, int m) {
      super (s, new FisherFDist(n, m));
      setParams (n, m);
      }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public FisherFGen (RandomStream s, FisherFDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getN(), dist.getM());
   }

   /**
    * Generates a variate from the *Fisher F* distribution with @f$n@f$
    * and @f$m@f$ degrees of freedom, using stream `s`.
    */
   public static double nextDouble (RandomStream s, int n, int m) {
      return FisherFDist.inverseF (n, m, 15, s.nextDouble());
   }

   /**
    * Returns the parameter @f$n@f$ of this object.
    */
   public int getN() {
      return n;
   }

   /**
    * Returns the parameter @f$p@f$ of this object.
    */
   public int getM() {
      return m;
   }

   /**
    * Sets the parameters @f$n@f$ and @f$m@f$ of this object.
    */
   protected void setParams (int n, int m) {
      if (m <= 0)
         throw new IllegalArgumentException ("m <= 0");
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      this.m = m;
      this.n = n;
   }
}