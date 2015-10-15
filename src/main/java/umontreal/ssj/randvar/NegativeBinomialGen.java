/*
 * Class:        NegativeBinomialGen
 * Description:  random variate generators for the negative binomial distribution
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
 * This class implements random variate generators having the *negative
 * binomial* distribution. Its mass function is
 * @anchor REF_randvar_NegativeBinomialGen_eq_fmass_negbin
 * @f[
 *   p(x) = \frac{\Gamma(\gamma+ x)}{x!\: \Gamma(\gamma)}\: p^{\gamma}(1 - p)^x, \qquad\mbox{for } x = 0, 1, 2, …\tag{fmass-negbin}
 * @f]
 * where @f$\Gamma@f$ is the gamma function, @f$\gamma> 0@f$ and
 * @f$0\le p\le1@f$. No local copy of the parameters @f$\gamma@f$ and
 * @f$p@f$ is maintained in this class. The (non-static) `nextInt` method
 * simply calls `inverseF` on the distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_discrete
 */
public class NegativeBinomialGen extends RandomVariateGenInt {
   protected double gamma;
   protected double p;

   /**
    * Creates a negative binomial random variate generator with parameters
    * @f$\gamma= @f$ `gamma` and @f$p@f$, using stream `s`.
    */
   public NegativeBinomialGen (RandomStream s, double gamma, double p) {
      super (s, new NegativeBinomialDist (gamma, p));
      setParams (gamma, p);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public NegativeBinomialGen (RandomStream s, NegativeBinomialDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getGamma(), dist.getP());
   }

   /**
    * Generates a new variate from the *negative binomial* distribution,
    * with parameters @f$\gamma= @f$&nbsp;`gamma` and @f$p =
    * @f$&nbsp;`p`, using stream `s`.
    */
   public static int nextInt (RandomStream s, double gamma, double p) {
      return NegativeBinomialDist.inverseF (gamma, p, s.nextDouble());
   }

   /**
    * Returns the parameter @f$\gamma@f$ of this object.
    */
   public double getGamma() {
      return gamma;
   }

   /**
    * Returns the parameter @f$p@f$ of this object.
    */
   public double getP() {
      return p;
   }

   /**
    * Sets the parameter @f$\gamma@f$ and @f$p@f$ of this object.
    */
   protected void setParams (double gamma, double p) {
      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in [0, 1]");
      if (gamma <= 0.0)
         throw new IllegalArgumentException ("gamma <= 0");
      this.p = p;
      this.gamma = gamma;
   }
}