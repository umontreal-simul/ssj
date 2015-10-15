/*
 * Class:        LognormalGen
 * Description:  random variate generator for the lognormal distribution
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
 * This class implements methods for generating random variates from the
 * *lognormal* distribution. Its density is
 * @anchor REF_randvar_LognormalGen_eq_flognormal
 * @f[
 *   f(x) = \frac{1}{\sqrt{2\pi}\sigma x} e^{-(\ln(x) - \mu)^2/(2\sigma^2)} \qquad\mbox{ for }x>0, \tag{flognormal}
 * @f]
 * where @f$\sigma>0@f$.
 *
 * The (non-static) `nextDouble` method simply calls `inverseF` on the
 * lognormal distribution object. One can also generate a lognormal random
 * variate @f$X@f$ via
 * @f[
 *   \mbox{\texttt{X = Math.exp (NormalGen.nextDouble (s, mu, sigma))}},
 * @f]
 * in which `NormalGen` can actually be replaced by any subclass of
 * `NormalGen`.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class LognormalGen extends RandomVariateGen {
   private double mu;
   private double sigma = -1.0;

   /**
    * Creates a lognormal random variate generator with parameters
    * @f$\mu=@f$ `mu` and @f$\sigma=@f$ `sigma`, using stream `s`.
    */
   public LognormalGen (RandomStream s, double mu, double sigma) {
      this (s, new LognormalDist(mu, sigma));
      setParams (mu, sigma);
   }

   /**
    * Creates a lognormal random variate generator with parameters
    * @f$\mu= 0@f$ and @f$\sigma= 1@f$, using stream `s`.
    */
   public LognormalGen (RandomStream s) {
      this (s, 0.0, 1.0);
   }

   /**
    * Create a random variate generator for the lognormal distribution
    * `dist` and stream `s`.
    */
   public LognormalGen (RandomStream s, LognormalDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getMu(), dist.getSigma());
   }

   /**
    * Generates a new variate from the *lognormal* distribution with
    * parameters @f$\mu= @f$&nbsp;`mu` and @f$\sigma= @f$&nbsp;`sigma`,
    * using stream `s`.
    */
   public static double nextDouble (RandomStream s, double mu, double sigma) {
      return LognormalDist.inverseF (mu, sigma, s.nextDouble());
   }

   /**
    * Returns the parameter @f$\mu@f$ of this object.
    */
   public double getMu() {
      return mu;
   }

   /**
    * Returns the parameter @f$\sigma@f$ of this object.
    */
   public double getSigma() {
      return sigma;
   }

   /**
    * Sets the parameters @f$\mu@f$ and @f$\sigma@f$ of this object.
    */
   protected void setParams (double mu, double sigma) {
      if (sigma <= 0)
         throw new IllegalArgumentException ("sigma <= 0");
      this.mu = mu;
      this.sigma = sigma;
   }
}