/*
 * Class:        HyperbolicSecantGen
 * Description:  random variate generators for the hyperbolic secant distribution
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
 * This class implements random variate generators for the *hyperbolic
 * secant* distribution with location parameter @f$\mu@f$ and scale
 * parameter @f$\sigma@f$. The density function of this distribution is
 * @anchor REF_randvar_HyperbolicSecantGen_eq_fHyperbolicSecant
 * @f[
 *   f(x) = \frac{1}{2 \sigma} \mbox{ sech}\left(\frac{\pi}{2} \frac{(x - \mu)}{\sigma}\right), \qquad-\infty<x < \infty. \tag{fHyperbolicSecant}
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class HyperbolicSecantGen extends RandomVariateGen {
   protected double mu;
   protected double sigma;

   /**
    * Creates a *hyperbolic secant* random variate generator with
    * parameters @f$\mu=@f$ `mu` and @f$\sigma=@f$ `sigma`, using stream
    * `s`.
    */
   public HyperbolicSecantGen (RandomStream s, double mu, double sigma) {
      super (s, new HyperbolicSecantDist(mu, sigma));
      setParams (mu, sigma);
   }

   /**
    * Creates a *hyperbolic secant* random variate generator with
    * parameters @f$\mu=0@f$ and @f$\sigma=1@f$, using stream `s`.
    */
   public HyperbolicSecantGen (RandomStream s) {
      this (s, 0.0, 1.0);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public HyperbolicSecantGen (RandomStream s, HyperbolicSecantDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getMu(), dist.getSigma());
   }

   /**
    * Generates a variate from the *hyperbolic secant* distribution with
    * location parameter @f$\mu@f$ and scale parameter @f$\sigma@f$.
    */
   public static double nextDouble (RandomStream s, double mu, double sigma) {
      return HyperbolicSecantDist.inverseF (mu, sigma, s.nextDouble());
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
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      this.mu = mu;
      this.sigma = sigma;
   }
}