/*
 * Class:        LaplaceGen
 * Description:  generator of random variates from the Laplace distribution
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
 * *Laplace* distribution. Its density is
 *  (see @cite tJOH95b&thinsp; (page 165))
 * @anchor REF_randvar_LaplaceGen_eq_flaplace
 * @f[
 *   f(x) = \frac{1}{2\beta}e^{-|x-\mu|/\beta} \qquad\mbox{for } -\infty< x < \infty, \tag{flaplace}
 * @f]
 * where @f$\beta> 0@f$.
 *
 * The (non-static) `nextDouble` method simply calls `inverseF` on the
 * distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class LaplaceGen extends RandomVariateGen {
   private double mu;
   private double beta;

   /**
    * Creates a Laplace random variate generator with parameters
    * @f$\mu@f$ = `mu` and @f$\beta@f$ = `beta`, using stream `s`.
    */
   public LaplaceGen (RandomStream s, double mu, double beta) {
      super (s, new LaplaceDist(mu, beta));
      setParams (mu, beta);
   }

   /**
    * Creates a Laplace random variate generator with parameters @f$\mu=
    * 0@f$ and @f$\beta= 1@f$, using stream `s`.
    */
   public LaplaceGen (RandomStream s) {
      this (s, 0.0, 1.0);
   }

   /**
    * Creates a new generator for the Laplace distribution `dist` and
    * stream `s`.
    */
   public LaplaceGen (RandomStream s, LaplaceDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getMu(), dist.getBeta());
   }

   /**
    * Generates a new variate from the Laplace distribution with
    * parameters @f$\mu= @f$&nbsp;`mu` and @f$\beta= @f$&nbsp;`beta`,
    * using stream `s`.
    */
   public static double nextDouble (RandomStream s, double mu, double beta) {
      return LaplaceDist.inverseF (mu, beta, s.nextDouble());
   }

   /**
    * Returns the parameter @f$\mu@f$.
    */
   public double getMu() {
      return mu;
   }

   /**
    * Returns the parameter @f$\beta@f$.
    */
   public double getBeta() {
      return beta;
   }

   /**
    * Sets the parameters @f$\mu@f$ and @f$\beta@f$ of this object.
    */
   protected void setParams (double mu, double beta) {
     if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      this.mu = mu;
      this.beta = beta;
   }
}