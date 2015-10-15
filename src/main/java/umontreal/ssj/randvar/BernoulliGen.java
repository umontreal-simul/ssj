/*
 * Class:        BernoulliGen
 * Description:  random variate generators for the Bernoulli distribution
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
import umontreal.ssj.probdist.BernoulliDist;
import umontreal.ssj.rng.RandomStream;

/**
 * This class implements random variate generators for the *Bernoulli*
 * distribution (see class  @ref umontreal.ssj.probdist.BernoulliDist ).
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_discrete
 */
public class BernoulliGen extends RandomVariateGenInt {
   protected double p;

   /**
    * Creates a Bernoulli random variate generator with parameter @f$p@f$,
    * using stream `s`.
    */
   public BernoulliGen (RandomStream s, double p) {
      super (s, new BernoulliDist (p));
      setParams (p);
   }

   /**
    * Creates a random variate generator for the *Bernoulli* distribution
    * `dist` and the random stream `s`.
    */
   public BernoulliGen (RandomStream s, BernoulliDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getP());
   }

   /**
    * Generates a new integer from the *Bernoulli* distribution with
    * parameter @f$p = @f$&nbsp;`p`, using the given stream `s`.
    */
   public static int nextInt (RandomStream s, double p) {
      return BernoulliDist.inverseF (p, s.nextDouble());
   }

   /**
    * Returns the parameter @f$p@f$ of this object.
    */
   public double getP() {
      return p;
   }

   /**
    * Sets the parameter @f$p@f$ of this object.
    */
   protected void setParams (double p) {
      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in range [0, 1]");
      this.p = p;
   }
}