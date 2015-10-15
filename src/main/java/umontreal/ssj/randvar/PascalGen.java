/*
 * Class:        PascalGen
 * Description:  Pascal random variate generators
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
 * Implements Pascal random variate generators, which is a special case of
 * the negative binomial generator with parameter @f$\gamma@f$ equal to a
 * positive integer. See  @ref NegativeBinomialGen for a description.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_discrete
 */
public class PascalGen extends RandomVariateGenInt {
   protected int    n;
   protected double p;

   /**
    * Creates a Pascal random variate generator with parameters @f$n@f$
    * and @f$p@f$, using stream `s`.
    */
   public PascalGen (RandomStream s, int n, double p) {
      super (s, new PascalDist (n, p));
      setParams (n, p);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public PascalGen (RandomStream s, PascalDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getN1(), dist.getP());
   }

   /**
    * Generates a new variate from the *Pascal* distribution, with
    * parameters @f$n = @f$&nbsp;`n` and @f$p = @f$&nbsp;`p`, using stream
    * `s`.
    */
   public static int nextInt (RandomStream s, int n, double p) {
      return PascalDist.inverseF (n, p, s.nextDouble());
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
   public double getP() {
      return p;
   }

   /**
    * Sets the parameter @f$n@f$ and @f$p@f$ of this object.
    */
   protected void setParams (int n, double p) {
      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in [0, 1]");
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      this.p = p;
      this.n = n;
   }
}