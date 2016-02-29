/*
 * Class:        GumbelGen
 * Description:  generator of random variates from the Gumbel distribution 
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package umontreal.ssj.randvar;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;

/**
 * This class implements methods for generating random variates from the
 * *Gumbel* distribution. Its density is given by
 * @anchor REF_randvar_GumbelGen_eq_densgumbel
 * @f[
 *   f (x) = \frac{e^{-z} e^{-e^{-z}}}{|\beta|}, \qquad\mbox{for } -\infty< x < \infty, 
 *    \tag{densgumbel}
 * @f]
 * where we use the notation @f$z = (x-\delta)/\beta@f$. The scale
 * parameter @f$\beta@f$ can be positive (for the Gumbel distribution) or
 * negative (for the reverse Gumbel distribution), but not 0.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class GumbelGen extends RandomVariateGen {
   private double delta;
   private double beta;

   /**
    * Creates a Gumbel random number generator with @f$\beta= 1@f$ and
    * @f$\delta= 0@f$ using stream `s`.
    */
   public GumbelGen (RandomStream s) {
      this (s, 1.0, 0.0);
   }

   /**
    * Creates a Gumbel random number generator with parameters
    * @f$\beta@f$ = `beta` and @f$\delta@f$ = `delta` using stream `s`.
    */
   public GumbelGen (RandomStream s, double beta, double delta) {
      super (s, new GumbelDist(beta, delta));
      setParams (beta, delta);
   }

   /**
    * Creates a new generator for the Gumbel distribution `dist` and
    * stream `s`.
    */
   public GumbelGen (RandomStream s, GumbelDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getBeta(), dist.getDelta());
   }

   /**
    * Generates a new variate from the Gumbel distribution with parameters
    * @f$\beta= @f$&nbsp;`beta` and @f$\delta= @f$&nbsp;`delta` using
    * stream `s`.
    */
   public static double nextDouble (RandomStream s, double beta, double delta) {
      return GumbelDist.inverseF (beta, delta, s.nextDouble());
   }

   /**
    * Returns the parameter @f$\beta@f$.
    */
   public double getBeta() {
      return beta;
   }

   /**
    * Returns the parameter @f$\delta@f$.
    */
   public double getDelta() {
      return delta;
   }

   /**
    * Sets the parameters @f$\beta@f$ and @f$\delta@f$ of this object.
    */
   protected void setParams (double beta, double delta) {
     if (beta == 0.0)
         throw new IllegalArgumentException ("beta = 0");
      this.delta = delta;
      this.beta = beta;
   }
}