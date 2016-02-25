/*
 * Class:        FrechetGen
 * Description:  generator of random variates from the Fréchet distribution
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
 * <em>Fréchet</em> distribution, with location parameter @f$\delta@f$,
 * scale parameter @f$\beta> 0@f$, and shape parameter @f$\alpha> 0@f$,
 * where we use the notation @f$z = (x-\delta)/\beta@f$. It has density
 * @f[
 *   f (x) = \frac{\alpha e^{-z^{-\alpha}}}{\beta z^{\alpha+1}}, \qquad\mbox{for } x > \delta.
 * @f]
 * The density is 0 for @f$x \le\delta@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class FrechetGen extends RandomVariateGen {
   private double delta;
   private double beta;
   private double alpha;

   /**
    * Creates a <em>Fréchet</em> random number generator with
    * @f$\alpha=@f$ `alpha`, @f$\beta= 1@f$ and @f$\delta= 0@f$ using
    * stream `s`.
    */
   public FrechetGen (RandomStream s, double alpha) {
      this (s, alpha, 1.0, 0.0);
   }

   /**
    * Creates a <em>Fréchet</em> random number generator with parameters
    * @f$\alpha@f$ = `alpha`, @f$\beta@f$ = `beta` and @f$\delta@f$ =
    * `delta` using stream `s`.
    */
   public FrechetGen (RandomStream s, double alpha, double beta,
                      double delta) {
      super (s, new FrechetDist (alpha, beta, delta));
      setParams (alpha, beta, delta);
   }

   /**
    * Creates a new generator for the <em>Fréchet</em> distribution `dist`
    * and stream `s`.
    */
   public FrechetGen (RandomStream s, FrechetDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getAlpha(), dist.getBeta(), dist.getDelta());
   }

   /**
    * Generates a new variate from the <em>Fréchet</em> distribution with
    * parameters @f$\alpha=@f$ `alpha`, @f$\beta= @f$&nbsp;`beta` and
    * @f$\delta= @f$&nbsp;`delta` using stream `s`.
    */
   public static double nextDouble (RandomStream s, double alpha,
                                    double beta, double delta) {
      return FrechetDist.inverseF (alpha, beta, delta, s.nextDouble());
   }

   /**
    * Returns the parameter @f$\alpha@f$.
    */
   public double getAlpha() {
      return alpha;
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
    * Sets the parameters @f$\alpha@f$, @f$\beta@f$ and @f$\delta@f$ of
    * this object.
    */
   protected void setParams (double alpha, double beta, double delta) {
     if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      this.delta = delta;
      this.beta = beta;
      this.alpha = alpha;
   }
}