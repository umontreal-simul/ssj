/*
 * Class:        LaplaceGen
 * Description:  generator of random variates from the Laplace distribution
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