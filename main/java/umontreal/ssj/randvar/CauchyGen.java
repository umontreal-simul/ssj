/*
 * Class:        CauchyGen
 * Description:  random variate generators for the Cauchy distribution
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
 * This class implements random variate generators for the *Cauchy*
 * distribution. The density is
 * (see, e.g., @cite tJOH95a&thinsp; p. 299):
 * @anchor REF_randvar_CauchyGen_eq_fcauchy
 * @f[
 *   f (x) = \frac{\beta}{\pi[(x-\alpha)^2 + \beta^2]}, \qquad\mbox{for } -\infty< x < \infty, \tag{fcauchy}
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
public class CauchyGen extends RandomVariateGen {
   protected double alpha;
   protected double beta;

   /**
    * Creates a Cauchy random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\beta@f$ = `beta`, using stream `s`.
    */
   public CauchyGen (RandomStream s, double alpha, double beta) {
      super (s, new CauchyDist(alpha, beta));
      setParams(alpha, beta);
   }

   /**
    * Creates a Cauchy random variate generator with parameters
    * @f$\alpha=0 @f$ and @f$\beta= 1@f$, using stream `s`.
    */
   public CauchyGen (RandomStream s) {
      this (s, 0.0, 1.0);
   }

   /**
    * Create a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public CauchyGen (RandomStream s, CauchyDist dist) {
      super (s, dist);
      if (dist != null)
         setParams(dist.getAlpha(), dist.getBeta());
   }

   /**
    * Generates a new variate from the *Cauchy* distribution with
    * parameters @f$\alpha= @f$&nbsp;`alpha` and @f$\beta=
    * @f$&nbsp;`beta`, using stream `s`.
    */
   public static double nextDouble (RandomStream s,
                                    double alpha, double beta) {
      return CauchyDist.inverseF (alpha, beta, s.nextDouble());
   }

   /**
    * Returns the parameter @f$\alpha@f$ of this object.
    */
   public double getAlpha() {
      return alpha;
   }

   /**
    * Returns the parameter @f$\beta@f$ of this object.
    */
   public double getBeta() {
      return beta;
   }


   protected void setParams (double alpha, double beta) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      this.alpha = alpha;
      this.beta = beta;
   }
}