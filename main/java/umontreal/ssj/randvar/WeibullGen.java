/*
 * Class:        WeibullGen
 * Description:  Weibull random number generator
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
 * This class implements random variate generators for the *Weibull*
 * distribution. Its density is
 * @anchor REF_randvar_WeibullGen_eq_weibull
 * @f[
 *   f(x) = \alpha\lambda^{\alpha} (x - \delta)^{\alpha- 1} \exp[-(\lambda(x-\delta))^{\alpha}] \qquad\mbox{ for } x>\delta, \tag{weibull}
 * @f]
 * and @f$f(x)=0@f$ elsewhere, where @f$\alpha> 0@f$, and @f$\lambda> 0@f$.
 *
 * The (non-static) `nextDouble` method simply calls `inverseF` on the
 * distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class WeibullGen extends RandomVariateGen {
   private double alpha = -1.0;
   private double lambda = -1.0;
   private double delta = -1.0;

   /**
    * Creates a Weibull random variate generator with parameters
    * @f$\alpha=@f$ `alpha`, @f$\lambda@f$ = `lambda` and @f$\delta@f$
    * = `delta`, using stream `s`.
    */
   public WeibullGen (RandomStream s, double alpha, double lambda,
                                      double delta) {
      super (s, new WeibullDist(alpha, lambda, delta));
      setParams (alpha, lambda, delta);
   }

   /**
    * Creates a Weibull random variate generator with parameters
    * @f$\alpha=@f$ `alpha`, @f$\lambda= 1@f$ and @f$\delta= 0@f$,
    * using stream `s`.
    */
   public WeibullGen (RandomStream s, double alpha) {
      this (s, alpha, 1.0, 0.0);
   }

   /**
    * Creates a new generator for the Weibull distribution `dist` and
    * stream `s`.
    */
   public WeibullGen (RandomStream s, WeibullDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getAlpha(), dist.getLambda(), dist.getDelta());
   }

   /**
    * Uses inversion to generate a new variate from the Weibull
    * distribution with parameters @f$\alpha= @f$&nbsp;`alpha`,
    * @f$\lambda= @f$&nbsp;`lambda`, and @f$\delta= @f$&nbsp;`delta`,
    * using stream `s`.
    */
   public static double nextDouble (RandomStream s, double alpha,
                                    double lambda, double delta) {
       return WeibullDist.inverseF (alpha, lambda, delta, s.nextDouble());
   }

   /**
    * Returns the parameter @f$\alpha@f$.
    */
   public double getAlpha() {
      return alpha;
   }

   /**
    * Returns the parameter @f$\lambda@f$.
    */
   public double getLambda() {
      return lambda;
   }

   /**
    * Returns the parameter @f$\delta@f$.
    */
   public double getDelta() {
      return delta;
   }

   /**
    * Sets the parameters @f$\alpha@f$, @f$\lambda@f$ and @f$\delta@f$
    * for this object.
    */
   public void setParams (double alpha, double lambda, double delta) {
      if (alpha <= 0.0)
        throw new IllegalArgumentException ("alpha <= 0");
      if (lambda <= 0.0)
        throw new IllegalArgumentException ("lambda <= 0");
      this.alpha  = alpha;
      this.lambda = lambda;
      this.delta  = delta;
   }
}