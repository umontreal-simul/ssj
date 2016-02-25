/*
 * Class:        LogisticGen
 * Description:  random variate generators for the logistic distribution
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
 * This class implements random variate generators for the *logistic*
 * distribution. Its parameters are @f$\alpha@f$ and @f$\lambda> 0@f$. Its
 * density function is
 * @anchor REF_randvar_LogisticGen_eq_flogistic
 * @f[
 *   f(x) = \frac{\lambda e^{-\lambda(x-\alpha)}}{\left(1 + e^{-\lambda(x-\alpha)}\right)^2} \qquad\mbox{ for } -\infty<x<\infty. \tag{flogistic}
 * @f]
 * The (non-static) `nextDouble` method simply calls `inverseF` on the
 * distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class LogisticGen extends RandomVariateGen {
   protected double alpha = -1.0;
   protected double lambda = -1.0;

   /**
    * Creates a logistic random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\lambda=@f$ `lambda`, using stream
    * `s`.
    */
   public LogisticGen (RandomStream s, double alpha, double lambda) {
      super (s, new LogisticDist(alpha, lambda));
      setParams (alpha, lambda);
   }

   /**
    * Creates a logistic random variate generator with parameters
    * @f$\alpha= 0@f$ and @f$\lambda=1@f$, using stream `s`.
    */
   public LogisticGen (RandomStream s) {
      this (s, 0.0, 1.0);
   }

   /**
    * Creates a new generator for the logistic distribution `dist` and
    * stream `s`.
    */
   public LogisticGen (RandomStream s, LogisticDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getAlpha(), dist.getLambda());
   }

   /**
    * Generates a new variate from the *logistic* distribution with
    * parameters @f$\alpha= @f$&nbsp;`alpha` and @f$\lambda=
    * @f$&nbsp;`lambda`, using stream `s`.
    */
   public static double nextDouble (RandomStream s,
                                    double alpha, double lambda) {
      return  LogisticDist.inverseF (alpha, lambda, s.nextDouble());
   }

   /**
    * Returns the parameter @f$\alpha@f$ of this object.
    */
   public double getAlpha() {
      return alpha;
   }

   /**
    * Returns the parameter @f$\lambda@f$ of this object.
    */
   public double getLambda() {
      return lambda;
   }

   /**
    * Sets the parameter @f$\alpha@f$ and @f$\lambda@f$ of this object.
    */
   protected void setParams (double alpha, double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      this.lambda = lambda;
      this.alpha = alpha;
   }
}