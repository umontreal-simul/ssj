/*
 * Class:        WeibullGen
 * Description:  Weibull random number generator
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