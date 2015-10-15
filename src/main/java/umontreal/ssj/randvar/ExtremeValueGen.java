/*
 * Class:        ExtremeValueGen
 * Description:  random variate generators for the Gumbel distribution
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
 * <strong>This class has been replaced by  @ref GumbelGen </strong>.
 *
 * This class implements random variate generators for the *Gumbel* (or
 * *extreme value*) distribution. Its density is
 * @anchor REF_randvar_ExtremeValueGen_eq_fextremevalue
 * @f[
 *   f(x) = \lambda e^{-e^{-\lambda(x - \alpha)}-\lambda(x-\alpha)} \tag{fextremevalue}
 * @f]
 * where @f$\lambda>0@f$.
 *
 * The (non-static) `nextDouble` method simply calls `inverseF` on the
 * distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
@Deprecated
public class ExtremeValueGen extends RandomVariateGen {
   protected double alpha = -1.0;
   protected double lambda = -1.0;

   /**
    * Creates an *extreme value* random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\lambda@f$ = `lambda`, using stream
    * `s`.
    */
   public ExtremeValueGen (RandomStream s, double alpha, double lambda) {
      super (s, new ExtremeValueDist(alpha, lambda));
      setParams (alpha, lambda);
   }

   /**
    * Creates an *extreme value* random variate generator with parameters
    * @f$\alpha= 0@f$ and @f$\lambda=1@f$, using stream `s`.
    */
   public ExtremeValueGen (RandomStream s) {
      this (s, 0.0, 1.0);
   }

   /**
    * Creates a new generator object for distribution `dist` and stream
    * `s`.
    */
   public ExtremeValueGen (RandomStream s, ExtremeValueDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getAlpha(), dist.getLambda());
   }

   /**
    * Uses inversion to generate a new variate from the extreme value
    * distribution with parameters @f$\alpha= @f$&nbsp;`alpha` and
    * @f$\lambda= @f$&nbsp;`lambda`, using stream `s`.
    */
   public static double nextDouble (RandomStream s, double alpha,
                                    double lambda) {
      return ExtremeValueDist.inverseF (alpha, lambda, s.nextDouble());
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