/*
 * Class:        LoglogisticGen
 * Description:  random variate generators for the log-logistic distribution
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
 * This class implements random variate generators for the *log-logistic*
 * distribution with shape parameter @f$\alpha> 0@f$ and scale parameter
 * @f$\beta> 0@f$. The density function of this distribution is
 * @anchor REF_randvar_LoglogisticGen_eq_floglogistic
 * @f[
 *   f(x) = \frac{\alpha(x / \beta)^{\alpha- 1}}{\beta[1 + (x / \beta)^{\alpha}]^2} \qquad\qquad\mbox{for } x > 0. \tag{floglogistic}
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class LoglogisticGen extends RandomVariateGen {
   protected double alpha;
   protected double beta;

   /**
    * Creates a log-logistic random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\beta=@f$ `beta`, using stream `s`.
    */
   public LoglogisticGen (RandomStream s, double alpha, double beta) {
      super (s, new LoglogisticDist(alpha, beta));
      setParams (alpha, beta);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public LoglogisticGen (RandomStream s, LoglogisticDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getAlpha(), dist.getBeta());
   }

   /**
    * Generates a variate from the *log-logistic* distribution with shape
    * parameter @f$\alpha> 0@f$ and scale parameter @f$\beta> 0@f$.
    */
   public static double nextDouble (RandomStream s,
                                    double alpha, double beta) {
      return LoglogisticDist.inverseF (alpha, beta, s.nextDouble());
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
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      this.alpha = alpha;
      this.beta = beta;
   }
}