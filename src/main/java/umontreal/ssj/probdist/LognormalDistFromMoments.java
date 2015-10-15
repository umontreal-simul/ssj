/*
 * Class:        LognormalDistFromMoments
 * Description:  lognormal distribution
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
package umontreal.ssj.probdist;

/**
 * Extends the  @ref LognormalDist class with a constructor accepting the
 * mean @f$m@f$ and the variance @f$v@f$ of the distribution as arguments.
 * The mean and variance of a lognormal random variable with parameters
 * @f$\mu@f$ and @f$\sigma@f$ are @f$e^{\mu+\sigma^2/2}@f$ and
 * @f$e^{2\mu+ \sigma^2}(e^{\sigma^2} - 1)@f$ respectively, so the
 * parameters are given by
 *  @f$\sigma=\sqrt{\ln(v/m^2+1)}@f$
 *   and @f$\mu=\ln(m) - \sigma^2/2@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class LognormalDistFromMoments extends LognormalDist {
   public LognormalDistFromMoments (double mean, double var) {
      super (getMu (mean, var), Math.sqrt (getSigma2 (mean, var)));
   }

   private static double getMu (double mean, double var) {
      final double sigma2 = getSigma2 (mean, var);
      return Math.log (mean) - sigma2 / 2.0;
   }

   private static double getSigma2 (double mean, double var) {
      if (mean <= 0)
         throw new IllegalArgumentException ("Mean must be positive");
      if (var <= 0)
         throw new IllegalArgumentException ("Variance must be positive");
      return Math.log (var / (mean * mean) + 1);
   }
}