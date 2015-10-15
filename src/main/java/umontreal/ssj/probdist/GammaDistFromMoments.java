/*
 * Class:        GammaDistFromMoments
 * Description:  gamma distribution
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
 * Extends the  @ref GammaDist distribution with constructors accepting the
 * mean @f$\mu@f$ and variance @f$\sigma^2@f$ as arguments instead of a
 * shape parameter @f$\alpha@f$ and a scale parameter @f$\lambda@f$. Since
 * @f$\mu=\alpha/ \lambda@f$, and @f$\sigma^2=\alpha/ \lambda^2@f$, the
 * shape and scale parameters are @f$\alpha=\mu^2 / \sigma^2@f$, and
 * @f$\lambda=\mu/ \sigma^2@f$, respectively.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class GammaDistFromMoments extends GammaDist {

   /**
    * Constructs a gamma distribution with mean `mean`, variance `var`,
    * and `d` decimal of precision.
    *  @param mean         the desired mean.
    *  @param var          the desired variance.
    *  @param d            the number of decimals of precision.
    */
   public GammaDistFromMoments (double mean, double var, int d) {
      super (mean * mean / var, mean / var, d);
   }

   /**
    * Constructs a gamma distribution with mean `mean`, and variance
    * `var`.
    *  @param mean         the desired mean.
    *  @param var          the desired variance.
    */
   public GammaDistFromMoments (double mean, double var) {
      super (mean * mean / var, mean / var);
   }
}