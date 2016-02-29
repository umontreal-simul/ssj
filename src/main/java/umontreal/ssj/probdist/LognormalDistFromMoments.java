/*
 * Class:        LognormalDistFromMoments
 * Description:  lognormal distribution
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