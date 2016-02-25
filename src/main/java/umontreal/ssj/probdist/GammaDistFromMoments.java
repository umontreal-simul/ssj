/*
 * Class:        GammaDistFromMoments
 * Description:  gamma distribution
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