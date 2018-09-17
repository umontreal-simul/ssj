/*
 * Class:        NormalInverseFromDensityGen
 * Description:  random variate generators using numerical inversion of
                 the normal density
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
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
 * This class implements *normal* random variate generators using numerical
 * inversion of the *normal* density as described in @cite rDER10a&thinsp;.
 * It makes use of the class
 * @ref umontreal.ssj.probdist.InverseDistFromDensity. A set of tables are
 * precomputed to speed up the generation of normal random variables by
 * numerical inversion. This will be useful if one wants to generate a large
 * number of random variables.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class NormalInverseFromDensityGen extends NormalGen {

   /**
    * Creates a normal random variate generator with parameters
    * @f$\mu=@f$ `mu` and @f$\sigma@f$ = `sigma`, using stream `stream`.
    * It uses numerical inversion with precomputed tables. The
    * @f$u@f$-resolution `ueps` is the desired absolute error in the
    * `cdf`, and `order` is the degree of the Newton interpolating
    * polynomial over each interval.
    */
   public NormalInverseFromDensityGen (RandomStream stream, double mu, 
                                       double sigma, double ueps, int order) {
      // dist is the normal distribution
      super (stream, mu, sigma);
      double xc = mu;

      // member (NormalDist) dist is replaced by 
      // (InverseDistFromDensity) dist
      dist = new InverseDistFromDensity ((ContinuousDistribution) dist,
                                         xc, ueps, order);
    }

   /**
    * Similar to the first constructor, with the normal distribution
    * `dist`.
    */
   public NormalInverseFromDensityGen (RandomStream stream, NormalDist dist,
                                       double ueps, int order) {
      super (stream, dist);
      double xc = mu;

      // member (NormalDist) dist is replaced by 
      // (InverseDistFromDensity) dist
      this.dist = new InverseDistFromDensity (dist, xc, ueps, order);
   }

   /**
    * Creates a new normal generator using the *normal* distribution
    * `dist` and stream `stream`. `dist` may be obtained by calling method
    * #getDistribution, after using one of the other constructors to
    * create the precomputed tables. This is useful when one needs many
    * generators using the same normal distribution. Precomputing tables
    * for numerical inversion is costly; thus using only one set of tables
    * for many generators is more efficient. The first
    * @ref NormalInverseFromDensityGen generator using the other
    * constructors creates the precomputed tables. Then all other streams
    * use this constructor with the same set of tables.
    */
   public NormalInverseFromDensityGen (RandomStream stream, 
                                       InverseDistFromDensity dist) {
      super (stream, null);
      mu = dist.getXc();
      this.dist = dist;
   }

   /**
    * Returns the @f$u@f$-resolution `ueps`.
    */
   public double getUepsilon() {
      return ((InverseDistFromDensity)dist).getEpsilon();
   }

   /**
    * Returns the order of the interpolating polynomial.
    */
   public int getOrder() {
      return ((InverseDistFromDensity)dist).getOrder();
   }
}