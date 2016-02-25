/*
 * Class:        RandomMultivariateGen
 * Description:  base class for multidimensional random variate generators
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
package umontreal.ssj.randvarmulti;

import umontreal.ssj.probdistmulti.ContinuousDistributionMulti;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.randvar.RandomVariateGen;

/**
 * This class is the multivariate counterpart of
 * @ref umontreal.ssj.randvar.RandomVariateGen. It is the base class for
 * general random variate generators over the @f$d@f$-dimensional real space
 * @f$\mathbb{R}^d@f$. It specifies the signature of the  #nextPoint method,
 * which is normally called to generate a random vector from a given
 * distribution. Contrary to univariate distributions and generators, here
 * the inversion method is not well defined, so we cannot construct a
 * multivariate generator simply by passing a multivariate distribution and a
 * stream; we must specify a generating method as well. For this reason, this
 * class is abstract. Generators can be constructed only by invoking the
 * constructor of a subclass. This is an important difference with
 * @ref umontreal.ssj.randvar.RandomVariateGen.
 *
 * <div class="SSJ-bigskip"></div>
 */
public abstract class RandomMultivariateGen {
   protected int dimension;
   // Careful here: there is also a RandomStream inside gen1. But only one
   // of these two is used in a given class.
   protected RandomStream stream;  // stream used to generate random numbers
   protected RandomVariateGen gen1; // 1-dim generator used to generate random variates
// This constructor is needed for subclasses with no associated distribution.
//   protected RandomMultivariateGen() {}

   /**
    * Generates a random point @f$p@f$ using the the stream contained in
    * this object.
    */
   abstract public void nextPoint (double[] p);

   /**
    * Generates @f$n@f$ random points. These points are stored in the
    * array `v`, starting at index `start`. Thus `v[start][i]` contains
    * coordinate @f$i@f$ of the first generated point. By default, this
    * method calls  #nextPoint @f$n@f$ times, but one can override it in
    * subclasses for better efficiency. The array argument `v[][d]` must
    * have @f$d@f$ elements reserved for each generated point before
    * calling this method.
    *  @param v            array in which the variates will be stored
    *  @param start        starting index, in `v`, of the new variates
    *  @param n            number of variates to generate
    */
   public void nextArrayOfPoints (double[][] v, int start, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n must be positive.");
      for (int i = 0; i < n; i++)
         nextPoint(v[start + i]);
   }

   /**
    * Returns the dimension of this multivariate generator (the dimension
    * of the random points).
    */
   public int getDimension() {
      return dimension;
   }

   /**
    * Returns the  @ref umontreal.ssj.rng.RandomStream used by this
    * object.
    *  @return the stream associated to this object
    */
   public RandomStream getStream() {
      if (null != gen1)
         return gen1.getStream();
      return stream;
   }

   /**
    * Sets the  @ref umontreal.ssj.rng.RandomStream used by this object to
    * `stream`.
    */
   public void setStream (RandomStream stream) {
      if (null != gen1)
         gen1.setStream(stream);
      else
         this.stream = stream;
   }

}