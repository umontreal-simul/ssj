/*
 * Class:        MixtureGen
 * Description:  random variate generators for a mixture of distributions
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

import umontreal.ssj.probdist.DiscreteDistribution;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.probdist.Distribution;

/**
 * This class implements random variate generators for a mixture distribution.
 * Consider a random variable @f$X@f$ defined by a mixture of @f$k@f$ distributions.
 * Let @f$f_i(x)@f$ be the density or probability mass function of the @f$i@f$-th distribution with associated
 * weight @f$w_i \in [0, 1]@f$.
 * The sum of the @f$k@f$ weights must equal to 1, i.e., @f$\sum_{i=1}^k w_i = 1@f$.
 * 
 * The density or probability mass function @f$f(x)@f$ of @f$X@f$ is then
 * @f[
 *   f(x) = \sum_{i=1}^k f_i(x).
 * @f]
 * 
 * In the class constructor, the user must give an array of size @f$k@f$ of @ref umontreal.ssj.probdist.Distribution and an array
 * of `double` of same length for the associated weights @f$w_i@f$.
 * To generate a random variate, the implementation of this class will first select
 * randomly (and proportionally to the weights @f$w_i@f$) a distribution @f$D_i@f$ from the mixture.
 * Then, it will generate a random variate by inversion by calling the method `Distribution.inverseF`.
 *
 * @ingroup randvar_continuous
 */
public class MixtureGen extends RandomVariateGen {

   /**
    * The different distributions that compose this mixture.
    */
   protected Distribution[] dists;
   
   /**
    * The discrete distribution that is used to select randomly a distribution from the mixture
    * when generating a random variate.
    */
   protected DiscreteDistribution weightsDist;  // it is used to select the distribution
   
   /**
    * The probability of each distribution that compose this mixture, see variable #dists .
    */
   protected double[] weights;

   /**
    * Creates a new mixture distribution generator defined by the given
    * @ref Distribution and weights, using the random stream `s`.
    * The list of distributions and associated weights are given by arguments `dists` and `weights`.
    * The weights must be non-negative, and they must sum to 1.
    * 
    * @param s the random stream that will be used to generate the random variates
    * @param dists the distributions that define the mixture
    * @param weights the probability to select each distribution in the mixture 
    */
   public MixtureGen (RandomStream s, Distribution[] dists, double[] weights) {
      super (s, null);
      if (dists.length != weights.length)
         throw new IllegalArgumentException ("The arrays dists and weigths must have the same length");
      
      this.dists = dists;
      this.weights = weights;
      
      initWeightDistribution();
   }

   /**
    * Creates the discrete distribution that will be used to select randomly the distribution from the mixture.
    */
   private void initWeightDistribution() {
      int[] idx = new int[weights.length];
      for (int i = 0; i < idx.length; i++)
         idx[i] = i;
      weightsDist = new DiscreteDistribution(idx, weights, idx.length);
   }

   /**
    * @name Mixture distribution parameters
    * @{
    */

   /**
    * This method is not supported, it will throw an `UnsupportedOperationException`. Use method #getDistributions instead.
    * 
    * @throws UnsupportedOperationException this method is not supported
    */
   @Override
   public Distribution getDistribution() {
      throw new UnsupportedOperationException("Use getDistributions method instead");
   }
   
   /**
    * Returns the distributions of this mixture.
    * See also #getWeights .
    */
   public Distribution[] getDistributions() {
      return dists;
   }
   
   /**
    * Returns the probability associated with each distribution of the mixture.
    * See also #getDistributions .
    */
   public double[] getWeights() {
      return weights;
   }


   /**
    * @}
    */

   @Override
   public double nextDouble() {
      // select distribution from mixture
      int idx = (int) Math.round(weightsDist.inverseF(stream.nextDouble()));
      
      // generate random variate from selected distribution
      return dists[idx].inverseF (stream.nextDouble());
   }
}

