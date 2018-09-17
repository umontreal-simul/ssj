/*
 * Class:        ConstantDist
 * Description:  constant distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Éric Buist
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
 * Represents a *constant* discrete distribution taking a single real value
 * with probability 1. Its mass function is
 * @anchor REF_probdist_ConstantDist_eq_fcons
 * @f[
 *   p(x) = \left\{\begin{array}{ll}
 *    1, 
 *    & 
 *    \qquad\mbox{for } x = c,
 *    \\ 
 *    0, 
 *    & 
 *    \qquad\mbox{elsewhere. } 
 *   \end{array}\right. \tag{fcons}
 * @f]
 * Its distribution function is
 * @anchor REF_probdist_ConstantDist_eq_cdfcons
 * @f[
 *   F(x) = \left\{\begin{array}{ll}
 *    0, 
 *    & 
 *    \qquad\mbox{ for } x < c
 *    \\ 
 *    1, 
 *    & 
 *    \qquad\mbox{ for } x \ge c. 
 *   \end{array}\right. \tag{cdfcons}
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_discretereal
 */
public class ConstantDist extends DiscreteDistribution {
   private double c;

   /**
    * Constructs a new constant distribution with probability 1 at `c`.
    */
   public ConstantDist (double c) {
      super (new double[] { c }, new double[] { 1.0 }, 1);
      this.c = c;
   }

   /**
    * Returns the mean @f$E[X] = c@f$.
    *  @return @f$c@f$
    */
   @Override
   public double getMean() {
      return c;
   }

   /**
    * Returns the variance @f$\mbox{Var}[X] = 0@f$.
    *  @return 0
    */
   @Override
   public double getVariance()  {
      return 0;
   }

   /**
    * Returns the standard deviation = 0.
    *  @return 0
    */
   @Override
   public double getStandardDeviation() {
      return 0;
   }

   /**
    * Returns the inverse distribution function @f$c = F^{-1}(u)@f$.
    *  @return @f$c@f$
    */
   @Override
   public double inverseF (double u) {
      return c;
   }

}