/*
 * Class:        Distribution
 * Description:  interface for all discrete and continuous distributions
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
 * This interface should be implemented by all classes supporting discrete
 * and continuous distributions. It specifies the signature of methods that
 * compute the distribution function @f$F(x)@f$, the complementary
 * distribution function @f$\bar{F}(x)@f$, and the inverse distribution
 * function @f$ F^{-1} (u)@f$. It also specifies the signature of methods
 * that returns the mean, the variance and the standard deviation.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_general
 */
public interface Distribution {

/**
 * Returns the distribution function @f$F(x)@f$.
 *  @param x            value at which the distribution function is evaluated
 *  @return distribution function evaluated at `x`
 */
public double cdf (double x);

   /**
    * Returns @f$\bar{F}(x) = 1 - F(x)@f$.
    *  @param x            value at which the complementary distribution
    *                      function is evaluated
    *  @return complementary distribution function evaluated at `x`
    */
   public double barF (double x);

   /**
    * Returns the inverse distribution function @f$F^{-1}(u)@f$, defined
    * in ( {@link REF_probdist_overview_eq_inverseF
    * inverseF} ).
    *  @param u            value in the interval @f$(0,1)@f$ for which the
    *                      inverse distribution function is evaluated
    *  @return the inverse distribution function evaluated at `u`
    */
   public double inverseF (double u);

   /**
    * Returns the mean of the distribution function.
    */
   public double getMean();

   /**
    * Returns the variance of the distribution function.
    */
   public double getVariance();

   /**
    * Returns the standard deviation of the distribution function.
    */
   public double getStandardDeviation();

   /**
    * Returns the parameters of the distribution function in the same
    * order as in the constructors.
    */
   public double[] getParams();

}