/*
 * Class:        Distribution
 * Description:  interface for all discrete and continuous distributions
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