/*
 * Class:        ContinuousDistribution2Dim
 * Description:  Mother class 2-dimensional continuous distributions
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
package umontreal.ssj.probdistmulti;
import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.util.Num;

/**
 * Classes implementing 2-dimensional continuous distributions should inherit
 * from this class. Such distributions are characterized by a *density*
 * function @f$f(x, y)@f$; thus the signature of a `density` method is
 * supplied here. This class also provides a default implementation of
 * @f$\overline{F}(x, y)@f$, the upper `CDF`. The inverse function
 * @f$F^{-1}(u)@f$ represents a curve @f$y = h(x)@f$ of constant @f$u@f$ and
 * it is not implemented.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdistmulti_general
 */
public abstract class ContinuousDistribution2Dim
                          extends ContinuousDistributionMulti {

/**
 * Defines the target number of decimals of accuracy when approximating a
 * distribution function, but there is *no guarantee* that this target is
 * always attained.
 */
public int decPrec = 15;


    // x infinity for some distributions
     protected static final double XINF = Double.MAX_VALUE;  

    // x infinity for some distributions                                       
    protected static final double XBIG = 1000.0;  

    // EPSARRAY[j]: Epsilon required for j decimal degits of precision
    protected static final double[] EPSARRAY = {
    0.5, 0.5E-1, 0.5E-2, 0.5E-3, 0.5E-4, 0.5E-5, 0.5E-6, 0.5E-7, 0.5E-8,
    0.5E-9, 0.5E-10, 0.5E-11, 0.5E-12, 0.5E-13, 0.5E-14, 0.5E-15, 0.5E-16,
    0.5E-17, 0.5E-18, 0.5E-19, 0.5E-20, 0.5E-21, 0.5E-22, 0.5E-23, 0.5E-24,
    0.5E-25, 0.5E-26, 0.5E-27, 0.5E-28, 0.5E-29, 0.5E-30, 0.5E-31, 0.5E-32,
    0.5E-33, 0.5E-34, 0.5E-35
    };

/**
 * Returns @f$f(x, y)@f$, the density of @f$(X, Y)@f$ evaluated at @f$(x,
 * y)@f$.
 *  @param x            value @f$x@f$ at which the density is evaluated
 *  @param y            value @f$y@f$ at which the density is evaluated
 *  @return density function evaluated at @f$(x, y)@f$
 */
public abstract double density (double x, double y);

   /**
    * Simply calls `density (x[0], x[1])`.
    *  @param x            point @f$(x[0], x[1])@f$ at which the density
    *                      is evaluated
    *  @return density function evaluated at @f$(x[0], x[1])@f$
    */
   public double density (double[] x) {
      if (x.length != 2)
         throw new IllegalArgumentException("x must be in dimension 2");

      return density (x[0], x[1]);
   }

   /**
    * Computes the distribution function @f$F(x, y)@f$:
    * @f[
    *   F(x, y) = P[X\le x, Y \le y] = \int_{-\infty}^x ds \int_{-\infty}^y dt  f(s, t).
    * @f]
    * @param x            value @f$x@f$ at which the distribution function
    *                      is evaluated
    *  @param y            value @f$y@f$ at which the distribution
    *                      function is evaluated
    *  @return distribution function evaluated at @f$(x, y)@f$
    */
   public abstract double cdf (double x, double y); 

   /**
    * Computes the upper cumulative distribution function
    * @f$\overline{F}(x, y)@f$:
    * @f[
    *   \overline{F}(x, y) = P[X\ge x, Y \ge y] = \int^{\infty}_x ds \int^{\infty}_y dt  f(s, t).
    * @f]
    * @param x            value @f$x@f$ at which the upper distribution is
    *                      evaluated
    *  @param y            value @f$y@f$ at which the upper distribution
    *                      is evaluated
    *  @return upper distribution function evaluated at @f$(x, y)@f$
    */
   public double barF (double x, double y) {
      double u = 1.0 + cdf (x, y) - cdf (XINF, y) - cdf (x, XINF);
      if (u <= 0.0) return 0.0;
      if (u >= 1.0) return 1.0;
      return u;
   }

   /**
    * Computes the cumulative probability in the square region
    * @f[
    *   P[a_1 \le X \le b_1,\: a_2 \le Y \le b_2] = \int_{a_1}^{b_1} dx \int_{a_2}^{b_2} dy  f(x, y).
    * @f]
    * @param a1           @f$x@f$ lower limit of the square
    *  @param a2           @f$y@f$ lower limit of the square
    *  @param b1           @f$x@f$ upper limit of the square
    *  @param b2           @f$y@f$ upper limit of the square
    *  @return the cumulative probability in the square region
    */
   public double cdf (double a1, double a2, double b1, double b2) {
      if (a1 >= b1 || a2 >= b2) return 0.0;
      return cdf (b1, b2) - cdf (a1, b2) - cdf (b1, a2) + cdf(a1, a2);
   }

}