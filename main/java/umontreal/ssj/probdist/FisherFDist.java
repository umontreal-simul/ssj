/*
 * Class:        FisherFDist
 * Description:  Fisher F-distribution
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
import umontreal.ssj.probdist.BetaDist;
import umontreal.ssj.util.*;

/**
 * Extends the class  @ref ContinuousDistribution for the *Fisher F*
 * distribution with @f$n_1@f$ and @f$n_2@f$ degrees of freedom, where
 * @f$n_1@f$ and @f$n_2@f$ are positive integers. Its density is
 * @anchor REF_probdist_FisherFDist_eq_FisherF
 * @f[
 *   f(x) = \frac{\Gamma(\frac{n_1 + n_2}{2})n_1^{\frac{n_1}{2}}n_2^{\frac{n_2}{2}}}{\Gamma(\frac{n_1}{2})\Gamma(\frac{n_2}{2})} \frac{x^{\frac{n_1 - 2}{2}}}{(n_2 + n_1x)^{\frac{n_1 + n_2}{2}}}, \qquad\mbox{for } x > 0 \tag{FisherF}
 * @f]
 * where @f$\Gamma(x)@f$ is the gamma function defined in (
 * {@link REF_probdist_GammaDist_eq_Gamma Gamma} ).
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class FisherFDist extends ContinuousDistribution {
   protected int n1;
   protected int n2;
   protected double C1;
   private static final int DECPREC = 15;    // decimal precision

   /**
    * Constructs a Fisher @f$F@f$ distribution with `n1` and `n2` degrees
    * of freedom.
    */
   public FisherFDist (int n1, int n2) {
      setParams (n1, n2);
   }


   public double density (double x) {
      if (x <= 0.0)
         return 0.0;
      return Math.exp (C1 + 0.5 * (n1 - 2) * Math.log (x) -
           (0.5 * (n1 + n2) * Math.log (n2 + n1 * x)));
   }

   public double cdf (double x) {
      return FisherFDist.cdf (n1, n2, x);
   }

   public double barF (double x) {
      return FisherFDist.barF (n1, n2, x);
   }

   public double inverseF (double u) {
      return FisherFDist.inverseF (n1, n2, u);
   }

   public double getMean() {
      return FisherFDist.getMean (n1, n2);
   }

   public double getVariance() {
      return FisherFDist.getVariance (n1, n2);
   }

   public double getStandardDeviation() {
      return FisherFDist.getStandardDeviation (n1, n2);
   }

/**
 * Computes the density function (
 * {@link REF_probdist_FisherFDist_eq_FisherF FisherF} ) for a
 * Fisher @f$F@f$ distribution with `n1` and `n2` degrees of freedom,
 * evaluated at @f$x@f$.
 */
public static double density (int n1, int n2, double x) {
      if (n1 <= 0)
         throw new IllegalArgumentException ("n1 <= 0");
      if (n2 <= 0)
         throw new IllegalArgumentException ("n2 <= 0");
      if (x <= 0.0)
         return 0.0;

      return Math.exp (((n1/2.0) * Math.log (n1) + (n2/2.0) * Math.log(n2) +
          ((n1 - 2) / 2.0) * Math.log (x)) -
          (Num.lnBeta (n1/2.0, n2/2.0) + 
          ((n1 + n2) / 2.0) * Math.log (n2 + n1 * x)));
   }

   /**
    * Computes the distribution function of the Fisher @f$F@f$
    * distribution with parameters `n1` and `n2`, evaluated at @f$x@f$,
    * with roughly @f$d@f$ decimal digits of precision.
    */
   @Deprecated
   public static double cdf (int n1, int n2, int d, double x) {
      if (n1 <= 0)
         throw new IllegalArgumentException ("n1 <= 0");
      if (n2 <= 0)
         throw new IllegalArgumentException ("n2 <= 0");
      if (x <= 0.0)
         return 0.0;
      return BetaDist.cdf (n1/2.0, n2/2.0, d, (n1*x)/(n1*x + n2));
   }

   /**
    * Computes the distribution function of the Fisher @f$F@f$
    * distribution with parameters `n1` and `n2`, evaluated at @f$x@f$.
    */
   public static double cdf (int n1, int n2, double x) {
       return cdf (n1, n2, DECPREC, x);
   }

   /**
    * Computes the complementary distribution function of the Fisher
    * @f$F@f$ distribution with parameters `n1` and `n2`, evaluated at
    * @f$x@f$, with roughly @f$d@f$ decimal digits of precision.
    */
   @Deprecated
   public static double barF (int n1, int n2, int d, double x) {
      if (n1 <= 0)
         throw new IllegalArgumentException ("n1 <= 0");
      if (n2 <= 0)
         throw new IllegalArgumentException ("n2 <= 0");
      if (x <= 0.0)
         return 1.0;
      return BetaDist.barF (n1/2.0, n2/2.0, d, (n1 * x) / (n1 * x + n2));
   }

   /**
    * Computes the complementary distribution function of the Fisher
    * @f$F@f$ distribution with parameters `n1` and `n2`, evaluated at
    * @f$x@f$.
    */
   public static double barF (int n1, int n2, double x) {
       return barF (n1, n2, DECPREC, x);      
   }

   /**
    * Computes the inverse of the Fisher @f$F@f$ distribution with
    * parameters `n1` and `n2`, evaluated at @f$u@f$, with roughly @f$d@f$
    * decimal digits of precision.
    */
   @Deprecated
   public static double inverseF (int n1, int n2, int d, double u) {
      if (n1 <= 0)
         throw new IllegalArgumentException ("n1 <= 0");
      if (n2 <= 0)
         throw new IllegalArgumentException ("n2 <= 0");
      if (u > 1.0 || u < 0.0)
         throw new IllegalArgumentException ("u < 0 or u > 1");
      if (u <= 0.0)
         return 0.0;
      if (u >= 1.0)
         return Double.POSITIVE_INFINITY;

      double z = BetaDist.inverseF (n1 / 2.0, n2 / 2.0, d, u);
      return ((n2 * z) / (n1 * (1 - z)));
   }

   /**
    * Computes the inverse of the Fisher @f$F@f$ distribution with
    * parameters `n1` and `n2`, evaluated at @f$u@f$.
    */
   public static double inverseF (int n1, int n2, double u) {
       return inverseF (n1, n2, DECPREC, u);      
   }

   /**
    * Computes and returns the mean @f$E[X] = n_2 / (n_2 - 2)@f$ of the
    * Fisher @f$F@f$ distribution with parameters `n1` and `n2`
    * @f$=n_2@f$.
    *  @return the mean of the Fisher @f$F@f$ distribution
    */
   public static double getMean (int n1, int n2) {
      if (n1 <= 0)
         throw new IllegalArgumentException ("n1 <= 0");
      if (n2 <= 2)
         throw new IllegalArgumentException ("n2 <= 2");

      return (n2 / (n2 - 2.0));
   }

   /**
    * Computes and returns the variance
    * @f[
    *   \mbox{Var}[X] = \frac{2n_2^2 (n_2 + n_1 - 2)}{n_1 (n_2 - 2)^2 (n_2 - 4)}
    * @f]
    * of the Fisher @f$F@f$ distribution with parameters `n1` @f$=n_1@f$
    * and `n2` @f$=n_2@f$.
    *  @return the variance of the Fisher @f$F@f$ distribution
    * @f$\mbox{Var}[X] = (2n2^2 (n2 + n1 - 2)) / (n1 (n2 - 2)^2 (n2 -
    * 4))@f$
    */
   public static double getVariance (int n1, int n2) {
      if (n1 <= 0)
         throw new IllegalArgumentException ("n1 <= 0");
      if (n2 <= 4)
         throw new IllegalArgumentException ("n2 <= 4");

      return ((2.0 * n2 * n2 * (n2 + n1 - 2)) / (n1 * (n2 - 2.0) * (n2 - 2.0) * (n2 - 4.0)));
   }

   /**
    * Computes and returns the standard deviation of the Fisher @f$F@f$
    * distribution with parameters `n1` and `n2`.
    *  @return the standard deviation of the Fisher @f$F@f$ distribution
    */
   public static double getStandardDeviation (int n1, int n2) {
      return Math.sqrt (FisherFDist.getVariance (n1, n2));
   }

   /**
    * Returns the parameter `n1` of this object.
    */
   @Deprecated
   public int getN() {
      return n1;
   }

   /**
    * Returns the parameter `n2` of this object.
    */
   @Deprecated
   public int getM() {
      return n2;
   }

   /**
    * Returns the parameter `n1` of this object.
    */
   public int getN1() {
      return n1;
   }

   /**
    * Returns the parameter `n2` of this object.
    */
   public int getN2() {
      return n2;
   }

   /**
    * Sets the parameters `n1` and `n2` of this object.
    */
   public void setParams (int n1, int n2) {
      if (n1 <= 0)
         throw new IllegalArgumentException ("n1 <= 0");
      if (n2 <= 0)
         throw new IllegalArgumentException ("n2 <= 0");

      this.n1 = n1;
      this.n2 = n2;
      supportA = 0;
      C1 = (n1 / 2.0) * Math.log (n1) + (n2 / 2.0) * Math.log (n2) -
           Num.lnBeta (n1 / 2.0, n2 / 2.0);
   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [<tt>n1</tt>,
    * <tt>n2</tt>].
    */
   public double[] getParams () {
      double[] retour = {n1, n2};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : n1 = " + n1 + ", n2 = " + n2;
   }

}