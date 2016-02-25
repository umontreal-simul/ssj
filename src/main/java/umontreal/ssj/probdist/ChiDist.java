/*
 * Class:        ChiDist
 * Description:  chi distribution
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
import umontreal.ssj.util.Num;
import umontreal.ssj.functions.MathFunction;

/**
 * Extends the class  @ref ContinuousDistribution for the *chi* distribution
 * @cite tJOH95a&thinsp; (page 417) with shape parameter @f$\nu > 0@f$,
 * where the number of degrees of freedom @f$\nu@f$ is a positive integer.
 * The density function is given by
 * @anchor REF_probdist_ChiDist_eq_Fchi
 * @f[
 *   f (x) = \frac{e^{-x^2 /2} x^{\nu-1}}{2^{(\nu/2)-1}\Gamma(\nu/2)}, \qquad\mbox{ for } x > 0, \tag{Fchi}
 * @f]
 * where @f$\Gamma(x)@f$ is the gamma function defined in (
 * {@link REF_probdist_GammaDist_eq_Gamma Gamma} ). The
 * distribution function is
 * @f[
 *   F (x) = \frac{1}{\Gamma(\nu/2)} \int_0^{x^2/2} t^{\nu/2-1}e^{-t} dt.
 * @f]
 * It is equivalent to the gamma distribution function with parameters
 * @f$\alpha=\nu/2@f$ and @f$\lambda=1@f$, evaluated at @f$x^2/2@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class ChiDist extends ContinuousDistribution {
   private int nu;
   private double C1;

   private static class Function implements MathFunction {
      protected int n;
      protected double sum;

      public Function (double s, int n)
      {
         this.n = n;
         this.sum = s;
      }

      public double evaluate (double k)
      {
         if (k < 1.0) return 1.0e200;
         return (sum + n * (Num.lnGamma (k / 2.0) - 0.5*(Num.LN2) - Num.lnGamma ((k + 1.0) / 2.0)));
      }
   }

   /**
    * Constructs a `ChiDist` object.
    */
   public ChiDist (int nu) {
      setNu (nu);
   }


   public double density (double x) {
       if (x <= 0.0)
         return 0.0;
      return Math.exp ((nu - 1)*Math.log (x) - x*x/2.0 - C1);
   }

   public double cdf (double x) {
      return cdf (nu, x);
   }

   public double barF (double x) {
      return barF (nu, x);
   }

   public double inverseF (double u) {
      return inverseF (nu, u);
   }

   public double getMean() {
      return ChiDist.getMean (nu);
   }

   public double getVariance() {
      return ChiDist.getVariance (nu);
   }

   public double getStandardDeviation() {
      return ChiDist.getStandardDeviation (nu);
   }

/**
 * Computes the density function.
 */
public static double density (int nu, double x) {
      if (nu <= 0)
         throw new IllegalArgumentException ("nu <= 0");
      if (x <= 0.0)
         return 0.0;
      return Math.exp ((nu - 1)*Math.log (x) - x*x/2.0
                         - (nu/2.0 - 1.0)*Num.LN2 - Num.lnGamma (nu/2.0));
   }

   /**
    * Computes the distribution function by using the gamma distribution
    * function.
    */
   public static double cdf (int nu, double x) {
      if (x <= 0.0)
         return 0.0;
      return GammaDist.cdf (nu/2.0, 15, x*x/2.0);
   }

   /**
    * Computes the complementary distribution.
    */
   public static double barF (int nu, double x) {
      if (x <= 0.0)
         return 1.0;
      return GammaDist.barF (nu/2.0, 15, x*x/2.0);
   }

   /**
    * Returns the inverse distribution function computed using the gamma
    * inversion.
    */
   public static double inverseF (int nu, double u) {
       double res =  GammaDist.inverseF (nu/2.0, 15, u);
       return Math.sqrt (2*res);
   }

   /**
    * Estimates the parameter @f$\nu@f$ of the chi distribution using the
    * maximum likelihood method, from the @f$n@f$ observations @f$x[i]@f$,
    * @f$i = 0, 1, …, n-1@f$. The estimate is returned in element 0 of the
    * returned array.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @return returns the parameter [@f$\hat{\nu}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      double[] parameters = new double[1];

      double mean = 0.0;
      for (int i = 0; i < n; i++)
         mean += x[i];
      mean /= (double) n;

      double var = 0.0;
      for (int i = 0; i < n; i++)
         var += ((x[i] - mean) * (x[i] - mean));
      var /= (double) n;

      double k = Math.round (var + mean * mean) - 5.0;
      if (k < 1.0)
         k = 1.0;

      double sum = 0.0;
      for (int i = 0; i < n; i++) {
         if (x[i] > 0.0)
            sum += Math.log (x[i]);
         else
            sum -= 709.0;
      }

      Function f = new Function (sum, n);
      while (f.evaluate(k) > 0.0)
         k++;
      parameters[0] = k;

      return parameters;
   }

   /**
    * Creates a new instance of a chi distribution with parameter
    * @f$\nu@f$ estimated using the maximum likelihood method based on
    * the @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static ChiDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new ChiDist ((int) parameters[0]);
   }

   /**
    * Computes and returns the mean
    * @f[
    *   E[X] = \frac{\sqrt{2} \Gamma( \frac{\nu+ 1}{2} )}{\Gamma(\frac{\nu}{2})}
    * @f]
    * of the chi distribution with parameter @f$\nu@f$.
    *  @return the mean of the chi distribution @f$E[X] =
    * \sqrt{2}\Gamma((\nu+ 1) / 2) / \Gamma(\nu/ 2)@f$
    */
   public static double getMean (int nu) {
      if (nu <= 0)
         throw new IllegalArgumentException ("nu <= 0");
      return  Num.RAC2 * Num.gammaRatioHalf(nu / 2.0);
   }

   /**
    * Computes and returns the variance
    * @f[
    *   \mbox{Var}[X] = \frac{2 \Gamma(\frac{\nu}{2}) \Gamma(1 + \frac{\nu}{2}) - \Gamma^2(\frac{\nu+ 1}{2})}{\Gamma(\frac{\nu}{2})}
    * @f]
    * of the chi distribution with parameter @f$\nu@f$.
    *  @return the variance of the chi distribution @f$\mbox{Var}[X] = 2 [
    * \Gamma(\nu/ 2) \Gamma(1 + \nu/ 2) - \Gamma^2(1/2 (\nu+ 1)) ] /
    * \Gamma(\nu/ 2)@f$
    */
   public static double getVariance (int nu) {
      if (nu <= 0)
         throw new IllegalArgumentException ("nu <= 0");
      double mean = ChiDist.getMean(nu);
      return (nu - (mean * mean));
   }

   /**
    * Computes and returns the standard deviation of the chi distribution
    * with parameter @f$\nu@f$.
    *  @return the standard deviation of the chi distribution
    */
   public static double getStandardDeviation (int nu) {
      return Math.sqrt (ChiDist.getVariance (nu));
   }

   /**
    * Returns the value of @f$\nu@f$ for this object.
    */
   public int getNu() {
      return nu;
   }

   /**
    * Sets the value of @f$\nu@f$ for this object.
    */
   public void setNu (int nu) {
      if (nu <= 0)
         throw new IllegalArgumentException ("nu <= 0");
      this.nu = nu;
      supportA = 0.0;
      C1 = (nu/2.0 - 1.0)*Num.LN2 + Num.lnGamma (nu/2.0);
   }

   /**
    * Return a table containing parameters of the current distribution.
    */
   public double[] getParams () {
      double[] retour = {nu};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : nu = " + nu;
   }

}