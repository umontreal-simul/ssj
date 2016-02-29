/*
 * Class:        LognormalDist
 * Description:  lognormal distribution
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
import umontreal.ssj.util.*;

/**
 * Extends the class  @ref ContinuousDistribution for the *lognormal*
 * distribution @cite tJOH95a&thinsp;. (See also the <em>Johnson
 * @f$S_L@f$</em> distribution `JohnsonSLDist` in this package.) It has scale
 * parameter @f$\mu@f$ and shape parameter @f$\sigma> 0@f$. The density is
 * @anchor REF_probdist_LognormalDist_eq_flognormal
 * @f[
 *   f(x) = \frac{1}{\sqrt{2\pi}\sigma x} e^{-(\ln(x) - \mu)^2/(2\sigma^2)} \qquad\mbox{for } x>0, \tag{flognormal}
 * @f]
 * and 0 elsewhere. The distribution function is
 * @f[
 *   F(x) = \Phi\left({(\ln(x) - \mu)/\sigma}\right) \qquad\mbox{for } x>0,
 * @f]
 * where @f$\Phi@f$ is the standard normal distribution function. Its
 * inverse is given by
 * @f[
 *   F^{-1}(u) = e^{\mu+ \sigma\Phi^{-1} (u)} \qquad\mbox{for } 0 \le u < 1.
 * @f]
 * If @f$\ln(Y)@f$ has a *normal* distribution, then @f$Y@f$ has a
 * *lognormal* distribution with the same parameters.
 *
 * This class relies on the methods  NormalDist.cdf01 and
 * NormalDist.inverseF01 of  @ref NormalDist to approximate @f$\Phi@f$ and
 * @f$\Phi^{-1}@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class LognormalDist extends ContinuousDistribution {
   private double mu;
   private double sigma;

   /**
    * Constructs a `LognormalDist` object with default parameters @f$\mu=
    * 0@f$ and @f$\sigma= 1@f$.
    */
   public LognormalDist() {
      setParams (0.0, 1.0);
   }

   /**
    * Constructs a `LognormalDist` object with parameters @f$\mu@f$ =
    * `mu` and @f$\sigma@f$ = `sigma`.
    */
   public LognormalDist (double mu, double sigma) {
      setParams (mu, sigma);
   }


   public double density (double x) {
      return density (mu, sigma, x);
   }

   public double cdf (double x) {
      return cdf (mu, sigma, x);
   }

   public double barF (double x) {
      return barF (mu, sigma, x);
   }

   public double inverseF (double u) {
      return inverseF (mu, sigma, u);
   }

   public double getMean() {
      return LognormalDist.getMean (mu, sigma);
   }

   public double getVariance() {
      return LognormalDist.getVariance (mu, sigma);
   }

   public double getStandardDeviation() {
      return LognormalDist.getStandardDeviation (mu, sigma);
   }

/**
 * Computes the lognormal density function @f$f(x)@f$ in (
 * {@link REF_probdist_LognormalDist_eq_flognormal flognormal}
 * ).
 */
public static double density (double mu, double sigma, double x) {
      if (sigma <= 0)
         throw new IllegalArgumentException ("sigma <= 0");
      if (x <= 0)
         return 0;
      double diff = Math.log (x) - mu;
      return Math.exp (-diff*diff/(2*sigma*sigma))/
              (Math.sqrt (2*Math.PI)*sigma*x);
   }

   /**
    * Computes the lognormal distribution function, using
    * NormalDist.cdf01.
    */
   public static double cdf (double mu, double sigma, double x) {
      if (sigma <= 0.0)
        throw new IllegalArgumentException ("sigma  <= 0");
      if (x <= 0.0)
         return 0.0;
      return NormalDist.cdf01 ((Math.log (x) - mu)/sigma);
   }

   /**
    * Computes the lognormal complementary distribution function
    * @f$\bar{F}(x)@f$, using  NormalDist.barF01.
    */
   public static double barF (double mu, double sigma, double x) {
      if (sigma <= 0.0)
        throw new IllegalArgumentException ("sigma  <= 0");
      if (x <= 0.0)
         return 1.0;
      return NormalDist.barF01 ((Math.log (x) - mu)/sigma);
   }

   /**
    * Computes the inverse of the lognormal distribution function, using
    * NormalDist.inverseF01.
    */
   public static double inverseF (double mu, double sigma, double u) {
        double t, v;

        if (sigma <= 0.0)
            throw new IllegalArgumentException ("sigma  <= 0");

        if (u > 1.0 || u < 0.0)
            throw new IllegalArgumentException ("u not in [0,1]");

        if (Num.DBL_EPSILON >= 1.0 - u)
            return Double.POSITIVE_INFINITY;

        if (u <= 0.0)
            return 0.0;

        t = NormalDist.inverseF01 (u);
        v = mu + sigma * t;

        if ((t >= XBIG) || (v >= Num.DBL_MAX_EXP * Num.LN2))
            return Double.POSITIVE_INFINITY;
        if ((t <= -XBIG) || (v <= -Num.DBL_MAX_EXP*Num.LN2))
            return 0.0;

        return Math.exp (v);
   }

   /**
    * Estimates the parameters @f$(\mu, \sigma)@f$ of the lognormal
    * distribution using the maximum likelihood method, from the @f$n@f$
    * observations @f$x[i]@f$, @f$i = 0, 1,…, n-1@f$. The estimates are
    * returned in a two-element array, in regular order: [@f$\mu@f$,
    * @f$\sigma@f$].  The maximum likelihood estimators are the values
    * @f$(\hat{\mu}, \hat{\sigma})@f$ that satisfy the equations:
    * @f{align*}{
    *    \hat{\mu} 
    *    & 
    *    = 
    *    \frac{1}{n} \sum_{i=1}^n \ln(x_i)
    *    \\ 
    *    \hat{\sigma} 
    *    & 
    *    = 
    *    \sqrt{\frac{1}{n} \sum_{i=1}^n (\ln(x_i) - \hat{\mu})^2}.
    * @f}
    * See @cite tJOH95a&thinsp; (page 220).
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{\mu}@f$,
    * @f$\hat{\sigma}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      final double LN_EPS = Num.LN_DBL_MIN - Num.LN2;
      double parameters[];
      parameters = new double[2];
      double sum = 0.0;
      for (int i = 0; i < n; i++) {
         if (x[i] > 0.0)
            sum += Math.log (x[i]);
         else
            sum += LN_EPS;       // log(DBL_MIN / 2)
      }
      parameters[0] = sum / n;

      double temp;
      sum = 0.0;
      for (int i = 0; i < n; i++) {
         if (x[i] > 0.0)
            temp = Math.log (x[i]) - parameters[0];
         else
            temp = LN_EPS - parameters[0];
         sum += temp * temp;
      }
      parameters[1] = Math.sqrt (sum / n);

      return parameters;
   }

   /**
    * Creates a new instance of a lognormal distribution with parameters
    * @f$\mu@f$ and @f$\sigma@f$ estimated using the maximum likelihood
    * method based on the @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1, …,
    * n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static LognormalDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new LognormalDist (parameters[0], parameters[1]);
   }

   /**
    * Computes and returns the mean @f$E[X] = e^{\mu+ \sigma^2/2}@f$ of
    * the lognormal distribution with parameters @f$\mu@f$ and
    * @f$\sigma@f$.
    *  @return the mean of the lognormal distribution
    */
   public static double getMean (double mu, double sigma) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");

      return (Math.exp(mu + (sigma * sigma) / 2.0));
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] = e^{2\mu+
    * \sigma^2}(e^{\sigma^2} - 1)@f$ of the lognormal distribution with
    * parameters @f$\mu@f$ and @f$\sigma@f$.
    *  @return the variance of the lognormal distribution
    */
   public static double getVariance (double mu, double sigma) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");

      return (Math.exp(2.0 * mu + sigma * sigma) * (Math.exp(sigma * sigma) - 1.0));
   }

   /**
    * Computes and returns the standard deviation of the lognormal
    * distribution with parameters @f$\mu@f$ and @f$\sigma@f$.
    *  @return the standard deviation of the lognormal distribution
    */
   public static double getStandardDeviation (double mu, double sigma) {
      return Math.sqrt (LognormalDist.getVariance (mu, sigma));
   }

   /**
    * Returns the parameter @f$\mu@f$ of this object.
    */
   public double getMu() {
      return mu;
   }

   /**
    * Returns the parameter @f$\sigma@f$ of this object.
    */
   public double getSigma() {
      return sigma;
   }

   /**
    * Sets the parameters @f$\mu@f$ and @f$\sigma@f$ of this object.
    */
   public void setParams (double mu, double sigma) {
      if (sigma <= 0)
         throw new IllegalArgumentException ("sigma <= 0");
      this.mu = mu;
      this.sigma = sigma;
      supportA = 0.0;
   }

   /**
    * Returns a table containing the parameters of the current
    * distribution, in the order: [@f$\mu@f$, @f$\sigma@f$].
    */
   public double[] getParams () {
      double[] retour = {mu, sigma};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : mu = " + mu + ", sigma = " + sigma;
   }

}