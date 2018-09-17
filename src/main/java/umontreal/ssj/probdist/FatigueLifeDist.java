/*
 * Class:        FatigueLifeDist
 * Description:  fatigue life distribution
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
import umontreal.ssj.probdist.NormalDist;
import optimization.*;

/**
 * Extends the class  @ref ContinuousDistribution for the *fatigue life*
 * distribution @cite tBIR69a&thinsp; with location parameter @f$\mu@f$,
 * scale parameter @f$\beta@f$ and shape parameter @f$\gamma@f$. Its
 * density is
 * @anchor REF_probdist_FatigueLifeDist_eq_fFatigueLife
 * @f[
 *   f(x) = \left[\frac{\sqrt{\frac{x - \mu}{\beta}} + \sqrt{\frac{\beta}{x - \mu}}}{2\gamma(x - \mu)}\right] \phi\left(\frac{\sqrt{\frac{x - \mu}{\beta}} - \sqrt{\frac{\beta}{x - \mu}}}{\gamma}\right), \qquad\mbox{for } x>\mu, \tag{fFatigueLife}
 * @f]
 * where @f$\phi@f$ is the probability density of the standard normal
 * distribution. The distribution function is given by
 * @anchor REF_probdist_FatigueLifeDist_eq_FFatigueLife
 * @f[
 *   F(x) = \Phi\left(\frac{\sqrt{\frac{x - \mu}{\beta}} - \sqrt{\frac{\beta}{x - \mu}}}{\gamma}\right), \qquad\mbox{for } x>\mu, \tag{FFatigueLife}
 * @f]
 * where @f$\Phi@f$ is the standard normal distribution function.
 * Restrictions: @f$\beta> 0@f$, @f$\gamma> 0@f$.
 *
 * The non-static versions of the methods `cdf`, `barF`, and `inverseF` call
 * the static version of the same name.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class FatigueLifeDist extends ContinuousDistribution {
   protected double mu;
   protected double beta;
   protected double gamma;

   private static class Optim implements Uncmin_methods
   {
      private int n;
      private double[] xi;
      private double mu;

      public Optim (double[] x, int n, double min)
      {
         this.n = n;
         this.mu = min;
         this.xi = new double[n];
         System.arraycopy (x, 0, this.xi, 0, n);
      }

      public double f_to_minimize (double[] p)
      {
         double sum = 0.0;

         if ((p[1] <= 0.0) || (p[2] <= 0.0))
            return 1e200;
      
         for (int i = 0; i < n; i++)
            sum -= Math.log (density (mu, p[1], p[2], xi[i]));

         return sum;
      }
      
      public void gradient (double[] x, double[] g)
      {
      }

      public void hessian (double[] x, double[][] h)
      {
      }
   }

   /**
    * Constructs a fatigue life distribution with parameters @f$\mu@f$,
    * @f$\beta@f$ and @f$\gamma@f$.
    */
   public FatigueLifeDist (double mu, double beta, double gamma) {
      setParams (mu, beta, gamma);
   }


   public double density (double x) {
      return FatigueLifeDist.density (mu, beta, gamma, x);
   }

   public double cdf (double x) {
      return FatigueLifeDist.cdf (mu, beta, gamma, x);
   }

   public double barF (double x) {
      return FatigueLifeDist.barF (mu, beta, gamma, x);
   }

   public double inverseF (double u) {
      return FatigueLifeDist.inverseF (mu, beta, gamma, u);
   }

   public double getMean() {
      return FatigueLifeDist.getMean (mu, beta, gamma);      
   }

   public double getVariance() {
      return FatigueLifeDist.getVariance (mu, beta, gamma);      
   }

   public double getStandardDeviation() {
      return FatigueLifeDist.getStandardDeviation (mu, beta, gamma);      
   }

/**
 * Computes the density (
 * {@link REF_probdist_FatigueLifeDist_eq_fFatigueLife
 * fFatigueLife} ) for the fatigue life distribution with parameters
 * @f$\mu@f$, @f$\beta@f$ and @f$\gamma@f$.
 */
public static double density (double mu, double beta, double gamma,
                                 double x) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (gamma <= 0.0)
         throw new IllegalArgumentException ("gamma <= 0");
      if (x <= mu)
         return 0.0;
      double y;
      y = (Math.sqrt ((x - mu) / beta) - Math.sqrt (beta / (x - mu))) / gamma;

      return (((Math.sqrt ((x - mu) / beta) + Math.sqrt (beta / (x - mu))) /
              (2 * gamma * (x - mu))) * NormalDist.density (0.0, 1.0, y));
   }

   /**
    * Computes the fatigue life distribution function with parameters
    * @f$\mu@f$, @f$\beta@f$ and @f$\gamma@f$.
    */
   public static double cdf (double mu, double beta, double gamma, double x) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (gamma <= 0.0)
         throw new IllegalArgumentException ("gamma <= 0");
      if (x <= mu)
         return 0.0;

      return NormalDist.cdf01 ((Math.sqrt ((x - mu) / beta) - Math.sqrt (beta / (x - mu))) / gamma);
   }

   /**
    * Computes the complementary distribution function of the fatigue life
    * distribution with parameters @f$\mu@f$, @f$\beta@f$ and
    * @f$\gamma@f$.
    */
   public static double barF (double mu, double beta, double gamma,
                              double x) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (gamma <= 0.0)
         throw new IllegalArgumentException ("gamma <= 0");
      if (x <= mu)
         return 1.0;

      return NormalDist.barF01 ((Math.sqrt ((x - mu) / beta) - Math.sqrt (beta / (x - mu))) / gamma);
   }

   /**
    * Computes the inverse of the fatigue life distribution with
    * parameters @f$\mu@f$, @f$\beta@f$ and @f$\gamma@f$.
    */
   public static double inverseF (double mu, double beta, double gamma,
                                  double u) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (gamma <= 0.0)
         throw new IllegalArgumentException ("gamma <= 0");
      if (u > 1.0 || u < 0.0)
          throw new IllegalArgumentException ("u not in [0,1]");
      if (u <= 0.0)    // if u == 0, in fact
          return mu;
      if (u >= 1.0)    // if u == 1, in fact
          return Double.POSITIVE_INFINITY;

      double w = gamma * NormalDist.inverseF01 (u);
      double sqrtZ = 0.5 * (w + Math.sqrt (w * w + 4.0));

      return (mu + sqrtZ * sqrtZ * beta);
   }

   /**
    * Estimates the parameters (@f$\mu@f$, @f$\beta@f$, @f$\gamma@f$)
    * of the fatigue life distribution using the maximum likelihood
    * method, from the @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1,…,
    * n-1@f$. The estimates are returned in a three-element array, in
    * regular order: [@f$\mu@f$, @f$\beta@f$, @f$\gamma@f$].  The
    * estimate of the parameters is given by maximizing numerically the
    * log-likelihood function, using the Uncmin package @cite iSCHa,
    * @cite iVERa&thinsp;.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @param mu           the location parameter
    *  @return returns the parameters [@f$\hat{\beta}@f$,
    * @f$\hat{\gamma}@f$]
    */
   public static double[] getMLE (double[] x, int n, double mu) {
      double sum = 0.0;
      
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double[] parameters = new double[3];
      double[] xpls = new double[3];
      double[] param = new double[3];
      double[] fpls = new double[3];
      double[] gpls = new double[3];
      int[] itrcmd = new int[2];
      double[][] h = new double[3][3];
      double[] udiag = new double[3];

      Optim system = new Optim (x, n, mu);

      double mean = 0.0;
      for (int i = 0; i < n; i++)
         mean += x[i];
      mean /= (double) n;

      double var = 0.0;
      for (int i = 0; i < n; i++)
         var += (x[i] - mean) * (x[i] - mean);
      var /= (double) n;

      double loc2 = (mean - mu) * (mean - mu);
      double a = 0.25 * (var - 5 * loc2);
      double b = (var - loc2);
      double c = var;

      double delta = b * b - 4.0 * a * c;

      double gamma2 = (- b - Math.sqrt (delta)) / (2.0 * a);
      param[2] = Math.sqrt (gamma2);
      param[1] = (mean - mu) / (1.0 + gamma2 / 2.0);

      Uncmin_f77.optif0_f77 (2, param, system, xpls, fpls, gpls, itrcmd, h, udiag);

      for (int i = 1; i < 3; i++)
         parameters[i] = xpls[i];
      parameters[0] = mu;
      return parameters;
   }

   /**
    * Computes and returns the mean @f$E[X] = \mu+ \beta(1 +
    * \gamma^2/2)@f$ of the fatigue life distribution with parameters
    * @f$\mu@f$, @f$\beta@f$ and @f$\gamma@f$.
    *  @return the mean of the fatigue life distribution
    */
   public static double getMean (double mu, double beta, double gamma) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (gamma <= 0.0)
         throw new IllegalArgumentException ("gamma <= 0");

      return (mu + beta * (1 + 0.5 * gamma * gamma));
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] = \beta^2
    * \gamma^2 (1 + 5 \gamma^2/4)@f$ of the fatigue life distribution
    * with parameters @f$\mu@f$, @f$\beta@f$ and @f$\gamma@f$.
    *  @return the variance of the fatigue life distribution
    */
   public static double getVariance (double mu, double beta, double gamma) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (gamma <= 0.0)
         throw new IllegalArgumentException ("gamma <= 0");

      return (beta * beta * gamma * gamma * (1.0 + 5.0/4.0 * gamma * gamma));
   }

   /**
    * Computes and returns the standard deviation of the fatigue life
    * distribution with parameters @f$\mu@f$, @f$\beta@f$ and
    * @f$\gamma@f$.
    *  @return the standard deviation of the fatigue life distribution
    */
   public static double getStandardDeviation (double mu, double beta,
                                              double gamma) {
      return Math.sqrt (FatigueLifeDist.getVariance (mu, beta, gamma));
   }

   /**
    * Returns the parameter @f$\beta@f$ of this object.
    */
   public double getBeta() {
      return beta;
   }

   /**
    * Returns the parameter @f$\gamma@f$ of this object.
    */
   public double getGamma() {
      return gamma;
   }

   /**
    * Returns the parameter @f$\mu@f$ of this object.
    */
   public double getMu() {
      return mu;
   }

   /**
    * Sets the parameters @f$\mu@f$, @f$\beta@f$ and @f$\gamma@f$ of
    * this object.
    */
   public void setParams (double mu, double beta, double gamma) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (gamma <= 0.0)
         throw new IllegalArgumentException ("gamma <= 0");
      
      this.mu = mu;
      this.beta = beta;
      this.gamma = gamma;
      supportA = mu;
   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$\mu@f$,
    * @f$\beta@f$, @f$\gamma@f$].
    */
   public double[] getParams () {
      double[] retour = {mu, beta, gamma};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : mu = " + mu + ", beta = " + beta + ", gamma = " + gamma;
   }

}