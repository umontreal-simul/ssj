/*
 * Class:        NormalInverseGaussianDist
 * Description:  normal inverse gaussian distribution
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
import umontreal.ssj.util.Num;

/**
 * Extends the class  @ref ContinuousDistribution for the *normal inverse
 * gaussian* distribution with location parameter @f$\mu@f$, scale parameter
 * @f$\delta> 0@f$, tail heavyness @f$\alpha> 0@f$, and asymmetry parameter
 * @f$\beta@f$ such that @f$0 \le|\beta| < \alpha@f$. Its density is
 * @anchor REF_probdist_NormalInverseGaussianDist_eq_fNormalInverseGaussian
 * @f[
 *   f(x) = \frac{\alpha\delta e^{\delta\gamma+ \beta(x - \mu)} K_1\left(\alpha\sqrt{\delta^2 + (x - \mu)^2}\right)}{\pi\sqrt{\delta^2 + (x - \mu)^2}}, \qquad\mbox{for } -\infty< x < \infty, \tag{fNormalInverseGaussian}
 * @f]
 * where @f$K_1@f$ is the modified Bessel function of the second kind of
 * order 1, and @f$\gamma= \sqrt{\alpha^2 - \beta^2}@f$.
 *
 * The distribution function is given by
 * @anchor REF_probdist_NormalInverseGaussianDist_eq_FNormalInverseGaussian
 * @f[
 *   F(x) = \int_{-\infty}^x dt f(t), \tag{FNormalInverseGaussian}
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class NormalInverseGaussianDist extends ContinuousDistribution {
   protected double alpha;
   protected double beta;
   protected double gamma;
   protected double delta;
   protected double mu;

   /**
    * Constructor for a *normal inverse gaussian* distribution with
    * parameters @f$\alpha@f$ = `alpha`, @f$\beta@f$ = `beta`,
    * @f$\mu@f$ = `mu` and @f$\delta@f$ = `delta`.
    */
   public NormalInverseGaussianDist (double alpha, double beta, double mu,
                                     double delta) {
      setParams (alpha, beta, mu, delta);
   }


   public double density (double x) {
      return density (alpha, beta, mu, delta, x);
   }

   public double cdf (double x) {
      return cdf (alpha, beta, mu, delta, x);
   }

   public double barF (double x) {
      return barF (alpha, beta, mu, delta, x);
   }

   public double getMean() {
      return getMean (alpha, beta, mu, delta);
   }

   public double getVariance() {
      return getVariance (alpha, beta, mu, delta);
   }

   public double getStandardDeviation() {
      return getStandardDeviation (alpha, beta, mu, delta);
   }

/**
 * Computes the density function (
 * {@link REF_probdist_NormalInverseGaussianDist_eq_fNormalInverseGaussian
 * fNormalInverseGaussian} ) for the *normal inverse gaussian* distribution
 * with parameters @f$\alpha@f$, @f$\beta@f$, @f$\mu@f$ and @f$\delta@f$,
 * evaluated at @f$x@f$.
 */
public static double density (double alpha, double beta, double mu,
                                 double delta, double x) {
      if (delta <= 0.0)
         throw new IllegalArgumentException ("delta <= 0");
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (Math.abs(beta) >= alpha)
         throw new IllegalArgumentException ("|beta| >= alpha");

      double gamma = Math.sqrt(alpha*alpha - beta*beta);
      double z = (x - mu)/delta;
      double w;
      if (Math.abs(z) <= 1.0e10)
         w = Math.sqrt (1.0 + z*z);
      else
         w = Math.abs(z);
      double y = alpha*delta*w;
      double v = delta*(gamma + beta*z);
      double R = Num.expBesselK1(v, y);
      return alpha * R / (Math.PI*w);
   }

   /**
    * NOT IMPLEMENTED. Computes the distribution function (
    * {@link REF_probdist_NormalInverseGaussianDist_eq_FNormalInverseGaussian
    * FNormalInverseGaussian} ) of the *normal inverse gaussian*
    * distribution with parameters @f$\alpha@f$, @f$\beta@f$, @f$\mu@f$
    * and @f$\delta@f$, evaluated at @f$x@f$.
    */
   public static double cdf (double alpha, double beta, double mu,
                             double delta, double x) {
      if (delta <= 0.0)
         throw new IllegalArgumentException ("delta <= 0");
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (Math.abs(beta) >= alpha)
         throw new IllegalArgumentException ("|beta| >= alpha");

      double gamma = Math.sqrt(alpha*alpha - beta*beta);
      double z = (x - mu)/delta;
      if (z > 0.0 && (gamma + (beta - alpha)*z >= XBIG))
         return 1.0;
      if (z < 0.0 && (gamma + (beta + alpha)*z <= -XBIGM))
         return 0.0;
 //     double w = Math.sqrt (1.0 + z*z);

      throw new UnsupportedOperationException
         ("NormalInverseGaussianDist:   cdf NOT IMPLEMENTED");
   }

   /**
    * NOT IMPLEMENTED. Computes the complementary distribution function of
    * the *normal inverse gaussian* distribution with parameters
    * @f$\alpha@f$, @f$\beta@f$, @f$\mu@f$ and @f$\delta@f$, evaluated
    * at @f$x@f$.
    */
   public static double barF (double alpha, double beta, double mu,
                              double delta, double x) {
      return 1.0 - cdf (alpha, beta, mu, delta, x);
   }

   /**
    * NOT IMPLEMENTED. Computes the inverse of the *normal inverse
    * gaussian* distribution with parameters @f$\alpha@f$, @f$\beta@f$,
    * @f$\mu@f$ and @f$\delta@f$.
    */
   public static double inverseF (double alpha, double beta, double mu,
                                  double delta, double u) {
      throw new UnsupportedOperationException(" Inversion NOT IMPLEMENTED");
   }

   /**
    * NOT IMPLEMENTED.
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{\alpha}@f$,
    * @f$\hat{\beta}@f$, @f$\hat{\mu}@f$, @f$\hat{\delta}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
/*
      double[] parameters = new double[4];
      double sum = 0;
      for (int i = 0; i < n; i++) {
         sum += x[i];
      }
      */
      throw new UnsupportedOperationException("getMLE is not implemented");

  //    return parameters;
   }

   /**
    * NOT IMPLEMENTED.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static NormalInverseGaussianDist getInstanceFromMLE (double[] x,
                                                               int n) {
      double parameters[] = getMLE (x, n);
      return new NormalInverseGaussianDist (parameters[0], parameters[1],
                                            parameters[2], parameters[3]);
   }

   /**
    * Returns the mean @f$E[X] = \mu+ \delta\beta/\gamma@f$ of the
    * *normal inverse gaussian* distribution with parameters
    * @f$\alpha@f$, @f$\beta@f$, @f$\mu@f$ and @f$\delta@f$.
    *  @return the mean of the normal inverse gaussian distribution
    * @f$E[X] = \mu+ \delta\beta/\gamma@f$
    */
   public static double getMean (double alpha, double beta, double mu,
                                 double delta) {
      if (delta <= 0.0)
         throw new IllegalArgumentException ("delta <= 0");
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (Math.abs(beta) >= alpha)
         throw new IllegalArgumentException ("|beta| >= alpha");

      double gamma = Math.sqrt(alpha*alpha - beta*beta);
      return mu + delta*beta/gamma;
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] =
    * \delta\alpha^2 / \gamma^3@f$ of the *normal inverse gaussian*
    * distribution with parameters @f$\alpha@f$, @f$\beta@f$, @f$\mu@f$
    * and @f$\delta@f$.
    *  @return the variance of the normal inverse gaussian distribution
    * @f$\mbox{Var}[X] = \delta\alpha^2 / \gamma^3@f$
    */
   public static double getVariance (double alpha, double beta, double mu,
                                     double delta) {
      if (delta <= 0.0)
         throw new IllegalArgumentException ("delta <= 0");
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (Math.abs(beta) >= alpha)
         throw new IllegalArgumentException ("|beta| >= alpha");

      double gamma = Math.sqrt(alpha*alpha - beta*beta);
      return delta*alpha*alpha / (gamma*gamma*gamma);
   }

   /**
    * Computes and returns the standard deviation of the *normal inverse
    * gaussian* distribution with parameters @f$\alpha@f$, @f$\beta@f$,
    * @f$\mu@f$ and @f$\delta@f$.
    *  @return the standard deviation of the normal inverse gaussian
    * distribution
    */
   public static double getStandardDeviation (double alpha, double beta,
                                              double mu, double delta) {
      return Math.sqrt (getVariance (alpha, beta, mu, delta));
   }

   /**
    * Returns the parameter @f$\alpha@f$ of this object.
    */
   public double getAlpha() {
      return alpha;
   }

   /**
    * Returns the parameter @f$\beta@f$ of this object.
    */
   public double getBeta() {
      return beta;
   }

   /**
    * Returns the parameter @f$\mu@f$ of this object.
    */
   public double getMu() {
      return mu;
   }

   /**
    * Returns the parameter @f$\delta@f$ of this object.
    */
   public double getDelta() {
      return delta;
   }

   /**
    * Sets the parameters @f$\alpha@f$, @f$\beta@f$, @f$\mu@f$ and
    * @f$\delta@f$ of this object.
    */
   public void setParams (double alpha, double beta, double mu,
                          double delta) {
      if (delta <= 0.0)
         throw new IllegalArgumentException ("delta <= 0");
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (Math.abs(beta) >= alpha)
         throw new IllegalArgumentException ("|beta| >= alpha");

      gamma = Math.sqrt(alpha*alpha - beta*beta);

      this.mu = mu;
      this.delta = delta;
      this.beta = beta;
      this.alpha = alpha;
   }

   /**
    * Returns a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$\alpha@f$,
    * @f$\beta@f$, @f$\mu@f$, @f$\delta@f$].
    */
   public double[] getParams () {
      double[] retour = {alpha, beta, mu, delta};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + ": alpha = " + alpha + ", beta = " + beta +
                  ", mu = " + mu + ", delta = " + delta;
   }

}