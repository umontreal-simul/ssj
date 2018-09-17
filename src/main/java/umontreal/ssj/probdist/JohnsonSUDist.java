/*
 * Class:        JohnsonSUDist
 * Description:  Johnson S_U distribution
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

/**
 * Extends the class  @ref ContinuousDistribution for the <em>Johnson
 * @f$S_U@f$</em> distribution (see @cite sLAW00a&thinsp; (page 316)). It has
 * shape parameters @f$\gamma@f$ and @f$\delta> 0@f$, location parameter
 * @f$\xi@f$, and scale parameter @f$\lambda> 0@f$. Denoting
 * @f$t=(x-\xi)/\lambda@f$ and @f$z = \gamma+ \delta\ln\left(t +
 * \sqrt{t^2 + 1}\right)@f$, the distribution has density
 * @f[
 *   f(x) = \frac{\delta e^{-z^2/2}}{\lambda\sqrt{2\pi(t^2 + 1)}},  \qquad\mbox{for } -\infty< x < \infty,
 * @f]
 * and distribution function
 * @f[
 *   F(x) = \Phi(z), \qquad\mbox{for } -\infty< x < \infty,
 * @f]
 * where @f$\Phi@f$ is the standard normal distribution function. The
 * inverse distribution function is
 * @f[
 *   F^{-1} (u) = \xi+ \lambda\sinh(v(u)), \qquad\mbox{for } 0 \le u \le1,
 * @f]
 * where
 * @f[
 *   v(u) = [\Phi^{-1}(u) - \gamma]/\delta.
 * @f]
 * This class relies on the methods  NormalDist.cdf01 and
 * NormalDist.inverseF01 of  @ref NormalDist to approximate @f$\Phi@f$ and
 * @f$\Phi^{-1}@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class JohnsonSUDist extends JohnsonSystem {

   private static double calcR (double a, double b, double x) {
      /*** cette fonction calcule
                 r = z + sqrt(z*z + 1)
           en utilisant un algorithme stable
       ***/

      double z = (x - a)/b;

      double s = Math.abs(z);
      if (s < 1.0e20)
         s = Math.sqrt (z * z + 1.0);

      // compute r = z + sqrt (z * z + 1)
      double r;
      if (z >= 0.0)
         r = s + z;
      else
         r = 1.0/(s - z);

      return r;
   }

   /**
    * Same as  {@link #JohnsonSUDist(double,double,double,double)
    * JohnsonSUDist(gamma, delta, 0, 1)}.
    */
   public JohnsonSUDist (double gamma, double delta) {
      this (gamma, delta, 0, 1);
   }

   /**
    * Constructs a `JohnsonSUDist` object with shape parameters
    * @f$\gamma@f$ and @f$\delta@f$, location parameter @f$\xi@f$, and
    * scale parameter @f$\lambda@f$.
    */
   public JohnsonSUDist (double gamma, double delta,
                         double xi, double lambda) {
      super (gamma, delta, xi, lambda);
   }
   public double density (double x) {
      return density (gamma, delta, xi, lambda, x);
   }

   public double cdf (double x) {
      return cdf (gamma, delta, xi, lambda, x);
   }

   public double barF (double x) {
      return barF (gamma, delta, xi, lambda, x);
   }

   public double inverseF (double u){
      return inverseF (gamma, delta, xi, lambda, u);
   }

   public double getMean() {
      return JohnsonSUDist.getMean (gamma, delta, xi, lambda);
   }

   public double getVariance() {
      return JohnsonSUDist.getVariance (gamma, delta, xi, lambda);
   }

   public double getStandardDeviation() {
      return JohnsonSUDist.getStandardDeviation (gamma, delta, xi, lambda);
   }

   /**
    * Returns the density function @f$f(x)@f$.
    */
   public static double density (double gamma, double delta,
                                 double xi, double lambda, double x) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (delta <= 0)
         throw new IllegalArgumentException ("delta <= 0");
      double r = calcR (xi, lambda, x);
      if (r <= 0.0)
         return 0.0;
      double z = gamma + delta*Math.log (r);
      double y = (x - xi)/lambda;
      if (z >= 1.0e10)
         return 0;
      return delta/(lambda*Math.sqrt (2.0*Math.PI)*Math.sqrt (y*y + 1.0))*
           Math.exp (-z*z/2.0);
   }

   /**
    * Returns the distribution function @f$F(x)@f$.
    */
   public static double cdf (double gamma, double delta,
                             double xi, double lambda, double x) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (delta <= 0)
         throw new IllegalArgumentException ("delta <= 0");
      double r = calcR (xi, lambda, x);
      if (r > 0.0)
         return NormalDist.cdf01 (gamma + delta*Math.log (r));
      else
         return 0.0;
   }

   /**
    * Returns the complementary distribution function @f$1-F(x)@f$.
    */
   public static double barF (double gamma, double delta,
                              double xi, double lambda, double x) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (delta <= 0)
         throw new IllegalArgumentException ("delta <= 0");

      double r = calcR (xi, lambda, x);
      if (r > 0.0)
         return NormalDist.barF01 (gamma + delta * Math.log (r));
      else
         return 1.0;
   }

   /**
    * Returns the inverse distribution function @f$F^{-1}(u)@f$.
    */
   public static double inverseF (double gamma, double delta,
                                  double xi, double lambda, double u) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (delta <= 0)
         throw new IllegalArgumentException ("delta <= 0");
      if (u > 1.0 || u < 0.0)
          throw new IllegalArgumentException ("u not in [0,1]");

      if (u >= 1.0)
         return Double.POSITIVE_INFINITY;
      if (u <= 0.0)
         return Double.NEGATIVE_INFINITY;

      double z = NormalDist.inverseF01 (u);
      double v = (z - gamma)/delta;
      if (v >= Num.DBL_MAX_EXP*Num.LN2)
         return Double.POSITIVE_INFINITY;
      if (v <= Num.LN2*Num.DBL_MIN_EXP)
         return Double.NEGATIVE_INFINITY;

      return xi + lambda * Math.sinh(v);
   }

   /**
    * Returns the mean
    * @f[
    *   E[X] = \xi- \lambda e^{1/(2\delta^2)} \sinh({\gamma}/{\delta})
    * @f]
    * of the Johnson @f$S_U@f$ distribution with parameters @f$\gamma@f$,
    * @f$\delta@f$, @f$\xi@f$ and @f$\lambda@f$.
    *  @return the mean of the Johnson @f$S_U@f$ distribution @f$E[X] =
    * \xi- \lambda\exp^{1 / (2\delta^2)} sinh(\gamma/ \delta)@f$
    */
   public static double getMean (double gamma, double delta,
                                 double xi, double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (delta <= 0.0)
         throw new IllegalArgumentException ("delta <= 0");

      return (xi - (lambda * Math.exp(1.0 / (2.0 * delta * delta)) *
                             Math.sinh(gamma / delta)));
   }

   /**
    * Returns the variance
    * @f[
    *   \mbox{Var}[X] = \frac{\lambda^2}{2} \left(e^{1/\delta^2} - 1\right)\left(e^{1/\delta^2} \cosh(2 {\gamma}/{\delta}) + 1\right)
    * @f]
    * of the Johnson @f$S_U@f$ distribution with parameters @f$\gamma@f$,
    * @f$\delta@f$, @f$\xi@f$ and @f$\lambda@f$.
    *  @return the variance of the Johnson @f$S_U@f$ distribution
    * @f$\mbox{Var}[X] = (\lambda^2/2)(\exp^{1/\delta^2} -
    * 1)(\exp^{1/\delta^2} cosh(2 \gamma/ \delta) + 1)@f$
    */
   public static double getVariance (double gamma, double delta,
                                     double xi, double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (delta <= 0.0)
         throw new IllegalArgumentException ("delta <= 0");

      double omega2 = Math.exp(1 / (delta * delta));
      return ((omega2 - 1) * (omega2 * Math.cosh(2 * gamma / delta) + 1) / 2 * lambda * lambda);
   }

   /**
    * Returns the standard deviation of the Johnson @f$S_U@f$ distribution
    * with parameters @f$\gamma@f$, @f$\delta@f$, @f$\xi@f$,
    * @f$\lambda@f$.
    *  @return the standard deviation of the Johnson @f$S_U@f$
    * distribution
    */
   public static double getStandardDeviation (double gamma, double delta,
                                              double xi, double lambda) {
      return Math.sqrt (JohnsonSUDist.getVariance (gamma, delta, xi, lambda));
   }

   /**
    * Sets the value of the parameters @f$\gamma@f$, @f$\delta@f$,
    * @f$\xi@f$ and @f$\lambda@f$ for this object.
    */
   public void setParams (double gamma, double delta,
                          double xi, double lambda) {
      setParams0(gamma, delta, xi, lambda);
   }

}