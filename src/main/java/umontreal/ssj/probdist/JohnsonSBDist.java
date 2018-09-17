/*
 * Class:        JohnsonSBDist
 * Description:  Johnson S_B distribution
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
 * @f$S_B@f$</em> distribution @cite tJOH49a, @cite sLAW00a,
 * @cite tFLY06a&thinsp; with shape parameters @f$\gamma@f$ and @f$\delta>
 * 0@f$, location parameter @f$\xi@f$, and scale parameter @f$\lambda>0@f$.
 * Denoting @f$t=(x-\xi)/\lambda@f$ and @f$z = \gamma+
 * \delta\ln(t/(1-t))@f$, the density is
 * @anchor REF_probdist_JohnsonSBDist_eq_JohnsonSB_density
 * @f[
 *   f(x) = \frac{\delta e^{-z^2/2}}{\lambda t(1 - t)\sqrt{2\pi}},  \qquad\mbox{ for } \xi< x < \xi+\lambda, \tag{JohnsonSB-density}
 * @f]
 * and 0 elsewhere. The distribution function is
 * @anchor REF_probdist_JohnsonSBDist_eq_JohnsonSB_dist
 * @f[
 *   F(x) = \Phi(z), \mbox{ for } \xi< x < \xi+\lambda, \tag{JohnsonSB-dist}
 * @f]
 * where @f$\Phi@f$ is the standard normal distribution function. The
 * inverse distribution function is
 * @anchor REF_probdist_JohnsonSBDist_eq_JohnsonSB_inverse
 * @f[
 *   F^{-1}(u) = \xi+ \lambda\left(1/\left(1+e^{-v(u)}\right)\right) \qquad\mbox{for } 0 \le u \le1, \tag{JohnsonSB-inverse}
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
public class JohnsonSBDist extends JohnsonSystem {
   // m_psi is used in computing the mean and the variance
   private double m_psi = -1.0e100;


   private static double getMeanPsi (double gamma, double delta,
                           double xi, double lambda, double[] tpsi) {
      // Returns the theoretical mean of t = (x - xi)/lambda;
      // also compute psi and returns it in tpsi[0], since
      // it is used in computing the mean and the variance

      final int NMAX = 10000;
      final double EPS = 1.0e-15;

      double a1 = 1.0/(2*delta*delta);
      double a2 = (1.0 - 2*delta*gamma)/(2*delta*delta);
      double a3 = (gamma - 1./delta)/delta;
      int n = 0;
      double tem = 1;
      double sum = 0;
      double v;
      while (Math.abs(tem) > EPS* Math.abs(sum) && n < NMAX) {
         ++n;
         v = Math.exp(-n*gamma/delta) + Math.exp(n*a3);
         tem = Math.exp(-n*n*a1) * v / (1 + Math.exp(-2*n*a1));
      //   tem = Math.exp(-n*n*a1) * Math.cosh(n*a2) / Math.cosh(n*a1);
         sum += tem;
      }
      if (n >= NMAX)
         System.err.println ("JohnsonSBDist:  possible lack of accuracy on mean");
      double A = (0.5 + sum) / (delta);

      a1 = Math.PI * Math.PI * delta * delta;
      a2 = Math.PI * delta * gamma;
      int j;
      n = 0;
      tem = 1;
      sum = 0;
      while (Math.abs(tem) > EPS*Math.abs(sum) && n < NMAX) {
         ++n;
         j = 2*n - 1;
         tem = Math.exp(-j*j*a1/2.0) * Math.sin(j*a2) / Math.sinh(j*a1);
         sum += tem;
      }
      if (n >= NMAX)
         System.err.println ("JohnsonSBDist:  possible lack of accuracy on mean");
      double B = 2.0* Math.PI * delta * sum;

      a1 = 2*Math.PI * Math.PI * delta * delta;
      a2 = 2*Math.PI * delta * gamma;
      n = 0;
      tem = 1;
      sum = 0;
      while (Math.abs(tem) > EPS*Math.abs(sum) && n < NMAX) {
         ++n;
         tem = Math.exp(-n*n*a1) * Math.cos(n*a2);
         sum += tem;
      }
      if (n >= NMAX)
         System.err.println ("JohnsonSBDist:  possible lack of accuracy on mean");
      double C = 1 + 2.0 * sum;

      double D = Math.sqrt(2*Math.PI) * Math.exp(gamma* gamma / 2.0);
      double tmean = (A - B) / (C*D);
      tpsi[0] = C*D;
      return tmean;
   }

   /**
    * Constructs a `JohnsonSBDist` object with shape parameters
    * @f$\gamma@f$ and @f$\delta@f$, location parameter @f$\xi@f$ and
    * scale parameter @f$\lambda@f$.
    */
   public JohnsonSBDist (double gamma, double delta,
                         double xi, double lambda) {
      super (gamma, delta, xi, lambda);
      setLastParams(xi, lambda);
   }
   private void setLastParams(double xi, double lambda) {
      supportA = xi;
      supportB = xi + lambda;
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

   public double inverseF (double u) {
      return inverseF (gamma, delta, xi, lambda, u);
   }

   public double getMean() {
      return JohnsonSBDist.getMean (gamma, delta, xi, lambda);
   }

   public double getVariance() {
      return JohnsonSBDist.getVariance (gamma, delta, xi, lambda);
   }

   public double getStandardDeviation() {
      return JohnsonSBDist.getStandardDeviation (gamma, delta, xi, lambda);
   }

   /**
    * Returns the density function (
    * {@link REF_probdist_JohnsonSBDist_eq_JohnsonSB_density
    * JohnsonSB-density} ).
    */
   public static double density (double gamma, double delta,
                                 double xi, double lambda, double x) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (delta <= 0)
         throw new IllegalArgumentException ("delta <= 0");
      if (x <= xi || x >= (xi+lambda))
         return 0.0;
      double y = (x - xi)/lambda;
      double z = gamma + delta*Math.log (y/(1.0 - y));
      return delta/(lambda*y*(1.0 - y)*Math.sqrt (2.0*Math.PI))*
           Math.exp (-z*z/2.0);
   }

   /**
    * Returns the distribution function (
    * {@link REF_probdist_JohnsonSBDist_eq_JohnsonSB_dist
    * JohnsonSB-dist} ).
    */
   public static double cdf (double gamma, double delta,
                             double xi, double lambda, double x) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (delta <= 0)
         throw new IllegalArgumentException ("delta <= 0");
      if (x <= xi)
         return 0.0;
      if (x >= xi+lambda)
         return 1.0;
      double y = (x - xi)/lambda;
      double z = gamma + delta*Math.log (y/(1.0 - y));
      return NormalDist.cdf01 (z);
   }

   /**
    * Returns the complementary distribution.
    */
   public static double barF (double gamma, double delta,
                              double xi, double lambda, double x) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (delta <= 0)
         throw new IllegalArgumentException ("delta <= 0");
      if (x <= xi)
         return 1.0;
      if (x >= xi+lambda)
         return 0.0;
      double y = (x - xi)/lambda;
      double z = gamma + delta*Math.log (y/(1.0 - y));
      return NormalDist.barF01 (z);
   }

   /**
    * Returns the inverse of the distribution (
    * {@link REF_probdist_JohnsonSBDist_eq_JohnsonSB_inverse
    * JohnsonSB-inverse} ).
    */
   public static double inverseF (double gamma, double delta,
                                  double xi, double lambda, double u) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (delta <= 0)
         throw new IllegalArgumentException ("delta <= 0");
      if (u > 1.0 || u < 0.0)
          throw new IllegalArgumentException ("u not in [0,1]");

      if (u >= 1.0)    // if u == 1, in fact
          return xi+lambda;
      if (u <= 0.0)    // if u == 0, in fact
          return xi;

      double z = NormalDist.inverseF01 (u);
      double v = (z - gamma)/delta;

      if (v >= Num.DBL_MAX_EXP*Num.LN2)
            return xi + lambda;
      if (v <= Num.DBL_MIN_EXP*Num.LN2)
            return xi;

      v = Math.exp (v);
      return (xi + (xi+lambda)*v)/(1.0 + v);
   }

   /**
    * Estimates the parameters @f$(\gamma,\delta)@f$ of the Johnson
    * @f$S_B@f$ distribution, using the maximum likelihood method, from
    * the @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1,…, n-1@f$.
    * Parameters @f$\xi= \mathtt{xi}@f$ and @f$\lambda=
    * \mathtt{lambda}@f$ are known. The estimated parameters are returned
    * in a two-element array in the order: [@f$\gamma@f$, @f$\delta@f$].
    * The maximum likelihood estimators are the values @f$(\hat{\gamma},
    * \hat{\delta})@f$ that satisfy the equations @cite tFLY06a&thinsp;
    * @f$\hat{\gamma}= -\bar{f} / s_f@f$ and @f$\hat{\delta}= 1/s_f@f$,
    * where @f$f = \ln(t/(1-t))@f$, @f$\bar{f}@f$ is the sample mean of
    * the @f$f_i@f$, and
    * @f[
    *   s_f = \sqrt{\frac{1}{n} \sum_{i=0}^{n-1} (f_i - \bar{f})^2},
    * @f]
    * with @f$f_i = \ln(t_i/(1-t_i))@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @param xi           parameter @f$\xi@f$
    *  @param lambda       parameter @f$\lambda@f$
    *  @return returns the parameters [@f$\hat{\gamma}@f$,
    * @f$\hat{\delta}@f$]
    */
   public static double[] getMLE (double[] x, int n,
                                  double xi, double lambda) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      final double LN_EPS = Num.LN_DBL_MIN - Num.LN2;
      double[] ftab = new double[n];
      double sum = 0.0;
      double t;

      for (int i = 0; i < n; i++) {
         t = (x[i] - xi) / lambda;
         if (t <= 0.)
            ftab[i] = LN_EPS;
         else if (t >= 1 - Num.DBL_EPSILON)
            ftab[i] = Math.log (1. / Num.DBL_EPSILON);
         else
            ftab[i] = Math.log (t/(1. - t));
         sum += ftab[i];
      }
      double empiricalMean = sum / n;

      sum = 0.0;
      for (int i = 0; i < n; i++) {
         t = ftab[i] - empiricalMean;
         sum += t * t;
      }
      double sigmaf = Math.sqrt(sum / n);

      double[] param = new double[2];
      param[0] = -empiricalMean / sigmaf;
      param[1] = 1.0 / sigmaf;

      return param;
   }

   /**
    * Creates a new instance of a `JohnsonSBDist` object using the maximum
    * likelihood method based on the @f$n@f$ observations @f$x[i]@f$, @f$i
    * = 0, 1, …, n-1@f$. Given the parameters @f$\xi= \mathtt{xi}@f$ and
    * @f$\lambda= \mathtt{lambda}@f$, the parameters @f$\gamma@f$ and
    * @f$\delta@f$ are estimated from the observations.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @param xi           parameter @f$\xi@f$
    *  @param lambda       parameter @f$\lambda@f$
    */
   public static JohnsonSBDist getInstanceFromMLE (double[] x, int n,
                                                   double xi, double lambda) {
      double parameters[] = getMLE (x, n, xi, lambda);
      return new JohnsonSBDist (parameters[0], parameters[1], xi, lambda);
   }

   /**
    * Returns the mean @cite tJOH49a&thinsp; of the Johnson @f$S_B@f$
    * distribution with parameters @f$\gamma@f$, @f$\delta@f$,
    * @f$\xi@f$ and @f$\lambda@f$.
    */
   public static double getMean (double gamma, double delta,
                                 double xi, double lambda) {
      if (lambda <= 0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (delta <= 0)
         throw new IllegalArgumentException ("delta <= 0");
      double[] tpsi = new double[1];
      double mu = getMeanPsi (gamma, delta, xi, lambda, tpsi);
      return xi + lambda * mu;
   }

   /**
    * Returns the variance @cite tFLY06a&thinsp; of the Johnson @f$S_B@f$
    * distribution with parameters @f$\gamma@f$, @f$\delta@f$,
    * @f$\xi@f$ and @f$\lambda@f$.
    *  @return the variance of the Johnson @f$S_B@f$ distribution.
    */
   public static double getVariance (double gamma, double delta,
                                     double xi, double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (delta <= 0.0)
         throw new IllegalArgumentException ("delta <= 0");

      final int NMAX = 10000;
      final double EPS = 1.0e-15;

      double a1 = 1.0/(2.0*delta*delta);
      double a2 = (1.0 - 2.0*delta*gamma)/(2.0*delta*delta);
      double a3 = (gamma - 1./delta)/delta;
      double v;
      int n = 0;
      double tem = 1;
      double sum = 0;
      while (Math.abs(tem) > EPS*Math.abs(sum) && n < NMAX) {
         ++n;
         v = Math.exp(-n*gamma/delta) - Math.exp(n*a3);
         tem = n*Math.exp(-n*n*a1) * v / (1 + Math.exp(-2*n*a1));
       //  tem = n*Math.exp(-n*n*a1) * Math.sinh(n*a2) / Math.cosh(n*a1);
         sum += tem;
      }
      if (n >= NMAX)
         System.err.println ("JohnsonSBDist:  possible lack of accuracy on variance");
      double A = -sum / (delta*delta);

      a1 = Math.PI * Math.PI * delta * delta;
      a2 = Math.PI * delta * gamma;
      int j;
      n = 0;
      tem = 1;
      sum = 0;
      while (Math.abs(tem) > EPS*Math.abs(sum) && n < NMAX) {
         ++n;
         j = 2*n - 1;
         tem = j*Math.exp(-j*j*a1/2.0) * Math.cos(j*a2) / Math.sinh(j*a1);
         sum += tem;
      }
      if (n >= NMAX)
         System.err.println ("JohnsonSBDist:  possible lack of accuracy on variance");
      double B = 2.0* a1 * sum;

      a1 = 2*Math.PI * Math.PI * delta * delta;
      a2 = 2*Math.PI * delta * gamma;
      n = 0;
      tem = 1;
      sum = 0;
      while (Math.abs(tem) > EPS*Math.abs(sum) && n < NMAX) {
         ++n;
         tem = n * Math.exp(-n*n*a1) * Math.sin(n*a2);
         sum += tem;
      }
      if (n >= NMAX)
         System.err.println ("JohnsonSBDist:  possible lack of accuracy on variance");
      double C = - 4.0 * Math.PI * delta * sum;

      double D = Math.sqrt(2*Math.PI) * Math.exp(0.5 * gamma* gamma);
      double[] tpsi = new double[1];
      double mu = getMeanPsi (gamma, delta, xi, lambda, tpsi);

      double tvar = mu *(1 - delta * gamma) - mu *mu +
                    delta / tpsi[0] * (A - B - mu * C * D);
      return lambda*lambda*tvar;

   }

  /**
   * Returns the standard deviation of the Johnson @f$S_B@f$ distribution
   * with parameters @f$\gamma@f$, @f$\delta@f$, @f$\xi@f$,
   * @f$\lambda@f$.
   *  @return the standard deviation of the Johnson @f$S_B@f$ distribution
   */
  public static double getStandardDeviation (double gamma, double delta,
                                             double xi, double lambda) {
      return Math.sqrt (JohnsonSBDist.getVariance (gamma, delta, xi, lambda));
   }

   /**
    * Sets the value of the parameters @f$\gamma@f$, @f$\delta@f$,
    * @f$\xi@f$ and @f$\lambda@f$ for this object.
    */
   public void setParams (double gamma, double delta,
                          double xi, double lambda) {
      setParams0(gamma, delta, xi, lambda);
      setLastParams(xi, lambda);
   }

}