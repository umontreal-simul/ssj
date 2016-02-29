/*
 * Class:        JohnsonSLDist
 * Description:  Johnson S_L distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        July 2012
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
import umontreal.ssj.functions.MathFunction;
import optimization.*;

/**
 * Extends the class  @ref ContinuousDistribution for the <em>Johnson
 * @f$S_L@f$</em> distribution (see @cite tJOH49a, @cite tJOH95a&thinsp;). It
 * has shape parameters @f$\gamma@f$ and @f$\delta> 0@f$, location
 * parameter @f$\xi@f$, and scale parameter @f$\lambda> 0@f$. Denoting
 * @f$t=(x-\xi)/\lambda@f$ and @f$z = \gamma+ \delta\ln(t)@f$, the
 * distribution has density
 * @f[
 *   f(x) = \frac{\delta e^{-z^2/2}}{\lambda t \sqrt{2\pi}},  \qquad\mbox{for } \xi< x < \infty,
 * @f]
 * and distribution function
 * @f[
 *   F(x) = \Phi(z), \qquad\mbox{for } \xi< x < \infty,
 * @f]
 * where @f$\Phi@f$ is the standard normal distribution function. The
 * inverse distribution function is
 * @f[
 *   F^{-1} (u) = \xi+ \lambda e^{v(u)}, \qquad\mbox{for } 0 \le u \le1,
 * @f]
 * where
 * @f[
 *   v(u) = [\Phi^{-1}(u) - \gamma]/\delta.
 * @f]
 * Without loss of generality, one may choose @f$\gamma= 0@f$ or
 * @f$\lambda=1@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class JohnsonSLDist extends JohnsonSystem {

   private static class Function implements MathFunction
   {
      // To find value of t > 0 in (Johnson 1949, eq. 16)
      protected double a;

      public Function (double sb1) {
         a = sb1;
      }

      public double evaluate (double t) {
         return (t*t*t - 3*t - a);
      }
   }


   private static double[] initPar (double[] x, int n, double xmin)
   {
      // Use moments to estimate initial values of params as input to MLE
      // (Johnson 1949, Biometrika 36, p. 149)
      int j;
      double sum = 0.0;
      for (j = 0; j < n; j++) {
         sum += x[j];
      }
      double mean = sum / n;

      double v;
      double sum3 = 0.0;
      sum = 0;
      for (j = 0; j < n; j++) {
         v = x[j] - mean;
         sum += v*v;
         sum3 += v*v*v;
      }
      double m2 = sum / n;
      double m3 = sum3 / n;

      v = m3 / Math.pow (m2, 1.5);
      Function f = new Function (v);
      double t0 = 0;
      double tlim = Math.cbrt (v);
      if (tlim <= 0) {
         t0 = tlim;
         tlim = 10;
      }
      double t = RootFinder.brentDekker (t0, tlim, f, 1e-5);
      if (t <= 0)
         throw new UnsupportedOperationException("t <= 0;   no MLE");
      double xi = mean - Math.sqrt(m2 / t);
      if (xi >= xmin)
         // throw new UnsupportedOperationException("xi >= xmin;   no MLE");
         xi = xmin - 1.0e-1;
      v = 1 + m2 / ((mean - xi)*(mean - xi));
      double delta = 1.0 / Math.sqrt((Math.log(v)));
      v = Math.sqrt(v);
      double lambda = (mean - xi) / v;
      double[] param = new double[3];
      param[0] = delta;
      param[1] = xi;
      param[2] = lambda;
      return param;
   }


   private static class Optim implements Uncmin_methods
   {
      // minimizes the loglikelihood function
      private int n;
      private double[] X;
      private double xmin;      // min{X_j}
      private static final double BARRIER = 1.0e100;

      public Optim (double[] X, int n, double xmin) {
         this.n = n;
         this.X = X;
         this.xmin = xmin;
      }

      public double f_to_minimize (double[] par) {
         // par = [0, delta, xi, lambda]
         // arrays in Uncmin starts at index 1; par[0] is unused
         double delta  = par[1];
         double xi     = par[2];
         double lambda = par[3];
         if (delta <= 0.0 || lambda <= 0.0)         // barrier at 0
            return BARRIER;
         if (xi >= xmin)
            return BARRIER;

         double loglam = Math.log(lambda);
         double v, z;
         double sumv = 0.0;
         double sumz = 0.0;
         for (int j = 0; j < n; j++) {
            v = Math.log (X[j] - xi);
            z = delta * (v - loglam);
            sumv += v;
            sumz += z*z;
         }

         return sumv + sumz / 2.0 - n *Math.log(delta);
      }

      public void gradient (double[] x, double[] g)
      {
      }

      public void hessian (double[] x, double[][] h)
      {
      }
   }

   /**
    * Same as  {@link #JohnsonSLDist(double,double,double,double)
    * JohnsonSLDist(gamma, delta, 0, 1)}.
    */
   public JohnsonSLDist (double gamma, double delta) {
      this (gamma, delta, 0, 1);
   }

   /**
    * Constructs a `JohnsonSLDist` object with shape parameters
    * @f$\gamma@f$ and @f$\delta@f$, location parameter @f$\xi@f$, and
    * scale parameter @f$\lambda@f$.
    */
   public JohnsonSLDist (double gamma, double delta,
                         double xi, double lambda) {
      super (gamma, delta, xi, lambda);
      setLastParams(xi);
   }
   private void setLastParams(double xi) {
      supportA = xi;
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
      return JohnsonSLDist.getMean (gamma, delta, xi, lambda);
   }

   public double getVariance() {
      return JohnsonSLDist.getVariance (gamma, delta, xi, lambda);
   }

   public double getStandardDeviation() {
      return JohnsonSLDist.getStandardDeviation (gamma, delta, xi, lambda);
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

      if (x <= xi)
         return 0;
      double y = (x - xi)/lambda;
      double z = gamma + delta*Math.log (y);
      if (z >= 1.0e10)
         return 0;
      return delta * Math.exp (-z*z/2.0)/(lambda*y*Math.sqrt (2.0*Math.PI));
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

      if (x <= xi)
         return 0.0;
      double y = (x - xi)/lambda;
      double z = gamma + delta*Math.log (y);
      return NormalDist.cdf01 (z);
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

      if (x <= xi)
         return 1.0;
      double y = (x - xi)/lambda;
      double z = gamma + delta*Math.log (y);
      return NormalDist.barF01 (z);
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
         return xi;

      double z = NormalDist.inverseF01 (u);
      double t = (z - gamma)/delta;
      if (t >= Num.DBL_MAX_EXP*Num.LN2)
         return Double.POSITIVE_INFINITY;
      if (t <= Num.DBL_MIN_EXP*Num.LN2)
         return xi;

      return xi + lambda * Math.exp(t);
   }

   /**
    * Estimates the parameters @f$(\gamma@f$, @f$\delta@f$, @f$\xi@f$,
    * @f$\lambda)@f$ of the <em>Johnson @f$S_L@f$</em> distribution using
    * the maximum likelihood method, from the @f$n@f$ observations
    * @f$x[i]@f$, @f$i = 0, 1,…, n-1@f$. The estimates are returned in a
    * 4-element array in the order [0, @f$\delta@f$, @f$\xi@f$,
    * @f$\lambda@f$] (with @f$\gamma@f$ always set to 0).
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @return returns the parameters [0, @f$\delta@f$, @f$\xi@f$,
    * @f$\lambda@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      int j;
      double xmin = Double.MAX_VALUE;
      for (j = 0; j < n; j++) {
         if (x[j] < xmin)
            xmin = x[j];
      }
      double[] paramIn = new double[3];
      paramIn = initPar(x, n, xmin);
      double[] paramOpt = new double[4];
      for (int i = 0; i < 3; i++)
         paramOpt[i+1] = paramIn[i];

      Optim system = new Optim (x, n, xmin);

      double[] xpls = new double[4];
      double[] fpls = new double[4];
      double[] gpls = new double[4];
      int[] itrcmd = new int[2];
      double[][] a = new double[4][4];
      double[] udiag = new double[4];

      Uncmin_f77.optif0_f77 (3, paramOpt, system, xpls, fpls, gpls,
                             itrcmd, a, udiag);

      double[] param = new double[4];
      param[0] = 0;
      for (int i = 1; i <= 3; i++)
         param[i] = xpls[i];
      return param;
   }

   /**
    * Creates a new instance of a <em>Johnson @f$S_L@f$</em> distribution
    * with parameters 0, @f$\delta@f$, @f$\xi@f$ and @f$\lambda@f$ over
    * the interval @f$[\xi,\infty]@f$ estimated using the maximum
    * likelihood method based on the @f$n@f$ observations @f$x[i]@f$, @f$i
    * = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static JohnsonSLDist getInstanceFromMLE (double[] x, int n) {
      double param[] = getMLE (x, n);
      return new JohnsonSLDist (0, param[0], param[1], param[2]);
   }

   /**
    * Returns the mean
    * @f[
    *   E[X] = \xi+ \lambda e^{1/2\delta^2 - \gamma/\delta}
    * @f]
    * of the Johnson @f$S_L@f$ distribution with parameters @f$\gamma@f$,
    * @f$\delta@f$, @f$\xi@f$ and @f$\lambda@f$.
    *  @return the mean of the Johnson @f$S_L@f$ distribution @f$E[X] =
    * \xi+ \lambda e^{1/2\delta^2 - \gamma/\delta}@f$
    */
   public static double getMean (double gamma, double delta,
                                 double xi, double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (delta <= 0.0)
         throw new IllegalArgumentException ("delta <= 0");

      double t = 1.0 / (2.0 * delta * delta) - gamma / delta;
      return (xi + lambda * Math.exp(t));
   }

   /**
    * Returns the variance
    * @f[
    *   \mbox{Var}[X] = \lambda^2 \left(e^{1/\delta^2} - 1\right) e^{1/\delta^2 - 2\gamma/\delta}
    * @f]
    * of the Johnson @f$S_L@f$ distribution with parameters @f$\gamma@f$,
    * @f$\delta@f$, @f$\xi@f$ and @f$\lambda@f$.
    *  @return the variance of the Johnson @f$S_L@f$ distribution
    * @f$\mbox{Var}[X] = \lambda^2 \left(e^{1/\delta^2} - 1\right)
    * e^{1/\delta^2 - 2\gamma/\delta}@f$
    */
   public static double getVariance (double gamma, double delta,
                                     double xi, double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (delta <= 0.0)
         throw new IllegalArgumentException ("delta <= 0");

      double t = 1.0 / (delta * delta) - 2 * gamma / delta;
      return lambda * lambda * Math.exp(t) * (Math.exp(1.0 / (delta * delta)) - 1);
   }

   /**
    * Returns the standard deviation of the Johnson @f$S_L@f$ distribution
    * with parameters @f$\gamma@f$, @f$\delta@f$, @f$\xi@f$,
    * @f$\lambda@f$.
    *  @return the standard deviation of the Johnson @f$S_L@f$
    * distribution
    */
   public static double getStandardDeviation (double gamma, double delta,
                                              double xi, double lambda) {
      return Math.sqrt (JohnsonSLDist.getVariance (gamma, delta, xi, lambda));
   }

   /**
    * Sets the value of the parameters @f$\gamma@f$, @f$\delta@f$,
    * @f$\xi@f$ and @f$\lambda@f$ for this object.
    */
   public void setParams (double gamma, double delta,
                          double xi, double lambda) {
      setParams0(gamma, delta, xi, lambda);
      setLastParams(xi);
   }

}