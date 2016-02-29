/*
 * Class:        GeometricDist
 * Description:  geometric distribution
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
import  umontreal.ssj.util.Num;

/**
 * Extends the class  @ref DiscreteDistributionInt for the *geometric*
 * distribution @cite sLAW00a&thinsp; (page 322) with parameter @f$p@f$,
 * where @f$0 < p < 1@f$. Its mass function is
 * @anchor REF_probdist_GeometricDist_eq_fgeom
 * @f[
 *   p (x) = p  (1-p)^x, \qquad\mbox{for } x = 0, 1, 2, …\tag{fgeom}
 * @f]
 * The distribution function is given by
 * @anchor REF_probdist_GeometricDist_eq_Fgeom
 * @f[
 *   F (x) = 1 - (1-p)^{x+1}, \qquad\mbox{for } x = 0, 1, 2, …\tag{Fgeom}
 * @f]
 * and its inverse is
 * @anchor REF_probdist_GeometricDist_eq_FInvgeom
 * @f[
 *   F^{-1}(u) = \left\lfloor\frac{\ln(1 - u)}{\ln(1 - p)} \right\rfloor, \qquad\mbox{for } 0 \le u < 1. \tag{FInvgeom}
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_discrete
 */
public class GeometricDist extends DiscreteDistributionInt {

   private double p;
   private double vp;

   /**
    * Constructs a geometric distribution with parameter @f$p@f$.
    */
   public GeometricDist (double p) {
      setP(p);
   }
   public double prob (int x) {
      return prob (p, x);
   }

   public double cdf (int x) {
      return cdf (p, x);
   }

   public double barF (int x) {
      return barF (p, x);
   }

   public int inverseFInt (double u) {
        if (u > 1.0 || u < 0.0)
            throw new IllegalArgumentException ("u not in [0,1]");

        if (p >= 1.0)
            return 0;
        if (u <= p)
            return 0;
        if (u >= 1.0 || p <= 0.0)
            return Integer.MAX_VALUE;

        return (int)Math.floor (Math.log1p(-u)/vp);
   }

   public double getMean() {
      return GeometricDist.getMean (p);
   }

   public double getVariance() {
      return GeometricDist.getVariance (p);
   }

   public double getStandardDeviation() {
      return GeometricDist.getStandardDeviation (p);
   }

   /**
    * Computes the geometric probability @f$p(x)@f$ given in (
    * {@link REF_probdist_GeometricDist_eq_fgeom fgeom} ) .
    */
   public static double prob (double p, int x) {
      if (p < 0 || p > 1)
         throw new IllegalArgumentException ("p not in range (0,1)");
      if (p <= 0)
         return 0;
      if (p >= 1)
         return 0;
      if (x < 0)
         return 0;
      return p*Math.pow (1 - p, x);
   }

   /**
    * Computes the distribution function @f$F(x)@f$.
    */
   public static double cdf (double p, int x) {
      if (p < 0.0 || p > 1.0)
        throw new IllegalArgumentException ("p not in [0,1]");
      if (x < 0)
         return 0.0;
      if (p >= 1.0)                  // In fact, p == 1
         return 1.0;
      if (p <= 0.0)                  // In fact, p == 0
         return 0.0;
      return 1.0 - Math.pow (1.0 - p, (double)x + 1.0);
   }

   /**
    * Computes the complementary distribution function. *WARNING:* The
    * complementary distribution function is defined as @f$\bar{F}(x) =
    * P[X \ge x]@f$.
    */
   public static double barF (double p, int x) {
      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in [0,1]");
      if (x < 0)
         return 1.0;
      if (p >= 1.0)                  // In fact, p == 1
         return 0.0;
      if (p <= 0.0)                  // In fact, p == 0
         return 1.0;

      return Math.pow (1.0 - p, x);
   }

   /**
    * Computes the inverse of the geometric distribution, given by (
    * {@link REF_probdist_GeometricDist_eq_FInvgeom
    * FInvgeom} ).
    */
   public static int inverseF (double p, double u) {
        if (p > 1.0 || p < 0.0)
            throw new IllegalArgumentException ( "p not in [0,1]");
        if (u > 1.0 || u < 0.0)
            throw new IllegalArgumentException ("u not in [0,1]");
        if (p >= 1.0)
            return 0;
        if (u <= p)
           return 0;
        if (u >= 1.0 || p <= 0.0)
            return Integer.MAX_VALUE;

         double v = Math.log1p (-p);
         return (int)Math.floor (Math.log1p (-u)/v);
   }

   /**
    * Estimates the parameter @f$p@f$ of the geometric distribution using
    * the maximum likelihood method, from the @f$n@f$ observations
    * @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$. The estimate is returned in
    * element 0 of the returned array.  The maximum likelihood estimator
    * @f$\hat{p}@f$ satisfies the equation (see @cite sLAW00a&thinsp;
    * (page 323))
    * @f{align*}{
    *    \hat{p} = \frac{1}{\bar{x}_n + 1\Rule{0.0pt}{11.0pt}{0.0pt}}
    * @f}
    * where @f$\bar{x}_n@f$ is the average of @f$x[0], …, x[n-1]@f$.
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameter [@f$\hat{p}@f$]
    */
   public static double[] getMLE (int[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double parameters[];
      parameters = new double[1];
      double sum = 0.0;
      for (int i = 0; i < n; i++) {
         sum += x[i];
      }

      parameters[0] = 1.0 / (((double) sum / (double) n) + 1.0);

      return parameters;
   }

   /**
    * Creates a new instance of a geometric distribution with parameter
    * @f$p@f$ estimated using the maximum likelihood method based on the
    * @f$n@f$ observations @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static GeometricDist getInstanceFromMLE (int[] x, int n) {

      double parameters[] = getMLE (x, n);

      return new GeometricDist (parameters[0]);
   }

   /**
    * Computes and returns the mean @f$E[X] = (1 - p)/p@f$ of the
    * geometric distribution with parameter @f$p@f$.
    *  @return the mean of the geometric distribution @f$E[X] = (1 - p) /
    * p@f$
    */
   public static double getMean (double p) {
      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in range (0,1)");

      return (1 - p) / p;
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] = (1 - p)/p^2@f$
    * of the geometric distribution with parameter @f$p@f$.
    *  @return the variance of the Geometric distribution @f$\mbox{Var}[X]
    * = (1 - p) / p^2@f$
    */
   public static double getVariance (double p) {
      if (p < 0.0 || p > 1.0)
         throw new IllegalArgumentException ("p not in range (0,1)");

      return ((1 - p) / (p * p));
   }

   /**
    * Computes and returns the standard deviation of the geometric
    * distribution with parameter @f$p@f$.
    *  @return the standard deviation of the geometric distribution
    */
   public static double getStandardDeviation (double p) {
      return Math.sqrt (GeometricDist.getVariance (p));
   }

   /**
    * Returns the @f$p@f$ associated with this object.
    */
   public double getP() {
      return p;
   }

   /**
    * Resets the value of @f$p@f$ associated with this object.
    */
   public void setP (double p) {
      if (p < 0 || p > 1)
         throw new IllegalArgumentException ("p not in range (0,1)");
      vp = Math.log1p (-p);
      this.p = p;
      supportA = 0;
   }

   /**
    * Return a table containing the parameters of the current
    * distribution.
    */
   public double[] getParams () {
      double[] retour = {p};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : p = " + p;
   }

}