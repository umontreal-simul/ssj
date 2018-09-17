/*
 * Class:        ChiSquareDist
 * Description:  chi-square distribution
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
import umontreal.ssj.functions.MathFunction;

/**
 * Extends the class  @ref ContinuousDistribution for the *chi-square*
 * distribution with @f$n@f$ degrees of freedom, where @f$n@f$ is a positive
 * integer @cite tJOH95a&thinsp; (page 416). Its density is
 * @anchor REF_probdist_ChiSquareDist_eq_Fchi2
 * @f[
 *   f(x) = \frac{x^{(n/2)-1}e^{-x/2}}{2^{n/2}\Gamma(n/2)},\qquad\mbox{for } x > 0 \tag{Fchi2}
 * @f]
 * where @f$\Gamma(x)@f$ is the gamma function defined in (
 * {@link REF_probdist_GammaDist_eq_Gamma Gamma} ). The
 * *chi-square* distribution is a special case of the *gamma* distribution
 * with shape parameter @f$n/2@f$ and scale parameter @f$1/2@f$. Therefore,
 * one can use the methods of  @ref GammaDist for this distribution.
 *
 * The non-static versions of the methods `cdf`, `barF`, and `inverseF` call
 * the static version of the same name.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class ChiSquareDist extends ContinuousDistribution {
   protected int n;
   protected double C1;

   private static class Function implements MathFunction {
      protected int n;
      protected double sumLog;

      public Function (double s, int n)
      {
         this.n = n;
         this.sumLog = s;
      }

      public double evaluate (double k)
      {
         if (k < 1.0) return 1.0e200;
         return (sumLog + n * (Num.lnGamma (k / 2.0) - 0.5*Num.LN2 - Num.lnGamma ((k + 1.0) / 2.0)));
      }
   }

   /**
    * Constructs a chi-square distribution with `n` degrees of freedom.
    */
   public ChiSquareDist (int n) {
      setN (n);
   }


   public double density (double x) {
      if (x <= 0)
         return 0.0;
      return Math.exp ((n/2.0 - 1)*Math.log (x) - x/2.0 - C1);
   }

   public double cdf (double x) {
      return cdf (n, decPrec, x);
   }

   public double barF (double x) {
      return barF (n, decPrec, x);
   }

   public double inverseF (double u) {
      return inverseF (n, u);
   }

   public double getMean() {
      return ChiSquareDist.getMean (n);
   }

   public double getVariance() {
      return ChiSquareDist.getVariance (n);
   }

   public double getStandardDeviation() {
      return ChiSquareDist.getStandardDeviation (n);
   }

/**
 * Computes the density function (
 * {@link REF_probdist_ChiSquareDist_eq_Fchi2 Fchi2} ) for a
 * *chi-square* distribution with @f$n@f$ degrees of freedom.
 */
public static double density (int n, double x) {
      if (x <= 0)
         return 0.0;
      return Math.exp ((n/2.0 - 1)*Math.log (x) - x/2.0 -
                  (n/2.0)*Num.LN2 - Num.lnGamma(n/2.0));
   }

   /**
    * Computes the chi-square distribution function with @f$n@f$ degrees
    * of freedom, evaluated at @f$x@f$. The method tries to return @f$d@f$
    * decimals digits of precision, but there is no guarantee.
    */
   public static double cdf (int n, int d, double x) {
      if (n <= 0)
        throw new IllegalArgumentException ("n <= 0");
      if (x <= 0.0)
         return 0.0;
      if (x >= XBIG*n)
         return 1.0;
      return GammaDist.cdf (n/2.0, d, x/2.0);
   }

   /**
    * Computes the complementary chi-square distribution function with
    * @f$n@f$ degrees of freedom, evaluated at @f$x@f$. The method tries
    * to return @f$d@f$ decimals digits of precision, but there is no
    * guarantee.
    */
   public static double barF (int n, int d, double x) {
      if (n <= 0)
        throw new IllegalArgumentException ("n <= 0");
      if (x <= 0.0)
         return 1.0;
      return GammaDist.barF (n/2.0, d, x/2.0);
   }

   /**
    * Computes an approximation of @f$F^{-1}(u)@f$, where @f$F@f$ is the
    * chi-square distribution with @f$n@f$ degrees of freedom. Uses the
    * approximation given in @cite tBES75a&thinsp; and in Figure L.23 of
    * @cite sBRA87a&thinsp;. It gives at least 6 decimal digits of
    * precision, except far in the tails (that is, for @f$u< 10^{-5}@f$ or
    * @f$u > 1 - 10^{-5}@f$) where the function calls the method
    * `GammaDist.inverseF (n/2, 7, u)` and multiplies the result by 2.0.
    * To get better precision, one may call `GammaDist.inverseF`, but this
    * method is slower than the current method, especially for large
    * @f$n@f$. For instance, for @f$n = @f$ 16, 1024, and 65536, the
    * `GammaDist.inverseF` method is 2, 5, and 8 times slower,
    * respectively, than the current method.
    */
   public static double inverseF (int n, double u) {
      /*
       * Returns an approximation of the inverse of Chi square cdf
       * with n degrees of freedom.
       * As in Figure L.23 of P.Bratley, B.L.Fox, and L.E.Schrage.
       *    A Guide to Simulation Springer-Verlag,
       *    New York, second edition, 1987.
       */
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u must be in [0,1]");
      if (u == 1.0)
         return Double.POSITIVE_INFINITY;
      if (u == 0.0)
         return 0.0;

      final double E = 0.5e-5;    // Precision of this approximation
      final double AA = 0.6931471805;
      double A, XX, X, C, G, CH, Q, P1, P2, T, B, S1, S2, S3, S4, S5, S6;

      if (u < 0.00001 || u > 1.0 - 1.0e-5)
         return 2.0 * GammaDist.inverseF (n / 2.0, 7, u);
      if (u >= 1.0)
         return n * XBIG;
      if (u >= 0.999998)
         return (n + 4.0 * Math.sqrt (2.0 * n));

      G = Num.lnGamma (n / 2.0);
      XX = 0.5 * n;
      C = XX - 1.0;

      if (n >= -1.24 * Math.log (u)) {
         X = NormalDist.inverseF01 (u);
         P1 = 0.222222 / n;
         Q = X * Math.sqrt (P1) + 1.0 - P1;
         CH = n * Q * Q * Q;
         if (CH > 2.2 * n + 6.0)
            CH = -2.0 * (Math.log1p (-u) - C * Math.log (0.5 * CH) + G);

      } else {
         CH = Math.pow (u * XX * Math.exp (G + XX * AA), 1.0 / XX);
         if (CH - E < 0)
            return CH;
      }

      Q = CH;
      P1 = 0.5 * CH;
      P2 = u - GammaDist.cdf (XX, 5, P1);
      if (GammaDist.cdf (XX, 5, P1) == -1.0)
         throw new IllegalArgumentException ("RESULT = -1");

      T = P2 * Math.exp (XX * AA + G + P1 - C * Math.log (CH));
      B = T / CH;
      A = 0.5 * T - B * C;
      S1 = (210.0 + A * (140.0 +
            A * (105.0 + A * (84.0 + A * (70.0 + 60.0 * A))))) / 420.0;
      S2 = (420.0 + A * (735.0 + A * (966.0 + A * (1141.0 + 1278.0 * A))))
         / 2520.0;
      S3 = (210.0 + A * (462.0 + A * (707.0 + 932.0 * A))) / 2520.0;
      S4 = (252.0 + A * (672.0 + 1182.0 * A) +
         C * (294.0 + A * (889.0 + 1740.0 * A))) / 5040.0;
      S5 = (84.0 + 264.0 * A + C * (175.0 + 606.0 * A)) / 2520.0;
      S6 = (120.0 + C * (346.0 + 127.0 * C)) / 5040.0;
      CH = CH + T * (1.0 + 0.5 * T * S1 - B * C * (S1 - B * (S2 -
               B * (S3 - B * (S4 - B * (S5 - B * S6))))));

      double temp;
      while (Math.abs (Q / CH - 1.0) > E) {
         Q = CH;
         P1 = 0.5 * CH;
         temp = GammaDist.cdf (XX, 6, P1);
         P2 = u - temp;

         if (temp == -1.0)
            return -1.0;

         T = P2 * Math.exp (XX * AA + G + P1 - C * Math.log (CH));
         B = T / CH;
         A = 0.5 * T - B * C;
         S1 = (210.0 + A * (140.0 + A * (105.0 + A * (84.0 +
                     A * (70.0 + 60.0 * A))))) / 420.0;
         S2 = (420.0 + A * (735.0 + A * (966.0 + A * (1141.0 +
                     1278.0 * A)))) / 2520.0;
         S3 = (210.0 + A * (462.0 + A * (707.0 + 932.0 * A))) / 2520.0;
         S4 = (252.0 + A * (672.0 + 1182.0 * A) +
            C * (294.0 + A * (889.0 + 1740.0 * A))) / 5040.0;
         S5 = (84.0 + 264.0 * A + C * (175.0 + 606.0 * A)) / 2520.0;
         S6 = (120.0 + C * (346.0 + 127.0 * C)) / 5040.0;
         CH = CH + T * (1.0 + 0.5 * T * S1 - B * C * (S1 - B * (S2 -
                  B * (S3 - B * (S4 - B * (S5 - B * S6))))));
      }
      return CH;
   }

   /**
    * Estimates the parameter @f$n@f$ of the chi-square distribution using
    * the maximum likelihood method, from the @f$m@f$ observations
    * @f$x[i]@f$, @f$i = 0, 1, …, m-1@f$. The estimate is returned in
    * element 0 of the returned array.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param m            the number of observations to use to evaluate
    *                      parameters
    *  @return returns the parameter [@f$\hat{n}@f$]
    */
   public static double[] getMLE (double[] x, int m) {
      if (m <= 0)
         throw new IllegalArgumentException ("m <= 0");

      double[] parameters;

      parameters = getMomentsEstimate (x, m);
      double k = Math.round (parameters[0]) - 5.0;
      if (k < 1.0)
         k = 1.0;

      double sum = 0.0;
      for (int i = 0; i < m; i++) {
         if (x[i] > 0.0)
            sum += 0.5*Math.log (x[i]);
         else
            sum -= 709.0;
      }

      Function f = new Function (sum, m);
      while (f.evaluate(k) > 0.0)
         k++;
      parameters[0] = k;

      return parameters;
   }

   /**
    * Creates a new instance of a chi-square distribution with parameter
    * @f$n@f$ estimated using the maximum likelihood method based on the
    * @f$m@f$ observations @f$x[i]@f$, @f$i = 0, 1, …, m-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param m            the number of observations to use to evaluate
    *                      parameters
    */
   public static ChiSquareDist getInstanceFromMLE (double[] x, int m) {
      double parameters[] = getMLE (x, m);
      return new ChiSquareDist ((int) parameters[0]);
   }

   /**
    * Computes and returns the mean @f$E[X] = n@f$ of the chi-square
    * distribution with parameter @f$n@f$.
    *  @return the mean of the Chi-square distribution @f$E[X] = n@f$
    */
   public static double getMean (int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("degrees of freedom " +
                              "must be non-null and positive.");

      return n;
   }

   /**
    * Estimates and returns the parameter [@f$\hat{n}@f$] of the
    * chi-square distribution using the moments method based on the
    * @f$m@f$ observations in table @f$x[i]@f$, @f$i = 0, 1, …, m-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param m            the number of observations to use to evaluate
    *                      parameters
    *  @return returns the parameter [@f$\hat{n}@f$]
    */
   public static double[] getMomentsEstimate (double[] x, int m) {
      double[] parameters = new double[1];

      double sum = 0.0;
      for (int i = 0; i < m; i++)
         sum += x[i];
      parameters[0] = sum / (double) m;

      return parameters;
   }

   /**
    * Returns the variance @f$\mbox{Var}[X] = 2n@f$ of the chi-square
    * distribution with parameter @f$n@f$.
    *  @return the variance of the chi-square distribution @f$\mbox{Var}X]
    * = 2n@f$
    */
   public static double getVariance (int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("degrees of freedom " +
                              "must be non-null and positive.");

      return (2 * n);
   }

   /**
    * Returns the standard deviation of the chi-square distribution with
    * parameter @f$n@f$.
    *  @return the standard deviation of the chi-square distribution
    */
   public static double getStandardDeviation (int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("degrees of freedom " +
                              "must be non-null and positive.");

      return Math.sqrt(2 * n);
   }

   /**
    * Returns the parameter @f$n@f$ of this object.
    */
   public int getN() {
      return n;
   }

   /**
    * Sets the parameter @f$n@f$ of this object.
    */
   public void setN (int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("degrees of freedom " +
                              "must be non-null and positive.");

      this.n = n;
      supportA = 0.0;
      C1 = 0.5 * n *Num.LN2 + Num.lnGamma(n/2.0);
   }

   /**
    * Return a table containing the parameters of the current
    * distribution.
    */
   public double[] getParams () {
      double[] retour = {n};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : n = " + n;
   }

}