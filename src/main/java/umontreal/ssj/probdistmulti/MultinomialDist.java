/*
 * Class:        MultinomialDist
 * Description:  multinomial distribution
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
package umontreal.ssj.probdistmulti;
import umontreal.ssj.util.Num;

/**
 * Implements the abstract class  @ref DiscreteDistributionIntMulti for the
 * *multinomial* distribution with parameters @f$n@f$ and (@f$p_1@f$,
 * …,@f$p_d@f$). The probability mass function is @cite tJOH69a&thinsp;
 * @anchor REF_probdistmulti_MultinomialDist_eq_fMultinomial
 * @f[
 *   P[X = (x_1,…,x_d)] = {n!} \prod_{i=1}^d\frac{p_i^{x_i}}{x_i!}, \tag{fMultinomial}
 * @f]
 * where @f$\sum_{i=1}^d x_i = n@f$ and @f$\sum_{i=1}^d p_i = 1@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdistmulti_discrete
 */
public class MultinomialDist extends DiscreteDistributionIntMulti {
   protected int n;
   protected double p[];

   /**
    * Creates a `MultinomialDist` object with parameters @f$n@f$ and
    * (@f$p_1@f$,…,@f$p_d@f$) such that @f$\sum_{i=1}^d p_i = 1@f$. We
    * have @f$p_i = @f$ `p[i-1]`.
    */
   public MultinomialDist (int n, double p[]) {
      setParams (n, p);
   }

   public double prob (int x[]) {
      return prob_ (n, p, x);
   }

   public double cdf (int x[]) {
      return cdf_ (n, p, x);
   }

   public double[] getMean() {
      return getMean_ (n, p);
   }

   public double[][] getCovariance() {
      return getCovariance_ (n, p);
   }

   public double[][] getCorrelation () {
      return getCorrelation_ (n, p);
   }

   private static void verifParam(int n, double p[]) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double sumPi = 0.0;
      for (int i = 0; i < p.length; i++) {
         if ((p[i] < 0) || (p[i] > 1))
            throw new IllegalArgumentException("p is not a probability vector");
         sumPi += p[i];
      }

      if (Math.abs(sumPi - 1.0) > 1.0e-15)
         throw new IllegalArgumentException ("p is not a probability vector");
   }

   private static double prob_ (int n, double p[], int x[]) {
      if (x.length != p.length)
         throw new IllegalArgumentException ("x and p must have the same dimension");

      double sumXFact = 0.0;
      int sumX = 0;
      double sumPX = 0.0;

      for (int i = 0; i < p.length; i++) {
         sumX += x[i];
         sumXFact += Num.lnFactorial (x[i]);
         sumPX += (x[i] * Math.log (p[i]));
      }

      if (sumX != n)
         return 0.0;
      else {
         return Math.exp (Num.lnFactorial (n) - sumXFact + sumPX);
      }
   }

/**
 * Computes the probability mass function (
 * {@link REF_probdistmulti_MultinomialDist_eq_fMultinomial
 * fMultinomial} ) of the multinomial distribution with parameters @f$n@f$
 * and (@f$p_1@f$,…,@f$p_d@f$) evaluated at @f$x@f$.
 */
public static double prob (int n, double p[], int x[]) {
      verifParam (n, p);
      return prob_ (n, p, x);
   }


   private static double cdf_ (int n, double p[], int x[]) {
      boolean end = false;
      double sum = 0.0;
      int j;

      if (x.length != p.length)
         throw new IllegalArgumentException ("x and p must have the same dimension");

      int is[] = new int[x.length];
      for (int i = 0; i < is.length; i++)
         is[i] = 0;

      sum = 0.0;
      while (! end) {
         sum += prob (n, p, is);
         is[0]++;

         if (is[0] > x[0]) {
            is[0] = 0;
            j = 1;
            while (j < x.length && is[j] == x[j])
               is[j++] = 0;

            if (j == x.length)
               end = true;
            else
               is[j]++;
         }
      }

      return sum;
   }

/**
 * Computes the function @f$F@f$ of the multinomial distribution with
 * parameters @f$n@f$ and (@f$p_1@f$,…,@f$p_d@f$) evaluated at @f$x@f$.
 */
public static double cdf (int n, double p[], int x[]) {
      verifParam (n, p);
      return cdf_ (n, p, x);
   }


   private static double[] getMean_ (int n, double[] p) {
      double mean[] = new double[p.length];

      for (int i = 0; i < p.length; i++)
         mean[i] = n * p[i];

      return mean;
   }

/**
 * Computes the mean @f$E[X_i] = np_i@f$ of the multinomial distribution with
 * parameters @f$n@f$ and (@f$p_1@f$,…,@f$p_d@f$).
 */
public static double[] getMean (int n, double[] p) {
      verifParam (n, p);

      return getMean_ (n, p);
   }


   private static double[][] getCovariance_ (int n, double[] p) {
      double cov[][] = new double[p.length][p.length];

      for (int i = 0; i < p.length; i++) {
         for (int j = 0; j < p.length; j++)
            cov[i][j] = -n * p[i] * p[j];

         cov[i][i] = n * p[i] * (1.0 - p[i]);
      }
      return cov;
   }

/**
 * Computes the covariance matrix of the multinomial distribution with
 * parameters @f$n@f$ and (@f$p_1@f$,…,@f$p_d@f$).
 */
public static double[][] getCovariance (int n, double[] p) {
      verifParam (n, p);
      return getCovariance_ (n, p);
   }


   private static double[][] getCorrelation_ (int n, double[] p) {
      double corr[][] = new double[p.length][p.length];

      for (int i = 0; i < p.length; i++) {
         for (int j = 0; j < p.length; j++)
            corr[i][j] = -Math.sqrt(p[i] * p[j] / ((1.0 - p[i]) * (1.0 - p[j])));
         corr[i][i] = 1.0;
      }
      return corr;
   }

/**
 * Computes the correlation matrix of the multinomial distribution with
 * parameters @f$n@f$ and (@f$p_1@f$,…,@f$p_d@f$).
 */
public static double[][] getCorrelation (int n, double[] p) {
      verifParam (n, p);
      return getCorrelation_ (n, p);
   }

   /**
    * Estimates and returns the parameters
    * [@f$\hat{p_i}@f$,…,@f$\hat{p_d}@f$] of the multinomial distribution
    * using the maximum likelihood method. It uses the @f$m@f$
    * observations of @f$d@f$ components in table @f$x[i][j]@f$, @f$i = 0,
    * 1, …, m-1@f$ and @f$j = 0, 1, …, d-1@f$. The equations of the
    * maximum likelihood are defined as
    * @f{align*}{
    *    \hat{p}_i = \frac{\bar{X_i}}{N}.
    * @f}
    * @param x            the list of observations used to evaluate
    *                      parameters
    *  @param m            the number of observations used to evaluate
    *                      parameters
    *  @param d            the dimension of each observation
    *  @param n            the number of independant trials for each
    *                      series
    *  @return returns the parameters [@f$\hat{p_i}@f$,…,@f$\hat{p_d}@f$]
    */
   public static double[] getMLE (int x[][], int m, int d, int n) {
      double parameters[] = new double[d];
      double xBar[] = new double[d];
      double N = 0.0;

      if (m <= 0)
         throw new IllegalArgumentException ("m <= 0");
      if (d <= 0)
         throw new IllegalArgumentException ("d <= 0");

      for (int i = 0; i < d; i++)
         xBar[i] = 0;

      for (int v = 0; v < m; v++)
         for (int c = 0; c < d; c++)
            xBar[c] += x[v][c];

      for (int i = 0; i < d; i++)
      {
         xBar[i] = xBar[i] / (double) n;
         N += xBar[i];
      }
      if (N != (double) n)
         throw new IllegalArgumentException("n is not correct");

      for (int i = 0; i < d; i++)
         parameters[i] = xBar[i] / (double) n;

      return parameters;
   }

   /**
    * Returns the parameter @f$n@f$ of this object.
    */
   public int getN() {
      return n;
   }

   /**
    * Returns the parameters (@f$p_1@f$,…,@f$p_d@f$) of this object.
    */
   public double[] getP() {
      return p;
   }

   /**
    * Sets the parameters @f$n@f$ and (@f$p_1@f$,…,@f$p_d@f$) of this
    * object.
    */
   public void setParams (int n, double p[]) {
      double sumP = 0.0;

      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      if (p.length < 2)
         throw new IllegalArgumentException ("p.length < 2");

      this.n = n;
      this.dimension = p.length;
      this.p = new double[dimension];
      for (int i = 0; i < dimension; i++)
      {
         if ((p[i] < 0) || (p[i] > 1))
            throw new IllegalArgumentException("p is not a probability vector");

         this.p[i] = p[i];
         sumP += p[i];
      }

      if (Math.abs(sumP - 1.0) > 1.0e-15)
         throw new IllegalArgumentException ("p is not a probability vector");
   }

}
