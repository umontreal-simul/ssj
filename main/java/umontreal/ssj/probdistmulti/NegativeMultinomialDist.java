/*
 * Class:        NegativeMultinomialDist
 * Description:  negative multinomial distribution
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
import umontreal.ssj.util.RootFinder;
import umontreal.ssj.functions.MathFunction;

/**
 * Implements the class  @ref DiscreteDistributionIntMulti for the *negative
 * multinomial* distribution with parameters @f$n > 0@f$ and (@f$p_1, …,
 * p_d@f$) such that all @f$0<p_i<1@f$ and @f$\sum_{i=1}^d p_i < 1@f$. The
 * probability mass function is @cite tJOH69a&thinsp;
 * @anchor REF_probdistmulti_NegativeMultinomialDist_eq_fNegativeMultinomial
 * @f[
 *   P[X = (x_1, …,x_d)] = \frac{\Gamma\left(n + \sum_{i=1}^d x_i\right) p_0^n}{\Gamma(n)} \prod_{i=1}^d \frac{p_i^{x_i}}{x_i!} \tag{fNegativeMultinomial}
 * @f]
 * where @f$p_0 = 1 - \sum_{i=1}^d p_i@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdistmulti_discrete
 */
public class NegativeMultinomialDist extends DiscreteDistributionIntMulti {
   protected double n;
   protected double p[];

   private static class Function implements MathFunction {
      protected double Fl[];
      protected int ups[];
      protected int k;
      protected int M;
      protected int sumUps;

      public Function (int k, int m, int ups[], double Fl[]) {
         this.k = k;
         this.M = m;

         this.Fl = new double[Fl.length];
         System.arraycopy (Fl, 0, this.Fl, 0, Fl.length);
         this.ups = new int[ups.length];
         System.arraycopy (ups, 0, this.ups, 0, ups.length);

         sumUps = 0;
         for (int i = 0; i < ups.length; i++)
            sumUps += ups[i];
      }

      public double evaluate (double gamma) {
         double sum = 0.0;
         for (int l = 0; l < M; l++)
            sum += (Fl[l] / (gamma + (double) l));
         return (sum - Math.log1p (sumUps / (k * gamma)));
      }
   }


   private static class FuncInv extends Function implements MathFunction {

      public FuncInv (int k, int m, int ups[], double Fl[]) {
         super (k, m, ups, Fl);
      }

      public double evaluate (double nu) {
         double sum = 0.0;
         for (int l = 0; l < M; l++)
            sum += Fl[l] / (1.0 + nu * l);
         return (sum *nu - Math.log1p (sumUps * nu / k));
      }
   }

   /**
    * Creates a `NegativeMultinomialDist` object with parameters @f$n@f$
    * and (@f$p_1@f$, …, @f$p_d@f$) such that @f$\sum_{i=1}^d p_i < 1@f$,
    * as described above. We have @f$p_i = @f$ `p[i-1]`.
    */
   public NegativeMultinomialDist (double n, double p[]) {
      setParams (n, p);
   }


   public double prob (int x[]) {
      return prob_ (n, p, x);
   }
/*
   public double cdf (int x[]) {
      throw new UnsupportedOperationException ("cdf not implemented");
   }
*/
   public double[] getMean() {
      return getMean_ (n, p);
   }

   public double[][] getCovariance() {
      return getCovariance_ (n, p);
   }

   public double[][] getCorrelation() {
      return getCorrelation_ (n, p);
   }

   private static void verifParam (double n, double p[]) {
      double sumPi = 0.0;

      if (n <= 0.0)
         throw new IllegalArgumentException("n <= 0");

      for (int i = 0; i < p.length;i++) {
         if ((p[i] < 0) || (p[i] >= 1))
            throw new IllegalArgumentException("p is not a probability vector");

         sumPi += p[i];
      }
      if (sumPi >= 1.0)
         throw new IllegalArgumentException("p is not a probability vector");
   }

   private static double prob_ (double n, double p[], int x[]) {
      double p0 = 0.0;
      double sumPi= 0.0;
      double sumXi= 0.0;
      double sumLnXiFact = 0.0;
      double sumXiLnPi = 0.0;

      if (x.length != p.length)
         throw new IllegalArgumentException ("x and p must have the same size");

      for (int i = 0; i < p.length;i++)
      {
         sumPi += p[i];
         sumXi += x[i];
         sumLnXiFact += Num.lnFactorial (x[i]);
         sumXiLnPi += x[i] * Math.log (p[i]);
      }
      p0 = 1.0 - sumPi;

      return Math.exp (Num.lnGamma (n + sumXi) - (Num.lnGamma (n) +
           sumLnXiFact) + n * Math.log (p0) + sumXiLnPi);
   }

/**
 * Computes the probability mass function (
 * {@link REF_probdistmulti_NegativeMultinomialDist_eq_fNegativeMultinomial
 * fNegativeMultinomial} ) of the negative multinomial distribution with
 * parameters @f$n@f$ and (@f$p_1@f$, …, @f$p_d@f$), evaluated at
 * @f$\mathbf{x}@f$.
 */
public static double prob (double n, double p[], int x[]) {
      verifParam (n, p);
      return prob_ (n, p, x);
   }


   private static double cdf_ (double n, double p[], int x[]) {
      throw new UnsupportedOperationException ("cdf not implemented");
   }

/**
 * Computes the cumulative probability function @f$F@f$ of the negative
 * multinomial distribution with parameters @f$n@f$ and (@f$p_1@f$, …,
 * @f$p_k@f$), evaluated at @f$\mathbf{x}@f$.
 */
public static double cdf (double n, double p[], int x[]) {
      verifParam (n, p);

      return cdf_ (n, p, x);
   }


   private static double[] getMean_ (double n, double p[]) {
      double p0 = 0.0;
      double sumPi= 0.0;
      double mean[] = new double[p.length];

      for (int i = 0; i < p.length;i++)
         sumPi += p[i];
      p0 = 1.0 - sumPi;

      for (int i = 0; i < p.length; i++)
         mean[i] = n * p[i] / p0;

      return mean;
   }

/**
 * Computes the mean @f$E[X] = n p_i / p_0@f$ of the negative multinomial
 * distribution with parameters @f$n@f$ and (@f$p_1@f$, …, @f$p_d@f$).
 */
public static double[] getMean (double n, double p[]) {
      verifParam (n, p);

      return getMean_ (n, p);
   }


   private static double[][] getCovariance_ (double n, double p[]) {
      double p0 = 0.0;
      double sumPi= 0.0;
      double cov[][] = new double[p.length][p.length];

      for (int i = 0; i < p.length;i++)
         sumPi += p[i];
      p0 = 1.0 - sumPi;

      for (int i = 0; i < p.length; i++)
      {
         for (int j = 0; j < p.length; j++)
            cov[i][j] = n * p[i] * p[j] / (p0 * p0);

         cov[i][i] = n * p[i] * (p[i] + p0) / (p0 * p0);
      }

      return cov;
   }

/**
 * Computes the covariance matrix of the negative multinomial distribution
 * with parameters @f$n@f$ and (@f$p_1@f$, …, @f$p_d@f$).
 */
public static double[][] getCovariance (double n, double p[]) {
      verifParam (n, p);

      return getCovariance_ (n, p);
   }


   private static double[][] getCorrelation_ (double n, double[] p) {
      double corr[][] = new double[p.length][p.length];
      double sumPi= 0.0;
      double p0;

      for (int i = 0; i < p.length;i++)
         sumPi += p[i];
      p0 = 1.0 - sumPi;

      for (int i = 0; i < p.length; i++) {
         for (int j = 0; j < p.length; j++)
            corr[i][j] = Math.sqrt(p[i] * p[j] /((p0 + p[i]) * (p0 + p[j])));
         corr[i][i] = 1.0;
      }
      return corr;
   }

/**
 * Computes the correlation matrix of the negative multinomial distribution
 * with parameters @f$n@f$ and (@f$p_1@f$, …, @f$p_d@f$).
 */
public static double[][] getCorrelation (double n, double[] p) {
      verifParam (n, p);
      return getCorrelation_ (n, p);
   }

   /**
    * Estimates and returns the parameters [@f$\hat{n}@f$,
    * @f$\hat{p}_1@f$, …, @f$\hat{p}_d@f$] of the negative multinomial
    * distribution using the maximum likelihood method. It uses the
    * @f$m@f$ observations of @f$d@f$ components in table
    * <tt>x[</tt>@f$i@f$<tt>][</tt>@f$j@f$<tt>]</tt>, @f$i = 0, 1, …,
    * m-1@f$ and @f$j = 0, 1, …, d-1@f$.  The equations of the maximum
    * likelihood are defined in @cite tJOH69a&thinsp;:
    * @f{align*}{
    *    \sum_{s=1}^M \frac{F_s}{(\hat{n} + s - 1)} 
    *    & 
    *    = 
    *    \ln\left(1 + \frac{1}{\hat{n} m} \sum_{j=1}^m \Upsilon_j \right)
    *    \\ 
    *    p_i 
    *    & 
    *    = 
    *    \frac{\lambda_i}{1 + \sum_{j=1}^d \lambda_j} \qquad\mbox{for } i=1, …,d
    * @f}
    * where
    * @f{align*}{
    *    \lambda_i 
    *    & 
    *    = 
    *    \frac{\sum_{j=1}^m X_{i,j}}{\hat{n} m} \qquad\mbox{for } i=1, …,d
    *    \\ 
    *   \Upsilon_j 
    *    & 
    *    = 
    *    \sum_{i=1}^d X_{i,j} \qquad\mbox{for } j=1, …,m
    *    \\ 
    *   F_s 
    *    & 
    *    = 
    *    \frac{1}{m} \sum_{j=1}^m \mbox{\textbf{1}}\{\Upsilon_j \ge s\} \qquad\mbox{for } s=1, …,M
    *    \\ 
    *   M 
    *    & 
    *    = 
    *    \max_j \{\Upsilon_j\}
    * @f}
    * @param x            the list of observations used to evaluate
    *                      parameters
    *  @param m            the number of observations used to evaluate
    *                      parameters
    *  @param d            the dimension of each vector
    *  @return returns the parameters [@f$\hat{n}@f$, @f$\hat{p}_1@f$, …,
    * @f$\hat{p}_d@f$]
    */
   public static double[] getMLE (int x[][], int m, int d) {
      int ups[] = new int[m];
      double mean[] = new double[d];

      int i, j, l;
      int M;
      int prop;

      // Initialization
      for (i = 0; i < d; i++)
         mean[i] = 0;

      // Ups_j = Sum_k x_ji
      // mean_i = Sum_m x_ji / m
      for (j = 0; j < m; j++) {
         ups[j] = 0;
         for (i = 0; i < d; i++) {
            ups[j] += x[j][i];
            mean[i] += x[j][i];
         }
      }
      for (i = 0; i < d; i++)
         mean[i] /= m;

/*
      double var = 0.0;
      if (d > 1) {
         // Calcule la covariance 0,1
         for (j = 0; j < m; j++)
            var += (x[j][0] - mean[0])*(x[j][1] - mean[1]);
         var /= m;
      } else {
         // Calcule la variance 0
         for (j = 0; j < m; j++)
            var += (x[j][0] - mean[0])*(x[j][0] - mean[0]);
         var /= m;
      }
*/

      // M = Max(Ups_j)
      M = ups[0];
      for (j = 1; j < m; j++)
         if (ups[j] > M)
            M = ups[j];

      if (M >= Integer.MAX_VALUE)
         throw new IllegalArgumentException("n/p_i too large");

      double Fl[] = new double[M];
      for (l = 0; l < M; l++) {
         prop = 0;
         for (j = 0; j < m; j++)
            if (ups[j] > l)
               prop++;

         Fl[l] = (double) prop / (double) m;
      }

/*
      // Estime la valeur initiale de n pour brentDekker: pourrait
      // accélérer brentDekker (gam0/1000, gam0*1000, f, 1e-5).
      // Reste à bien tester.
      if (d > 1) {
         double gam0 = mean[0] * mean[1] / var;
         System.out.println ("gam0 = " + gam0);
      } else {
         double t = var/mean[0] - 1.0;
         double gam0 = mean[0] / t;
         System.out.println ("gam0 = " + gam0);
      }
*/
      double parameters[] = new double[d + 1];
      Function f = new Function (m, (int)M, ups, Fl);
      parameters[0] = RootFinder.brentDekker (1e-9, 1e9, f, 1e-5);

      double lambda[] = new double[d];
      double sumLambda = 0.0;
      for (i = 0; i < d; i++) {
         lambda[i] = mean[i] / parameters[0];
         sumLambda += lambda[i];
      }

      for (i = 0; i < d; i++) {
         parameters[i + 1] = lambda[i] / (1.0 + sumLambda);
         if (parameters[i + 1] > 1.0)
            throw new IllegalArgumentException("p_i > 1");
      }

      return parameters;
   }

   /**
    * Estimates and returns the parameter @f$\nu= 1/\hat{n}@f$ of the
    * negative multinomial distribution using the maximum likelihood
    * method. It uses the @f$m@f$ observations of @f$d@f$ components in
    * table <tt>x[</tt>@f$i@f$<tt>][</tt>@f$j@f$<tt>]</tt>, @f$i = 0, 1,
    * …, m-1@f$ and @f$j = 0, 1, …, d-1@f$.  The equation of maximum
    * likelihood is defined as:
    * @f[
    *   \sum_{s=1}^M \frac{\nu F_s}{(1 + \nu(s - 1))} = \ln\left(1 + \frac{\nu}{m} \sum_{j=1}^m \Upsilon_j \right)\\
    * @f]
    * where the symbols are defined as in  #getMLE(int[][],int,int).
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param m            the number of observations used to evaluate
    *                      parameters
    *  @param d            the dimension of each vector
    *  @return returns the parameter @f$1/\hat{n}@f$
    */
   public static double getMLEninv (int x[][], int m, int d) {
      int ups[] = new int[m];
      double mean[] = new double[d];
      int i, j, l;
      int M;
      int prop;

      // Initialization
      for (i = 0; i < d; i++)
         mean[i] = 0;

      // Ups_j = Sum_k x_ji
      // mean_i = Sum_m x_ji / m
      for (j = 0; j < m; j++) {
         ups[j] = 0;
         for (i = 0; i < d; i++) {
            ups[j] += x[j][i];
            mean[i] += x[j][i];
         }
      }
      for (i = 0; i < d; i++)
         mean[i] /= m;

      // M = Max(Ups_j)
      M = ups[0];
      for (j = 1; j < m; j++)
         if (ups[j] > M)
            M = ups[j];

      if (M >= Integer.MAX_VALUE)
         throw new IllegalArgumentException("n/p_i too large");

      double Fl[] = new double[M];
      for (l = 0; l < M; l++) {
         prop = 0;
         for (j = 0; j < m; j++)
            if (ups[j] > l)
               prop++;

         Fl[l] = (double) prop / (double) m;
      }

      FuncInv f = new FuncInv (m, M, ups, Fl);
      double nu = RootFinder.brentDekker (1.0e-8, 1.0e8, f, 1e-5);
      return nu;
   }

   /**
    * Returns the parameter @f$n@f$ of this object.
    */
   public double getGamma() {
      return n;
   }

   /**
    * Returns the parameters (@f$p_1@f$, …, @f$p_d@f$) of this object.
    */
   public double[] getP() {
      return p;
   }

   /**
    * Sets the parameters @f$n@f$ and (@f$p_1@f$, …, @f$p_d@f$) of this
    * object.
    */
   public void setParams (double n, double p[]) {
      if (n <= 0.0)
         throw new IllegalArgumentException("n <= 0");

      this.n = n;
      this.dimension = p.length;
      this.p = new double[dimension];

      double sumPi = 0.0;
      for (int i = 0; i < dimension; i++) {
         if ((p[i] < 0) || (p[i] >= 1))
            throw new IllegalArgumentException("p is not a probability vector");

         sumPi += p[i];
         this.p[i] = p[i];
      }

      if (sumPi >= 1.0)
         throw new IllegalArgumentException("p is not a probability vector");
   }

}