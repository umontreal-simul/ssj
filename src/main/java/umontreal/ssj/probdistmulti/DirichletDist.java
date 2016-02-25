/*
 * Class:        DirichletDist
 * Description:  Dirichlet distribution
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
import optimization.*;

/**
 * Implements the abstract class  @ref ContinuousDistributionMulti for the
 * *Dirichlet* distribution with parameters
 * (@f$\alpha_1@f$,…,@f$\alpha_d@f$), @f$\alpha_i > 0@f$. The probability
 * density is
 * @anchor REF_probdistmulti_DirichletDist_eq_fDirichlet
 * @f[
 *   f(x_1,…, x_d) = \frac{\Gamma(\alpha_0) \prod_{i=1}^d x_i^{\alpha_i - 1}}{\prod_{i=1}^d \Gamma(\alpha_i)} \tag{fDirichlet}
 * @f]
 * where @f$x_i \ge0@f$, @f$\sum_{i=1}^d x_i = 1@f$, @f$\alpha_0 =
 * \sum_{i=1}^d \alpha_i@f$, and @f$\Gamma@f$ is the Gamma function.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdistmulti_continuous
 */
public class DirichletDist extends ContinuousDistributionMulti {
   private static final double LOGMIN = -709.1;    // Log(MIN_DOUBLE/2)
   protected double[] alpha;

   private static class Optim implements Uncmin_methods
   {
      double[] logP;
      int n;
      int k;

      public Optim (double[] logP, int n) {
         this.n = n;
         this.k = logP.length;
         this.logP = new double[k];
         System.arraycopy (logP, 0, this.logP, 0, k);
      }

      public double f_to_minimize (double[] alpha) {
         double sumAlpha = 0.0;
         double sumLnGammaAlpha = 0.0;
         double sumAlphaLnP = 0.0;

         for (int i = 1; i < alpha.length; i++) {
            if (alpha[i] <= 0.0)
               return 1.0e200;

            sumAlpha += alpha[i];
            sumLnGammaAlpha += Num.lnGamma (alpha[i]);
            sumAlphaLnP += ((alpha[i] - 1.0) * logP[i - 1]);
         }

         return (- n * (Num.lnGamma (sumAlpha) - sumLnGammaAlpha + sumAlphaLnP));
      }

      public void gradient (double[] alpha, double[] g)
      {
      }

      public void hessian (double[] alpha, double[][] h)
      {
      }
   }
   public DirichletDist (double[] alpha) {
      setParams (alpha);
   }


   public double density (double[] x) {
      return density_ (alpha, x);
   }

   public double[] getMean() {
      return getMean_ (alpha);
   }

   public double[][] getCovariance() {
      return getCovariance_ (alpha);
   }

   public double[][] getCorrelation () {
      return getCorrelation_ (alpha);
   }

   private static void verifParam (double[] alpha) {

      for (int i = 0; i < alpha.length;i++)
      {
         if (alpha[i] <= 0)
            throw new IllegalArgumentException("alpha[" + i + "] <= 0");
      }
   }

   private static double density_ (double[] alpha, double[] x) {
      double alpha0 = 0.0;
      double sumLnGamma = 0.0;
      double sumAlphaLnXi = 0.0;

      if (alpha.length != x.length)
         throw new IllegalArgumentException ("alpha and x must have the same dimension");

      for (int i = 0; i < alpha.length; i++) {
         alpha0 += alpha[i];
         sumLnGamma += Num.lnGamma (alpha[i]);
         if (x[i] <= 0.0 || x[i] >= 1.0)
            return 0.0;
         sumAlphaLnXi += (alpha[i] - 1.0) * Math.log (x[i]);
      }

      return Math.exp (Num.lnGamma (alpha0) - sumLnGamma + sumAlphaLnXi);
   }

/**
 * Computes the density (
 * {@link REF_probdistmulti_DirichletDist_eq_fDirichlet
 * fDirichlet} ) of the Dirichlet distribution with parameters
 * (@f$\alpha_1@f$, …, @f$\alpha_d@f$).
 */
public static double density (double[] alpha, double[] x) {
      verifParam (alpha);
      return density_ (alpha, x);
   }


   private static double[][] getCovariance_ (double[] alpha) {
      double[][] cov = new double[alpha.length][alpha.length];
      double alpha0 = 0.0;

      for (int i =0; i < alpha.length; i++)
         alpha0 += alpha[i];

      for (int i = 0; i < alpha.length; i++) {
         for (int j = 0; j < alpha.length; j++)
            cov[i][j] = - (alpha[i] * alpha[j]) / (alpha0 * alpha0 * (alpha0 + 1.0));

         cov[i][i] = (alpha[i] / alpha0) * (1.0 - alpha[i] / alpha0) / (alpha0 + 1.0);
      }

      return cov;
   }

/**
 * Computes the covariance matrix of the Dirichlet distribution with
 * parameters (@f$\alpha_1@f$, …, @f$\alpha_d@f$).
 */
public static double[][] getCovariance (double[] alpha) {
      verifParam (alpha);

      return getCovariance_ (alpha);
   }


   private static double[][] getCorrelation_ (double[] alpha) {
      double corr[][] = new double[alpha.length][alpha.length];
      double alpha0 = 0.0;

      for (int i =0; i < alpha.length; i++)
         alpha0 += alpha[i];

      for (int i = 0; i < alpha.length; i++) {
         for (int j = 0; j < alpha.length; j++)
            corr[i][j] = - Math.sqrt ((alpha[i] * alpha[j]) /
                                      ((alpha0 - alpha[i]) * (alpha0 - alpha[j])));
         corr[i][i] = 1.0;
      }
      return corr;
   }

/**
 * Computes the correlation matrix of the Dirichlet distribution with
 * parameters (@f$\alpha_1@f$, …, @f$\alpha_d@f$).
 */
public static double[][] getCorrelation (double[] alpha) {
      verifParam (alpha);
      return getCorrelation_ (alpha);
   }

   /**
    * Estimates the parameters [@f$\hat{\alpha_1},…,\hat{\alpha_d}@f$]
    * of the Dirichlet distribution using the maximum likelihood method.
    * It uses the @f$n@f$ observations of @f$d@f$ components in table
    * @f$x[i][j]@f$, @f$i = 0, 1, …, n-1@f$ and @f$j = 0, 1, …, d-1@f$.
    * The equations of the maximum likelihood are defined in
    * @cite ccAVR04a&thinsp; (Technical appendix)
    * @f{align*}{
    *    L(\hat{\alpha}_1,\hat{\alpha}_2,…,\hat{\alpha}_k) 
    *    & 
    *    = 
    *    n \left( G(\alpha_0) - \sum_{i=1}^k G(\hat{\alpha}_i) \right) + \sum_{i=1}^k (\hat{\alpha}_i - 1) Z_i
    * @f}
    * where @f$G@f$ is the logarithm of the gamma function and
    * @f{align*}{
    *    \alpha_0 
    *    & 
    *    = 
    *    \sum_{i=1}^k \hat{\alpha}_i
    *    \\ 
    *   Z_i 
    *    & 
    *    = 
    *    \sum_{j=1}^n \ln(X_{i,j}) \qquad\mbox{for }i=1,…,k.
    * @f}
    * @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @param d            the dimension of each vector
    *  @return returns the parameter
    * [@f$\hat{\alpha_1},…,\hat{\alpha_d}@f$]
    */
   public static double[] getMLE (double[][] x, int n, int d) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      if (d <= 0)
         throw new IllegalArgumentException ("d <= 0");

      double[] logP = new double[d];
      double mean[] = new double[d];
      double var[] = new double[d];
      int i;
      int j;
      for (i = 0; i < d; i++) {
         logP[i] = 0.0;
         mean[i] = 0.0;
      }

      for (i = 0; i < n; i++) {
         for (j = 0; j < d; j++) {
            if (x[i][j] > 0.)
               logP[j] += Math.log (x[i][j]);
            else
               logP[j] += LOGMIN;
            mean[j] += x[i][j];
         }
      }

      for (i = 0; i < d; i++) {
         logP[i] /= (double) n;
         mean[i] /= (double) n;
      }

      double sum = 0.0;
      for (j = 0; j < d; j++) {
         sum = 0.0;
         for (i = 0; i < n; i++)
            sum += (x[i][j] - mean[j]) * (x[i][j] - mean[j]);
         var[j] = sum / (double) n;
      }

      double alpha0 = (mean[0] * (1.0 - mean[0])) / var[0] - 1.0;
      Optim system = new Optim (logP, n);

      double[] parameters = new double[d];
      double[] xpls = new double[d + 1];
      double[] alpha = new double[d + 1];
      double[] fpls = new double[d + 1];
      double[] gpls = new double[d + 1];
      int[] itrcmd = new int[2];
      double[][] a = new double[d + 1][d + 1];
      double[] udiag = new double[d + 1];

      for (i = 1; i <= d; i++)
         alpha[i] = mean[i - 1] * alpha0;

      Uncmin_f77.optif0_f77 (d, alpha, system, xpls, fpls, gpls, itrcmd, a, udiag);

      for (i = 0; i < d; i++)
         parameters[i] = xpls[i+1];

      return parameters;
   }


   private static double[] getMean_ (double[] alpha) {
      double alpha0 = 0.0;
      double[] mean = new double[alpha.length];

      for (int i = 0; i < alpha.length;i++)
         alpha0 += alpha[i];

      for (int i = 0; i < alpha.length; i++)
         mean[i] = alpha[i] / alpha0;

      return mean;
   }

/**
 * Computes the mean @f$E[X] = \alpha_i / \alpha_0@f$ of the Dirichlet
 * distribution with parameters (@f$\alpha_1@f$, …, @f$\alpha_d@f$), where
 * @f$\alpha_0 = \sum_{i=1}^d \alpha_i@f$.
 */
public static double[] getMean (double[] alpha) {
      verifParam (alpha);
      return getMean_ (alpha);
   }

   /**
    * Returns the parameters (@f$\alpha_1@f$, …, @f$\alpha_d@f$) of this
    * object.
    */
   public double[] getAlpha() {
      return alpha;
   }

   /**
    * Returns the @f$i@f$th component of the alpha vector.
    */
   public double getAlpha (int i) {
      return alpha[i];
   }

   /**
    * Sets the parameters (@f$\alpha_1@f$, …, @f$\alpha_d@f$) of this
    * object.
    */
   public void setParams (double[] alpha) {
      this.dimension = alpha.length;
      this.alpha = new double[dimension];
      for (int i = 0; i < dimension; i++) {
         if (alpha[i] <= 0)
            throw new IllegalArgumentException("alpha[" + i + "] <= 0");
         this.alpha[i] = alpha[i];
      }
   }

}