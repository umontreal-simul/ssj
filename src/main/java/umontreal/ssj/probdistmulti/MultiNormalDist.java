/*
 * Class:        MultiNormalDist
 * Description:  multinormal distribution
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
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

/**
 * Implements the abstract class  @ref ContinuousDistributionMulti for the
 * *multinormal* distribution with mean vector @f$\boldsymbol{\mu}@f$ and
 * covariance matrix @f$\boldsymbol{\Sigma}@f$. The probability density is
 * @anchor REF_probdistmulti_MultiNormalDist_eq_fMultinormal
 * @f[
 *   f(\mathbf{x}) = \frac{1}{\sqrt{(2\pi)^d \det\boldsymbol{\Sigma}}} \exp\left(-\frac{1}{2}(\mathbf{x}- \boldsymbol{\mu})^T \boldsymbol{\Sigma}^{-1} (\mathbf{x}- \boldsymbol{\mu})\right) \tag{fMultinormal}
 * @f]
 * where @f$\mathbf{x}= (x_1,…,x_d)@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdistmulti_continuous
 */
public class MultiNormalDist extends ContinuousDistributionMulti {
   protected int dim;
   protected double[] mu;
   protected DoubleMatrix2D sigma;
   protected DoubleMatrix2D invSigma;

   protected static Algebra algebra = new Algebra();
   public MultiNormalDist (double[] mu, double[][] sigma) {
      setParams (mu, sigma);
   }


   public double density (double[] x) {
      double sum = 0.0;

      if (invSigma == null)
         invSigma = algebra.inverse(sigma);

      double[] temp = new double[mu.length];
      for (int i = 0; i < mu.length; i++)
      {
         sum = 0.0;
         for (int j = 0; j < mu.length; j++)
            sum += ((x[j] - mu[j]) * invSigma.getQuick (j, i));
         temp[i] = sum;
      }

      sum = 0.0;
      for (int i = 0; i < mu.length; i++)
         sum += temp[i] * (x[i] - mu[i]);

      return (Math.exp(-0.5 * sum) / Math.sqrt (Math.pow (2 * Math.PI, mu.length) * algebra.det (sigma)));
   }

   public double[] getMean() {
      return mu;
   }

   public double[][] getCovariance() {
      return sigma.toArray();
   }

   public double[][] getCorrelation () {
      return getCorrelation_ (mu, sigma.toArray());
   }

/**
 * Computes the density (
 * {@link REF_probdistmulti_MultiNormalDist_eq_fMultinormal
 * fMultinormal} ) of the multinormal distribution with parameters
 * @f$\boldsymbol{\mu}=@f$ `mu` and @f$\boldsymbol{\Sigma}=@f$ `sigma`,
 * evaluated at `x`.
 */
public static double density (double[] mu, double[][] sigma, double[] x) {
      double sum = 0.0;
      DoubleMatrix2D sig;
      DoubleMatrix2D inv;

      if (sigma.length != sigma[0].length)
         throw new IllegalArgumentException ("sigma must be a square matrix");
      if (mu.length != sigma.length)
         throw new IllegalArgumentException ("mu and sigma must have the same dimension");

      sig = new DenseDoubleMatrix2D (sigma);
      inv = algebra.inverse (sig);

      double[] temp = new double[mu.length];
      for (int i = 0; i < mu.length; i++)
      {
         sum = 0.0;
         for (int j = 0; j < mu.length; j++)
            sum += ((x[j] - mu[j]) * inv.getQuick (j, i));
         temp[i] = sum;
      }

      sum = 0.0;
      for (int i = 0; i < mu.length; i++)
         sum += temp[i] * (x[i] - mu[i]);

      return (Math.exp(-0.5 * sum) / Math.sqrt (Math.pow (2 * Math.PI, mu.length) * algebra.det (sig)));
   }

   /**
    * Returns the dimension @f$d@f$ of the distribution.
    */
   public int getDimension() {
      return dim;
   }

   /**
    * Returns the mean @f$E[\mathbf{X}] = \boldsymbol{\mu}@f$ of the
    * multinormal distribution with parameters @f$\boldsymbol{\mu}@f$ and
    * @f$\boldsymbol{\Sigma}@f$.
    */
   public static double[] getMean (double[] mu, double[][] sigma) {
      if (sigma.length != sigma[0].length)
         throw new IllegalArgumentException ("sigma must be a square matrix");
      if (mu.length != sigma.length)
         throw new IllegalArgumentException ("mu and sigma must have the same dimension");

      return mu;
   }

   /**
    * Computes the covariance matrix of the multinormal distribution with
    * parameters @f$\boldsymbol{\mu}@f$ and @f$\boldsymbol{\Sigma}@f$.
    */
   public static double[][] getCovariance (double[] mu, double[][] sigma) {
      if (sigma.length != sigma[0].length)
         throw new IllegalArgumentException ("sigma must be a square matrix");
      if (mu.length != sigma.length)
         throw new IllegalArgumentException ("mu and sigma must have the same dimension");

      return sigma;
   }


   private static double[][] getCorrelation_ (double[] mu, double[][] sigma) {
      double corr[][] = new double[mu.length][mu.length];

      for (int i = 0; i < mu.length; i++) {
         for (int j = 0; j < mu.length; j++)
            corr[i][j] = - sigma[i][j] / Math.sqrt (sigma[i][i] * sigma[j][j]);
         corr[i][i] = 1.0;
      }
      return corr;
   }

/**
 * Computes the correlation matrix of the multinormal distribution with
 * parameters @f$\boldsymbol{\mu}@f$ and @f$\boldsymbol{\Sigma}@f$).
 */
public static double[][] getCorrelation (double[] mu, double[][] sigma) {
      if (sigma.length != sigma[0].length)
         throw new IllegalArgumentException ("sigma must be a square matrix");
      if (mu.length != sigma.length)
         throw new IllegalArgumentException ("mu and sigma must have the same dimension");

      return getCorrelation_ (mu, sigma);
   }

   /**
    * Estimates the parameters @f$\boldsymbol{\mu}@f$ of the multinormal
    * distribution using the maximum likelihood method. It uses the
    * @f$n@f$ observations of @f$d@f$ components in table @f$x[i][j]@f$,
    * @f$i = 0, 1, …, n-1@f$ and @f$j = 0, 1, …, d-1@f$.
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @param d            the dimension of each observation
    *  @return returns the parameters
    * [@f$\boldsymbol{\mu}_1@f$,…,@f$\boldsymbol{\mu}_d@f$]
    */
   public static double[] getMLEMu (double[][] x, int n, int d) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      if (d <= 0)
         throw new IllegalArgumentException ("d <= 0");

      double[] parameters = new double[d];
      for (int i = 0; i < parameters.length; i++)
         parameters[i] = 0.0;

      for (int i = 0; i < n; i++)
         for (int j = 0; j < d; j++)
            parameters[j] += x[i][j];

      for (int i = 0; i < parameters.length; i++)
         parameters[i] = parameters[i] / (double) n;

      return parameters;
   }

   /**
    * Estimates the parameters @f$\boldsymbol{\Sigma}@f$ of the
    * multinormal distribution using the maximum likelihood method. It
    * uses the @f$n@f$ observations of @f$d@f$ components in table
    * @f$x[i][j]@f$, @f$i = 0, 1, …, n-1@f$ and @f$j = 0, 1, …, d-1@f$.
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @param d            the dimension of each observation
    *  @return returns the covariance matrix @f$\boldsymbol{\Sigma}@f$
    */
   public static double[][] getMLESigma (double[][] x, int n, int d) {
      double sum = 0.0;

      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      if (d <= 0)
         throw new IllegalArgumentException ("d <= 0");

      double[] mean = getMLEMu (x, n, d);
      double[][] parameters = new double[d][d];
      for (int i = 0; i < parameters.length; i++)
         for (int j = 0; j < parameters.length; j++)
            parameters[i][j] = 0.0;

      for (int i = 0; i < parameters.length; i++)
      {
         for (int j = 0; j < parameters.length; j++)
         {
            sum = 0.0;
            for (int t = 0; t < n; t++)
               sum += (x[t][i] - mean[i]) * (x[t][j] - mean[j]);
            parameters[i][j] = sum / (double) n;
         }
      }

      return parameters;
   }

   /**
    * Returns the parameter @f$\boldsymbol{\mu}@f$ of this object.
    */
   public double[] getMu() {
      return mu;
   }

   /**
    * Returns the @f$i@f$-th component of the parameter
    * @f$\boldsymbol{\mu}@f$ of this object.
    */
   public double getMu (int i) {
      return mu[i];
   }

   /**
    * Returns the parameter @f$\boldsymbol{\Sigma}@f$ of this object.
    */
   public double[][] getSigma() {
      return sigma.toArray();
   }

   /**
    * Sets the parameters @f$\boldsymbol{\mu}@f$ and
    * @f$\boldsymbol{\Sigma}@f$ of this object.
    */
   public void setParams (double[] mu, double[][] sigma) {
      if (sigma.length != sigma[0].length)
         throw new IllegalArgumentException ("sigma must be a square matrix");
      if (mu.length != sigma.length)
         throw new IllegalArgumentException ("mu and sigma must have the same dimension");

      this.mu = new double[mu.length];
      this.dimension = mu.length;
      System.arraycopy(mu, 0, this.mu, 0, mu.length);
      this.sigma = new DenseDoubleMatrix2D (sigma);

      invSigma = null;
   }

}