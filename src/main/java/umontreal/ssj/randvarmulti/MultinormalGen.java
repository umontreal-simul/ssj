/*
 * Class:        MultinormalGen
 * Description:  multivariate normal random variable generator
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
package umontreal.ssj.randvarmulti;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.rng.RandomStream;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.CholeskyDecomposition;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

/**
 * Extends  @ref RandomMultivariateGen for a *multivariate normal* (or
 * *multinormal*) distribution @cite tJOH72a&thinsp;. The @f$d@f$-dimensional
 * multivariate normal distribution with mean vector
 * @f$\boldsymbol{\mu}\in\mathbb{R}^d@f$ and (symmetric positive-definite)
 * covariance matrix @f$\boldsymbol{\Sigma}@f$, denoted
 * @f$N(\boldsymbol{\mu}, \boldsymbol{\Sigma})@f$, has density
 * @f[
 *   f(\mathbf{X})=\frac{1}{\sqrt{(2\pi)^d\det(\boldsymbol{\Sigma})}} \exp\left(-(\mathbf{X}- \boldsymbol{\mu})^{\!\mathsf{t}}\boldsymbol{\Sigma}^{-1}(\mathbf{X}- \boldsymbol{\mu})/2\right),
 * @f]
 * for all @f$\mathbf{X}\in\mathbb{R}^d@f$, and
 * @f$\mathbf{X}^{\mathsf{t}}@f$ is the transpose vector of @f$\mathbf{X}@f$.
 * If @f$\mathbf{Z}\sim N(\boldsymbol{0}, \mathbf{I})@f$ where
 * @f$\mathbf{I}@f$ is the identity matrix, @f$\mathbf{Z}@f$ is said to have
 * the *standard multinormal* distribution.
 *
 * For the special case @f$d=2@f$, if the random vector @f$\mathbf{X}= (X_1,
 * X_2)^{\mathsf{t}}@f$ has a bivariate normal distribution, then it has mean
 * @f$\boldsymbol{\mu}= (\mu_1, \mu_2)^{\mathsf{t}}@f$, and covariance
 * matrix
 * @f[
 *   \boldsymbol{\Sigma}= \left[\begin{array}{cc}
 *    \sigma_1^2 
 *    & 
 *    \rho\sigma_1\sigma_2 
 *    \\ 
 *   \rho\sigma_1\sigma_2 
 *    & 
 *   \sigma_2^2 
 *   \end{array}\right]
 * @f]
 * if and only if @f$\mathrm{Var}[X_1] = \sigma_1^2@f$, @f$\mathrm{Var}[X_2]
 * = \sigma_2^2@f$, and the linear correlation between @f$X_1@f$ and
 * @f$X_2@f$ is @f$\rho@f$, where @f$\sigma_1 > 0@f$, @f$\sigma_2 > 0@f$,
 * and @f$-1 \le\rho\le1@f$.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class MultinormalGen extends RandomMultivariateGen {
   protected double[] mu;
   protected DoubleMatrix2D sigma;
   protected DoubleMatrix2D sqrtSigma;
   protected double[] temp;
   protected static final double MYINF = 37.54;


   private void initMN (NormalGen gen1, double[] mu, int d) {
      if (gen1 == null)
         throw new NullPointerException ("gen1 is null");

      if (gen1.getMu() != 0.0)
         throw new IllegalArgumentException ("mu != 0");
      if (gen1.getSigma() != 1.0)
         throw new IllegalArgumentException ("sigma != 1");
/*
      NormalDist dist = (NormalDist) gen1.getDistribution();
      if (dist != null) {
         if (dist.getMu() != 0.0)
            throw new IllegalArgumentException ("mu != 0");
         if (dist.getSigma() != 1.0)
            throw new IllegalArgumentException ("sigma != 1");
      }
      dist = null;
*/
      this.gen1 = gen1;

      if (mu == null) {    // d is the dimension
         dimension = d;
         this.mu = new double[d];
      } else {      // d is unused
         dimension = mu.length;
         this.mu = (double[])mu.clone();
      }
      temp = new double[dimension];
     }

   /**
    * Constructs a generator with the standard multinormal distribution
    * (with @f$\boldsymbol{\mu}=\boldsymbol{0}@f$ and
    * @f$\boldsymbol{\Sigma}= \mathbf{I}@f$) in @f$d@f$ dimensions. Each
    * vector @f$\mathbf{Z}@f$ will be generated via @f$d@f$ successive
    * calls to `gen1`, which must be a *standard normal* generator.
    *  @param gen1         the one-dimensional generator
    *  @param d            the dimension of the generated vectors
    *  @exception IllegalArgumentException if the one-dimensional normal
    * generator uses a normal distribution with @f$\mu@f$ not equal to 0,
    * or @f$\sigma@f$ not equal to 1.
    *  @exception IllegalArgumentException if `d` is negative.
    *  @exception NullPointerException if `gen1` is `null`.
    */
   public MultinormalGen (NormalGen gen1, int d) {
      initMN (gen1, null, d);
      sigma = new DenseDoubleMatrix2D (d, d);
      sqrtSigma = new DenseDoubleMatrix2D (d, d);
      for (int i = 0; i < d; i++) {
         sigma.setQuick (i, i, 1.0);
         sqrtSigma.setQuick (i, i, 1.0);
      }
   }

   /**
    * Constructs a multinormal generator with mean vector `mu` and
    * covariance matrix `sigma`. The mean vector must have the same length
    * as the dimensions of the covariance matrix, which must be symmetric
    * and positive-definite. If any of the above conditions is violated,
    * an exception is thrown. The vector @f$\mathbf{Z}@f$ is generated by
    * calling @f$d@f$ times the generator `gen1`, which must be *standard
    * normal*.
    *  @param gen1         the one-dimensional generator
    *  @param mu           the mean vector.
    *  @param sigma        the covariance matrix.
    *  @exception NullPointerException if any argument is `null`.
    *  @exception IllegalArgumentException if the length of the mean
    * vector is incompatible with the dimensions of the covariance matrix.
    */
   protected MultinormalGen (NormalGen gen1, double[] mu,
                             DoubleMatrix2D sigma) {
      initMN (gen1, mu, -1);
      this.sigma = sigma.copy();
   }

   /**
    * Equivalent to
    * {@link #MultinormalGen(NormalGen,double[],DoubleMatrix2D)
    * MultinormalGen(gen1, mu, new DenseDoubleMatrix2D (sigma))}.
    */
   protected MultinormalGen (NormalGen gen1, double[] mu, double[][] sigma) {
      initMN (gen1, mu, -1);
      this.sigma = new DenseDoubleMatrix2D (sigma);
   }

   /**
    * Returns the mean vector used by this generator.
    *  @return the current mean vector.
    */
   public double[] getMu() {
      return mu;
   }

   /**
    * Returns the @f$i@f$-th component of the mean vector for this
    * generator.
    *  @param i            the index of the required component.
    *  @return the value of @f$\mu_i@f$.
    *
    *  @exception ArrayIndexOutOfBoundsException if `i` is negative or
    * greater than or equal to  #getDimension.
    */
   public double getMu (int i) {
      return mu[i];
   }

   /**
    * Sets the mean vector to `mu`.
    *  @param mu           the new mean vector.
    *  @exception NullPointerException if `mu` is `null`.
    *  @exception IllegalArgumentException if the length of `mu` does not
    * correspond to  #getDimension.
    */
   public void setMu (double[] mu) {
      if (mu.length != this.mu.length)
         throw new IllegalArgumentException
            ("Incompatible length of mean vector");
      this.mu = mu;
   }

   /**
    * Sets the @f$i@f$-th component of the mean vector to `mui`.
    *  @param i            the index of the modified component.
    *  @param mui          the new value of @f$\mu_i@f$.
    *  @exception ArrayIndexOutOfBoundsException if `i` is negative or
    * greater than or equal to  #getDimension.
    */
   public void setMu (int i, double mui) {
      mu[i] = mui;
   }

   /**
    * Returns the covariance matrix @f$\boldsymbol{\Sigma}@f$ used by this
    * generator.
    *  @return the used covariance matrix.
    */
   public DoubleMatrix2D getSigma() {
      return sigma.copy();
   }

   /**
    * Generates a point from this multinormal distribution.
    *  @param p            the array to be filled with the generated point
    */
   public void nextPoint (double[] p) {
      int n = dimension;
      for (int i = 0; i < n; i++) {
         temp[i] = gen1.nextDouble();
         if (temp[i] == Double.NEGATIVE_INFINITY)
            temp[i] = -MYINF;
         if (temp[i] == Double.POSITIVE_INFINITY)
            temp[i] = MYINF;
      }
      for (int i = 0; i < n; i++) {
         p[i] = 0;
         for (int c = 0; c < n; c++)
            p[i] += sqrtSigma.getQuick (i, c)*temp[c];
         p[i] += mu[i];
      }
   }

}