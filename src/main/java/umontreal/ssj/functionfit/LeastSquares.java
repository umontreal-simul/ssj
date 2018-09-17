/*
 * Class:        LeastSquares
 * Description:  General linear regression with the least squares method
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2013  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        April 2013
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
package umontreal.ssj.functionfit;

import java.io.Serializable;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.QRDecomposition;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import cern.colt.matrix.linalg.Algebra;

/**
 * This class implements different *linear regression* models, using the
 * least squares method to estimate the regression coefficients. Given input
 * data @f$x_{ij}@f$ and response @f$y_i@f$, one want to find the
 * coefficients @f$\beta_j@f$ that minimize the residuals of the form (using
 * matrix notation)
 * @f[
 *   r = \min_{\beta}\| Y - X\beta\|_2,
 * @f]
 * where the @f$L_2@f$ norm is used. Particular cases are
 * @f[
 *   r = \min_{\beta}\sum_i \left(y_i - \beta_0 - \sum_{j=1}^k \beta_j x_{ij}\right)^2.
 * @f]
 * for @f$k@f$ regressor variables @f$x_j@f$. The well-known case of the
 * single variable @f$x@f$ is
 * @f[
 *   r = \min_{\alpha,\beta} \sum_i \left(y_i - \alpha- \beta x_i\right)^2.
 * @f]
 * Sometimes, one wants to use a basis of general functions @f$\psi_j(t)@f$
 * with a minimization of the form
 * @f[
 *   r = \min_{\beta}\sum_i \left(y_i - \sum_{j=1}^k \beta_j\psi_j(t_i)\right)^2.
 * @f]
 * For example, we could have @f$\psi_j(t) = e^{-\lambda_j t}@f$ or some
 * other functions. In that case, one has to choose the points @f$t_i@f$ at
 * which to compute the basis functions, and use a method below with
 * @f$x_{ij} = \psi_j(t_i)@f$.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class LeastSquares {

   private static double[] solution (DoubleMatrix2D X, DoubleMatrix2D Y, int k) {
      // Solve X * Beta = Y for Beta
      // Only the first column of Y is used
      // k is number of beta coefficients

      QRDecomposition qr = new QRDecomposition(X);

      if (qr.hasFullRank()) {
         DoubleMatrix2D B = qr.solve(Y);
         return B.viewColumn(0).toArray();

      } else {
         DoubleMatrix1D Y0 = Y.viewColumn(0);   // first column of Y
         SingularValueDecomposition svd = new SingularValueDecomposition(X);
         DoubleMatrix2D S = svd.getS();
         DoubleMatrix2D V = svd.getV();
         DoubleMatrix2D U = svd.getU();
         Algebra alg = new Algebra();
         DoubleMatrix2D Ut = alg.transpose(U);
         DoubleMatrix1D g = alg.mult(Ut, Y0);    // Ut*Y0

         for (int j = 0; j < k; j++) {
            // solve S*p = g for p;  S is a diagonal matrix
            double x = S.getQuick(j, j);
            if (x > 0.) {
               x = g.getQuick(j) / x;   // p[j] = g[j]/S[j]
               g.setQuick(j, x);        // overwrite g by p
            } else
               g.setQuick(j, 0.);
         }
         DoubleMatrix1D beta = alg.mult(V, g);   // V*p
         return beta.toArray();
      }
}

   /**
    * Computes the regression coefficients using the least squares method.
    * This is a simple linear regression with 2 regression coefficients,
    * @f$\alpha@f$ and @f$\beta@f$. The model is
    * @f[
    *   y = \alpha+ \beta x.
    * @f]
    * Given the @f$n@f$ data points @f$(X_i, Y_i)@f$, @f$i=0,1,…,(n-1)@f$,
    * the method computes and returns the array @f$[\alpha, \beta]@f$.
    *  @param X            the regressor variables
    *  @param Y            the response
    *  @return the regression coefficients
    */
   public static double[] calcCoefficients (double[] X, double[] Y) {
      if (X.length != Y.length)
         throw new IllegalArgumentException ("Lengths of X and Y are not equal");
      final int n = X.length;
      double[][] Xa = new double[n][1];
      for (int i = 0; i < n; i++)
         Xa[i][0] = X[i];

      return calcCoefficients0 (Xa, Y);
   }

   /**
    * Computes the regression coefficients using the least squares method.
    * This is a linear regression with a polynomial of degree `deg` @f$=
    * k@f$ and @f$k+1@f$ regression coefficients @f$\beta_j@f$. The model
    * is
    * @f[
    *   y = \beta_0 + \sum_{j=1}^k \beta_j x^j.
    * @f]
    * Given the @f$n@f$ data points @f$(X_i, Y_i)@f$, @f$i=0,1,…,(n-1)@f$,
    * the method computes and returns the array @f$[\beta_0, \beta_1, …,
    * \beta_k]@f$. Restriction: @f$n > k@f$.
    *  @param X     the regressor variables
    *  @param Y     the response
    *  @param deg   degree of the function    
    *  @return the regression coefficients
    */
   public static double[] calcCoefficients (double[] X, double[] Y, int deg) {
      final int n = X.length;
      if (n != Y.length)
         throw new IllegalArgumentException ("Lengths of X and Y are not equal");
      if (n < deg + 1)
         throw new IllegalArgumentException ("Not enough points");

      final double[] xSums = new double[2 * deg + 1];
      final double[] xySums = new double[deg + 1];
      xSums[0] = n;
      for (int i = 0; i < n; i++) {
         double xv = X[i];
         xySums[0] += Y[i];
         for (int j = 1; j <= 2 * deg; j++) {
            xSums[j] += xv;
            if (j <= deg)
               xySums[j] += xv * Y[i];
            xv *= X[i];
         }
      }
      final DoubleMatrix2D A = new DenseDoubleMatrix2D (deg + 1, deg + 1);
      final DoubleMatrix2D B = new DenseDoubleMatrix2D (deg + 1, 1);
      for (int i = 0; i <= deg; i++) {
         for (int j = 0; j <= deg; j++) {
            final int d = i + j;
            A.setQuick (i, j, xSums[d]);
         }
         B.setQuick (i, 0, xySums[i]);
      }

      return solution(A, B, deg + 1);
   }

   /**
    * Computes the regression coefficients using the least squares method.
    * This is a model for multiple linear regression. There are @f$k+1@f$
    * regression coefficients @f$\beta_j@f$, and @f$k@f$ regressors
    * variables @f$x_j@f$. The model is
    * @f[
    *   y = \beta_0 + \sum_{j=1}^k \beta_j x_j.
    * @f]
    * There are @f$n@f$ data points @f$Y_i@f$, @f$X_{ij}@f$,
    * @f$i=0,1,…,(n-1)@f$, and each @f$X_i@f$ is a @f$k@f$-dimensional
    * point. Given the response `Y[i]` and the regressor variables
    * `X[i][j]`, @f$\mathtt{i} =0,1,…,(n-1)@f$, @f$\mathtt{j}
    * =0,1,…,(k-1)@f$, the method computes and returns the array
    * @f$[\beta_0, \beta_1, …, \beta_k]@f$. Restriction: @f$n > k+1@f$.
    *  @param X            the regressor variables
    *  @param Y            the response
    *  @return the regression coefficients
    */
   public static double[] calcCoefficients0 (double[][] X, double[] Y) {
      if (X.length != Y.length)
         throw new IllegalArgumentException ("Lengths of X and Y are not equal");
      if (Y.length <= X[0].length + 1)
         throw new IllegalArgumentException ("Not enough points");

      final int n = Y.length;
      final int k = X[0].length;

      DoubleMatrix2D Xa = new DenseDoubleMatrix2D(n, k+1);
      DoubleMatrix2D Ya = new DenseDoubleMatrix2D(n, 1);

      for (int i = 0; i < n; i++) {
         Xa.setQuick (i, 0, 1.);
         for (int j = 1; j <= k; j++) {
            Xa.setQuick (i, j, X[i][j-1]);
         }
         Ya.setQuick (i, 0, Y[i]);
      }

      return solution(Xa, Ya, k + 1);
   }

   /**
    * Computes the regression coefficients using the least squares method.
    * This is a model for multiple linear regression. There are @f$k@f$
    * regression coefficients @f$\beta_j@f$, @f$j=0,1,…,(k-1)@f$ and
    * @f$k@f$ regressors variables @f$x_j@f$. The model is
    * @f[
    *   y = \sum_{j=0}^{k-1} \beta_j x_j.
    * @f]
    * There are @f$n@f$ data points @f$Y_i@f$, @f$X_{ij}@f$,
    * @f$i=0,1,…,(n-1)@f$, and each @f$X_i@f$ is a @f$k@f$-dimensional
    * point. Given the response `Y[i]` and the regressor variables
    * `X[i][j]`, @f$\mathtt{i} =0,1,…,(n-1)@f$, @f$\mathtt{j}
    * =0,1,…,(k-1)@f$, the method computes and returns the array
    * @f$[\beta_0, \beta_1, …, \beta_{k-1}]@f$. Restriction: @f$n >
    * k@f$.
    *  @param X            the regressor variables
    *  @param Y            the response
    *  @return the regression coefficients
    */
   public static double[] calcCoefficients (double[][] X, double[] Y) {
      if (X.length != Y.length)
         throw new IllegalArgumentException ("Lengths of X and Y are not equal");
      if (Y.length <= X[0].length + 1)
         throw new IllegalArgumentException ("Not enough points");

      final int n = Y.length;
      final int k = X[0].length;

      DoubleMatrix2D Xa = new DenseDoubleMatrix2D(n, k);
      DoubleMatrix2D Ya = new DenseDoubleMatrix2D(n, 1);

      for (int i = 0; i < n; i++) {
         for (int j = 0; j < k; j++) {
            Xa.setQuick (i, j, X[i][j]);
         }
         Ya.setQuick (i, 0, Y[i]);
      }

      return solution(Xa, Ya, k);
   }

}