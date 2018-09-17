/*
 * Class:        DMatrix
 * Description:  Methods for matrix calculations with double numbers
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
package umontreal.ssj.util;
   import cern.colt.matrix.*;
   import cern.colt.matrix.impl.*;
   import cern.colt.matrix.linalg.*;
   import cern.jet.math.Functions;

/**
 * This class implements a few methods for matrix calculations with `double`
 * numbers.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class DMatrix {
   private double [][] mat;        // matrix of double's
   private int r, c;               // number of rows, columns

   // [9/9] Pade numerator coefficients for exp(x);
   private static double[] cPade = {17643225600d, 8821612800d, 2075673600d,
       302702400, 30270240,  2162160, 110880, 3960, 90, 1};

   // [7/7] Pade numerator coefficients for (exp(x) - 1) / x
   private static double[] p_em1 = {1.0, 1.0/30, 1.0/30, 1.0/936,
       1.0/4680, 1.0/171600, 1.0/3603600, 1.0/259459200};

   // [7/7] Pade denominator coefficients for (exp(x) - 1) / x
   private static double[] q_em1 ={1.0, -(7.0/15), 1.0/10, -(1.0/78),
       1.0/936, -(1.0/17160), 1.0/514800, -(1.0/32432400)};


   //======================================================================

   /*
    Matrix multiplication $C = AB$. All three matrices are square, banded,
    and upper triangular. $A$ has a non-zero diagonal, \texttt{sa} non-zero
    superdiagonals, and thus a bandwidth of \texttt{sa + 1}. The non-zero
    elements of $A_{ij}$ are those for which $j - s_a \le i \le j$.
    Similarly for $B$ which has a bandwidth of \texttt{sb + 1}.
    The resulting matrix $C$ has \texttt{sa + sb} non-zero superdiagonals,
    and a bandwidth of \texttt{sa + sb + 1}.
   */
   private static void innerMultBand (DoubleMatrix2D A0, int sa,
                                      DoubleMatrix2D B0, int sb,
                                      DoubleMatrix2D C) {
      DoubleMatrix2D A, B;
      if (A0 == C)
         A = A0.copy ();
      else
         A = A0;
      if (B0 == C)
         B = B0.copy ();
      else
         B = B0;
      C.assign(0.);
      final int n = A.rows();
      int kmin, kmax;
      double x, y, z;
      for (int i = 0; i < n; ++i) {
         int jmax = Math.min(i + sa + sb, n - 1);
         for (int j = i; j <= jmax; ++j) {
            kmin = Math.max(i, j - sb);
            kmax = Math.min(i + sa, j);
            z = 0;
            for (int k = kmin; k <= kmax; ++k) {
               x = A.getQuick (i, k);
               y = B.getQuick (k, j);
               z += x * y;
            }
            C.setQuick (i, j, z);
         }
      }
   }


   //======================================================================

   private static int getScale (final DoubleMatrix2D A, double theta)
   {
      // assumes A is an upper bidiagonal matrix
      final double norm = norm1bidiag(A) / theta;
      int s;
      if (norm > 1)
         s = (int) Math.ceil(Num.log2(norm));
      else
         s = 0;
      return s;
   }


   //======================================================================

   private static DoubleMatrix2D m_taylor (final DoubleMatrix2D A)
   {
      // Compute and returns (e^A - I), using the Taylor series around 0

      final double EPS = 1.0e-12;
      final int k = A.rows();
      final int JMAX = 2 * k + 100;
      DoubleMatrix2D Sum = A.copy();
      DoubleMatrix2D Term = A.copy();

      Functions F = Functions.functions;    // alias F
      Algebra alge = new Algebra();
      double normS, normT;

      for (int j = 2; j <= JMAX; ++j) {
         Term = alge.mult(A, Term);            // Term <-- A*Term
         Term.assign (F.mult(1.0 / j));        // Term <-- Term/j
         Sum.assign (Term, F.plus);            // Sum <-- Sum + Term
         if (j > k + 5) {
            normS = alge.normInfinity(Sum);
            normT = alge.normInfinity(Term);
            if (normT <= normS * EPS)
               break;
         }
      }
      return Sum;
   }

   //======================================================================

   private static DoubleMatrix1D m_taylor (final DoubleMatrix2D A, final DoubleMatrix1D b)
   {
      // Compute and returns (exp(A) - I)b, using the Taylor series around 0

      final double EPS = 1.0e-12;
      final int k = A.rows();
      final int JMAX = 2 * k + 100;
      DoubleFactory1D factory = DoubleFactory1D.dense;
      DoubleMatrix1D Term = b.copy();
      DoubleMatrix1D Sum = factory.make(k);

      Functions F = Functions.functions;    // alias F
      Algebra alge = new Algebra();
      double Snorm, Tnorm;

      for (int j = 1; j <= JMAX; ++j) {
         Term = alge.mult(A, Term);            // Term <-- A*Term
         Term.assign (F.mult(1.0 / j));        // Term <-- Term/j
         Sum.assign (Term, F.plus);            // Sum  <-- Sum + Term
         if (j > k + 5) {
            Tnorm = alge.norm1(Term);
            Snorm = alge.norm1(Sum);
            if (Tnorm <= Snorm * EPS)
               break;
         }
      }
      return Sum;
   }


   //======================================================================

   private static DoubleMatrix2D m_expmiBidiag (final DoubleMatrix2D A)
   {
      // Use the diagonal Pade approximant of order [7/7] for (exp(A) - I)/A:
      // See Higham J.H., Functions of matrices, SIAM, 2008, p. 262.
      // This method scale A to a matrix with small norm B = A/2^s,
      // compute U = (exp(B) - I) with Pade approximants, and returns U.
      // Returns also the scale s in scale[0].

      final int n = A.rows();
      DoubleMatrix2D B = A.copy();
      DoubleFactory2D fac = DoubleFactory2D.dense;
      final DoubleMatrix2D I = fac.identity(n);
      final DoubleMatrix2D B2 = fac.make(n, n);
      final DoubleMatrix2D B4 = fac.make(n, n);
      final DoubleMatrix2D B6 = fac.make(n, n);
      DMatrix.multBand(B, 1, B, 1, B2); // B^2
      DMatrix.multBand(B2, 2, B2, 2, B4); // B^4
      DMatrix.multBand(B4, 4, B2, 2, B6); // B^6

      DoubleMatrix2D V = B6.copy();
      DoubleMatrix2D U = B4.copy();
      DMatrix.addMultBand(p_em1[4], U, 4, p_em1[2], B2, 2); // U <-- U*p_em1[4] + B2*p_em1[2]
      DMatrix.addMultBand(p_em1[6], V, 6, p_em1[0], I, 0); // V <-- V*p_em1[6] + I*p_em1[0]
      DMatrix.addMultBand(1.0, V, 6, 1.0, U, 6); // V <-- V + U

      DoubleMatrix2D W = B6.copy();
      U = B4.copy();
      DMatrix.addMultBand(p_em1[5], U, 4, p_em1[3], B2, 2); // U <-- U*p_em1[5] + B2*p_em1[3]
      DMatrix.addMultBand(p_em1[7], W, 6, p_em1[1], I, 0); // W <-- W*p_em1[7] + I*p_em1[1]
      DMatrix.addMultBand(1.0, W, 6, 1.0, U, 6); // W <-- W + U
      DMatrix.multBand(W, 6, B, 1, U);   // U <-- W*B

      DMatrix.addMultBand(1.0, V, 6, 1.0, U, 7); // V <-- V + U
      DoubleMatrix2D N = V.copy();   // numerator Pade

      V = B6.copy();
      U = B4.copy();
      DMatrix.addMultBand(q_em1[4], U, 4, q_em1[2], B2, 2); // U <-- U*q_em1[4] + B2*q_em1[2]
      DMatrix.addMultBand(q_em1[6], V, 6, q_em1[0], I, 0); // V <-- V*q_em1[6] + I*q_em1[0]
      DMatrix.addMultBand(1.0, V, 6, 1.0, U, 6); // V <-- V + U

      W = B6.copy();
      U = B4.copy();
      DMatrix.addMultBand(q_em1[5], U, 4, q_em1[3], B2, 2); // U <-- U*q_em1[5] + B2*q_em1[3]
      DMatrix.addMultBand(q_em1[7], W, 6, q_em1[1], I, 0); // W <-- W*q_em1[7] + I*q_em1[1]
      DMatrix.addMultBand(1.0, W, 6, 1.0, U, 6); // W <-- W + U
      DMatrix.multBand(W, 6, B, 1, U);   // U <-- W*B

      DMatrix.addMultBand(1.0, V, 6, 1.0, U, 7); // V <-- V + U, denominator Pade

      // Compute Pade approximant W = N/V for (exp(B) - I)/B
      DMatrix.solveTriangular(V, N, W);
      DMatrix.multBand(B, 1, W, n - 1, U);   // (exp(B) - I) = U <-- B*W

      //  calcDiagm1 (B, U);
      return U;
   }


   static void addMultTriang (final DoubleMatrix2D A, DoubleMatrix1D b, double h) {
      /* Multiplies the upper triangular matrix A by vector b multiplied by h.
         Put the result back in b.
      */
      final int n = A.rows();
      double z;
      for (int i = 0; i < n; ++i) {
         for (int j = i; j < n; ++j) {
            z = A.getQuick (i, j) * b.getQuick (j);
            b.setQuick (i, h*z);
         }
      }
   }


   //======================================================================
   /*
    * Compute the 1-norm for matrix B, which is bidiagonal. The only non-zero
    * elements are on the diagonal and the first superdiagonal.
    *
    * @param B matrix
    * @return the norm
    */
   private static double norm1bidiag (DoubleMatrix2D B)
   {
      final int n = B.rows();
      double x;
      double norm = Math.abs(B.getQuick(0, 0));
      for (int i = 1; i < n; ++i) {
         x = Math.abs(B.getQuick(i - 1, i)) + Math.abs(B.getQuick(i, i));
         if (x > norm)
            norm = x;
      }
      return norm;
   }

   /**
    * Creates a new `DMatrix` with `r` rows and `c` columns.
    *  @param r            the number of rows
    *  @param c            the number of columns
    */
   public DMatrix (int r, int c) {
      mat = new double[r][c];
      this.r = r;
      this.c = c;
   }

   /**
    * Creates a new `DMatrix` with `r` rows and `c` columns using the data
    * in `data`.
    *  @param data         the data of the new `DMatrix`
    *  @param r            the number of rows
    *  @param c            the number of columns
    */
   public DMatrix (double[][] data, int r, int c) {
      this (r, c);
      for(int i = 0; i < r; i++)
         for(int j = 0; j < c; j++)
            mat[i][j] = data[i][j];
   }

   /**
    * Copy constructor.
    *  @param that         the `DMatrix` to copy
    */
   public DMatrix (DMatrix that) {
      this (that.mat, that.r, that.c);
   }

   /**
    * Given a symmetric positive-definite matrix @f$M@f$, performs the
    * Cholesky decomposition of @f$M@f$ and returns the result as a lower
    * triangular matrix @f$L@f$, such that @f$M = L L^T@f$.
    *  @param M            the input matrix
    *  @param L            the Cholesky lower triangular matrix
    */
   public static void CholeskyDecompose (double[][] M, double[][] L) {
      int d = M.length;
      DoubleMatrix2D MM = new DenseDoubleMatrix2D (M);
      DoubleMatrix2D LL = new DenseDoubleMatrix2D (d, d);
      CholeskyDecomposition decomp = new CholeskyDecomposition (MM);
      LL = decomp.getL();
      for(int i = 0; i < L.length; i++)
         for(int j = 0; j <= i; j++)
            L[i][j] = LL.get(i,j);
      for(int i = 0; i < L.length; i++)
         for(int j = i + 1; j < L.length; j++)
            L[i][j] = 0.0;
    }

   /**
    * Given a symmetric positive-definite matrix @f$M@f$, performs the
    * Cholesky decomposition of @f$M@f$ and returns the result as a lower
    * triangular matrix @f$L@f$, such that @f$M = L L^T@f$.
    *  @param M            the input matrix
    *  @return the Cholesky lower triangular matrix
    */
   public static DoubleMatrix2D CholeskyDecompose (DoubleMatrix2D M) {
      CholeskyDecomposition decomp = new CholeskyDecomposition (M);
      return decomp.getL();
    }

   /**
    * Computes the principal components decomposition @f$M@f$ =
    * @f$U\Lambda U^{\mathsf{t}}@f$ by using the singular value
    * decomposition of matrix @f$M@f$. Puts the eigenvalues, which are the
    * diagonal elements of matrix @f$\Lambda@f$, sorted by decreasing
    * size, in vector `lambda`, and puts matrix @f$A = U\sqrt{\Lambda}@f$
    * in `A`.
    *  @param M            input matrix
    *  @param A            matrix square root of M
    *  @param lambda       the eigenvalues
    */
   public static void PCADecompose (double[][] M, double[][] A,
                                    double[] lambda) {
      int d = M.length;
      DoubleMatrix2D MM = new DenseDoubleMatrix2D (M);
      DoubleMatrix2D AA = new DenseDoubleMatrix2D (d, d);
      AA = PCADecompose(MM, lambda);

      for(int i = 0; i < d; i++)
         for(int j = 0; j < d; j++)
            A[i][j] = AA.get(i,j);
    }

   /**
    * Computes the principal components decomposition @f$M@f$ =
    * @f$U\Lambda U^{\mathsf{t}}@f$ by using the singular value
    * decomposition of matrix @f$M@f$. Puts the eigenvalues, which are the
    * diagonal elements of matrix @f$\Lambda@f$, sorted by decreasing
    * size, in vector `lambda`. Returns matrix @f$A = U\sqrt{\Lambda}@f$.
    *  @param M            input matrix
    *  @param lambda       the eigenvalues
    *  @return matrix square root of M
    */
   public static DoubleMatrix2D PCADecompose (DoubleMatrix2D M,
                                              double[] lambda) {
      // L'objet SingularValueDecomposition permet de recuperer la matrice
      // des valeurs propres en ordre decroissant et celle des vecteurs propres de
      // sigma (pour une matrice symetrique et definie-positive seulement).

      SingularValueDecomposition sv = new SingularValueDecomposition(M);
      // D contient les valeurs propres sur la diagonale
      DoubleMatrix2D D = sv.getS ();

      for (int i = 0; i < D.rows(); ++i)
         lambda[i] = D.getQuick (i, i);

      // Calculer la racine carree des valeurs propres
      for (int i = 0; i < D.rows(); ++i)
         D.setQuick (i, i, Math.sqrt (lambda[i]));
      DoubleMatrix2D P = sv.getV();
      // Multiplier par la matrice orthogonale (ici P)
      return P.zMult (D, null);
   }

   /**
    * Solves the matrix equation @f$Ax = b@f$ using LU decomposition.
    * @f$A@f$ is a square matrix, @f$b@f$ and @f$x@f$ are vectors. Returns
    * the solution @f$x@f$.
    *  @param A            square matrix
    *  @param b            right side vector
    *  @return the solution vector
    */
   public static double[] solveLU (double[][] A, double[] b) {
      DoubleMatrix2D M = new DenseDoubleMatrix2D(A);
      DoubleMatrix1D c = new DenseDoubleMatrix1D(b);
      LUDecompositionQuick lu = new LUDecompositionQuick();
      lu.decompose(M);
      lu.solve(c);
      return c.toArray();
   }

   /**
    * Solve the triangular matrix equation @f$UX = B@f$ for @f$X@f$.
    * @f$U@f$ is a square upper triangular matrix. @f$B@f$ and @f$X@f$
    * must have the same number of columns.
    *  @param U            input matrix
    *  @param B            right-hand side matrix
    *  @param X            output matrix
    */
   public static void solveTriangular (DoubleMatrix2D U, DoubleMatrix2D B,
                                       DoubleMatrix2D X) {
      final int n = U.rows();
      final int m = B.columns();
      double y, z;
      X.assign(0.);
      for (int j = 0; j < m; ++j) {
         for (int i = n - 1; i >= 0; --i) {
            z = B.getQuick(i, j);
            for (int k = i + 1; k < n; ++k)
               z -= U.getQuick(i, k) * X.getQuick(k, j);
            z /= U.getQuick(i, i);
            X.setQuick(i, j, z);
         }
      }
   }

   /**
    * Similar to \ref #exp(final DoubleMatrix2D).
    *  @param A            input matrix
    *  @return the exponential of @f$A@f$
    */
   public static double[][] exp (double[][] A) {
      DoubleMatrix2D V = new DenseDoubleMatrix2D(A);
      DoubleMatrix2D R = exp(V);
      return R.toArray();
   }

   /**
    * Returns @f$e^A@f$, the exponential of the square matrix @f$A@f$. The
    * scaling and squaring method @cite mHIG09a is used with Padé
    * approximants to compute the exponential.
    *  @param A            input matrix
    *  @return the exponential of @f$A@f$
    */
   public static DoubleMatrix2D exp (final DoubleMatrix2D A) {
      /*
       * Use the diagonal Pade approximant of order [9/9] for exp:
       * See Higham J.H., Functions of matrices, SIAM, 2008.
       */
      DoubleMatrix2D B = A.copy();
      int n = B.rows();
      Algebra alge = new Algebra();
      final double mu = alge.trace(B) / n;
      double x;

      // B <-- B - mu*I
      for (int i = 0; i < n; ++i) {
         x = B.getQuick (i, i);
         B.setQuick (i, i, x - mu);
      }
      /*
      int bal = 0;
      if (bal > 0) {
         throw new UnsupportedOperationException ("   balancing");
      } */

      final double THETA9 = 2.097847961257068;   // in Higham
      int s = getScale (B, THETA9);

      Functions F = Functions.functions;    // alias F
      // B <-- B/2^s
      double v = 1.0 / Math.pow(2.0, s);
      if (v <= 0)
          throw new IllegalArgumentException ("   v <= 0");
      B.assign (F.mult(v));

      DoubleFactory2D fac = DoubleFactory2D.dense;
      final DoubleMatrix2D B0 = fac.identity(n);    // B^0 = I
      final DoubleMatrix2D B2 = alge.mult(B, B);    // B^2
      final DoubleMatrix2D B4 = alge.mult(B2, B2);  // B^4

      DoubleMatrix2D T = B2.copy();          // T = work matrix
      DoubleMatrix2D W = B4.copy();          // W = work matrix
      W.assign (F.mult(cPade[9]));           // W <-- W*cPade[9]
      W.assign (T, F.plusMult(cPade[7]));    // W <-- W + T*cPade[7]
      DoubleMatrix2D U = alge.mult(B4, W);   // U <-- B4*W

      // T = B2.copy();
      W = B4.copy();
      W.assign (F.mult(cPade[5]));           // W <-- W*cPade[5]
      W.assign (T, F.plusMult(cPade[3]));    // W <-- W + T*cPade[3]
      W.assign (B0, F.plusMult(cPade[1]));   // W <-- W + B0*cPade[1]
      U.assign (W, F.plus);                  // U <-- U + W
      U = alge.mult(B, U);                   // U <-- B*U

      // T = B2.copy();
      W = B4.copy();
      W.assign (F.mult(cPade[8]));           // W <-- W*cPade[8]
      W.assign (T, F. plusMult(cPade[6]));   // W <-- W + T*cPade[6]
      DoubleMatrix2D V = alge.mult(B4, W);   // V <-- B4*W

      // T = B2.copy();
      W = B4.copy();
      W.assign (F.mult(cPade[4]));           // W <-- W*cPade[4]
      W.assign (T, F.plusMult(cPade[2]));    // W <-- W + T*cPade[2]
      W.assign (B0, F.plusMult(cPade[0]));   // W <-- W + B0*cPade[0]
      V.assign (W, F.plus);                  // V <-- V + W

      W = V.copy();
      W.assign(U, F.plus);                   // W = V + U, Pade numerator
      T = V.copy();
      T.assign(U, F.minus);                  // T = V - U, Pade denominator

      // Compute Pade approximant for exponential = W / T
      LUDecomposition lu = new LUDecomposition(T);
      B = lu.solve(W);

      if (false) {
         // This overflows for large |mu|
         // B <-- B^(2^s)
         for(int i = 0; i < s; i++)
            B = alge.mult(B, B);
         /*
         if (bal > 0) {
            throw new UnsupportedOperationException ("   balancing");
         } */
         v = Math.exp(mu);
         B.assign (F.mult(v));               // B <-- B*e^mu

      } else {
         // equivalent to B^(2^s) * e^mu, but only if no balancing
         double r = mu * v;
         r = Math.exp(r);
         B.assign (F.mult(r));
         for (int i = 0; i < s; i++)
            B = alge.mult(B, B);
      }

      return B;
   }

   /**
    * Returns @f$e^A@f$, the exponential of the *bidiagonal* square matrix
    * @f$A@f$. The only non-zero elements of @f$A@f$ are on the diagonal
    * and on the first superdiagonal. This method is faster than
    * \ref #exp(final DoubleMatrix2D) because of the special form of
    * @f$A@f$.
    *  @param A            bidiagonal matrix
    *  @return @f$e^A@f$
    */
   public static DoubleMatrix2D expBidiagonal (final DoubleMatrix2D A) {
      // Use the diagonal Pade approximant of order [9/9] for exp:
      // See Higham J.H., Functions of matrices, SIAM, 2008.
      // This method scale A to a matrix with small norm B = A/2^s,
      // compute U = exp(B) with Pade approximants, and returns U.

      DoubleMatrix2D B = A.copy();
      final int n = B.rows();
      Algebra alge = new Algebra();
      final double mu = alge.trace(B) / n;
      double x;

      // B <-- B - mu*I
      for (int i = 0; i < n; ++i) {
         x = B.getQuick(i, i);
         B.setQuick(i, i, x - mu);
      }

      final double THETA9 = 2.097847961257068; // in Higham
      int s = getScale (B, THETA9);
      final double v = 1.0 / Math.pow(2.0, s);
      if (v <= 0)
         throw new IllegalArgumentException("   v <= 0");
      DMatrix.multBand(B, 1, v); // B <-- B/2^s

      DoubleFactory2D fac = DoubleFactory2D.dense;
      DoubleMatrix2D T = fac.make(n, n);
      DoubleMatrix2D B4 = fac.make(n, n);
      DMatrix.multBand(B, 1, B, 1, T); // B^2
      DMatrix.multBand(T, 2, T, 2, B4); // B^4

      DoubleMatrix2D W = B4.copy(); // W = work matrix
      DMatrix.addMultBand(cPade[9], W, 4, cPade[7], T, 2); // W <-- W*cPade[9] + T*cPade[7]
      DoubleMatrix2D U = fac.make(n, n);
      DMatrix.multBand(W, 4, B4, 4, U); // U <-- B4*W

      W = B4.copy();
      DMatrix.addMultBand(cPade[5], W, 4, cPade[3], T, 2); // W <-- W*cPade[5] + T*cPade[3]
      for (int i = 0; i < n; ++i) {   // W <-- W + I*cPade[1]
         x = W.getQuick(i, i);
         W.setQuick(i, i, x + cPade[1]);
      }
      DMatrix.addMultBand(1.0, U, 8, 1.0, W, 4); // U <-- U + W
      DMatrix.multBand(B, 1, U, 8, U); // U <-- B*U

      W = B4.copy();
      DMatrix.addMultBand(cPade[8], W, 4, cPade[6], T, 2); // W <-- W*cPade[8] + T*cPade[6]
      DoubleMatrix2D V = B;
      DMatrix.multBand(W, 4, B4, 4, V); // V <-- B4*W

      W = B4.copy();
      DMatrix.addMultBand(cPade[4], W, 4, cPade[2], T, 2); // W <-- W*cPade[4] + T*cPade[2]
      for (int i = 0; i < n; ++i) {   // W <-- W + I*cPade[0]
         x = W.getQuick(i, i);
         W.setQuick(i, i, x + cPade[0]);
      }
      DMatrix.addMultBand(1.0, V, 8, 1.0, W, 4); // V <-- V + W

      W = V.copy();
      DMatrix.addMultBand(1.0, W, 9, 1.0, U, 9); // W = V + U, Pade numerator
      T = V.copy();
      DMatrix.addMultBand(1.0, T, 9, -1.0, U, 9); // T = V - U, Pade denominator

      // Compute Pade approximant B = W/T for exponential
      DMatrix.solveTriangular(T, W, B);

      // equivalent to B^(2^s) * e^mu
      double r = mu * v;
      r = Math.exp(r);
      DMatrix.multBand(B, n - 1, r); // B <-- B*r

      T.assign(0.);

      for (int i = 0; i < s; i++) {
         DMatrix.multBand(B, n - 1, B, n - 1, T);
         B = T.copy();
      }

      return B;
	}

   /**
    * Computes @f$c = e^A b@f$, where @f$e^A@f$ is the exponential of the
    * *bidiagonal* square matrix @f$A@f$. The only non-zero elements of
    * @f$A@f$ are on the diagonal and on the first superdiagonal. Uses the
    * scaling and squaring method @cite mHIG09a with Padé
    * approximants. Returns @f$c@f$.
    *  @param A            bidiagonal matrix
    *  @param b            vector
    *  @return @f$c@f$
    */
   public static DoubleMatrix1D expBidiagonal (final DoubleMatrix2D A,
                                               final DoubleMatrix1D b) {
      // This is probably not efficient;
      DoubleMatrix2D U = expBidiagonal (A);   // U = exp(A)
      Algebra alge = new Algebra();
      return alge.mult(U, b);
   }

   /**
    * Computes @f$e^A - I@f$, where @f$e^A@f$ is the exponential of the
    * *bidiagonal* square matrix @f$A@f$. The only non-zero elements of
    * @f$A@f$ are on the diagonal and on the first superdiagonal. Uses the
    * scaling and squaring method @cite mSKA09a, @cite mHIG09a&thinsp;
    * with Padé approximants. Returns @f$e^A - I@f$.
    *  @param A            bidiagonal matrix
    *  @return @f$(e^A - I)b@f$
    */
   public static DoubleMatrix2D expmiBidiagonal (final DoubleMatrix2D A) {
      // Use the diagonal Pade approximant of order [7/7] for (exp(A) - I)/A:
      // See Higham J.H., Functions of matrices, SIAM, 2008, p. 262.

      DoubleMatrix2D B = A.copy();
      final double THETA = 1.13; // theta_{7,1}
      int s = getScale (B, THETA);
      final double v = 1.0 / Math.pow(2.0, s);
      if (v <= 0)
         throw new IllegalArgumentException("   v <= 0");
      DMatrix.multBand(B, 1, v);   // B <-- B/2^s
      DoubleMatrix2D U = m_expmiBidiag (B);   // U = exp(B) - I

      DoubleMatrix2D N = U.copy();
      addIdentity (N);                     // N <-- exp(B)
      DoubleMatrix2D V = N.copy();         // V <-- exp(B)

      // Now undo scaling of B = A/2^s using
      // (exp(A) - I) = (exp(B) - I)(exp(B) + I)(exp(B^2) + I) ... (exp(B^(2^(s-1))) + I)

      Algebra alge = new Algebra();
      for (int i = 1; i <= s; i++) {
         addIdentity (N);        // N <-- exp(B) + I
         U = alge.mult(N, U);    // U <-- N*U
         if (i < s) {
            V = alge.mult(V, V);    // V <-- V*V
            N = V.copy();
         }
      }

      return U;
      }

   /**
    * Computes @f$c = (e^A - I)b@f$, where @f$e^A@f$ is the exponential of
    * the *bidiagonal* square matrix @f$A@f$. The only non-zero elements
    * of @f$A@f$ are on the diagonal and on the first superdiagonal. Uses
    * the scaling and squaring method @cite mSKA09a, @cite mHIG09a&thinsp;
    * with a Taylor expansion. Returns @f$c@f$.
    *  @param A            bidiagonal matrix
    *  @param b            vector
    *  @return @f$c@f$
    */
   public static DoubleMatrix1D expmiBidiagonal (final DoubleMatrix2D A,
                                                 final DoubleMatrix1D b) {
      DoubleMatrix2D F = A.copy();
      int s = getScale (F, 1.0 / 16.0);
      final double v = 1.0 / Math.pow(2.0, s);
      if (v <= 0)
         throw new IllegalArgumentException("   v <= 0");
      DMatrix.multBand(F, 1, v);   // F <-- F/2^s

      DoubleMatrix2D U = expBidiagonal (F);   // U = exp(F)
      DoubleMatrix2D N = U.copy();
      DoubleMatrix1D C = m_taylor (F, b);   // C = (exp(F) - I)b
      DoubleMatrix1D D = C.copy();
      Algebra alge = new Algebra();

      // Now undo scaling of F = A/2^s using
      //   (exp(A) - I)b = (exp(F^(2^(s-1))) + I)...(exp(F^2) + I)(exp(F) + I)(exp(F) - I)b
      for (int i = 1; i <= s; i++) {
         addIdentity (N);        // N <-- exp(F) + I
         C = alge.mult(N, C);    // C <-- N*C
         if (i < s) {
            U = alge.mult(U, U);    // U <-- U*U
            N = U.copy();
         }
      }

      return C;
   }
private static void addIdentity (DoubleMatrix2D A) {
      // add identity matrix to matrix A:  A <-- A + I
      final int n = A.rows();
      double x;
      for (int i = 0; i < n; ++i) {
         x = A.getQuick(i, i);
         A.setQuick(i, i, x + 1.0);
      }
   }


   private static void calcDiagm1 (DoubleMatrix2D A, DoubleMatrix2D R) {
      // calc diagonal of expm1 of triangular matrix A:  exp(A) - I
      final int n = A.rows();
      double x, v;
      for (int i = 0; i < n; ++i) {
         x = A.getQuick(i, i);
         v = Math.expm1(x);      // exp(x) - 1
         R.setQuick(i, i, v);
      }
   }

/**
 * Matrix multiplication @f$C = AB@f$. All three matrices are square, banded,
 * and upper triangular. @f$A@f$ has a non-zero diagonal, `sa` non-zero
 * superdiagonals, and thus a bandwidth of `sa + 1`. The non-zero elements of
 * @f$A_{ij}@f$ are those for which @f$j - s_a \le i \le j@f$. Similarly for
 * @f$B@f$ which has a bandwidth of `sb + 1`. The resulting matrix @f$C@f$
 * has `sa + sb` non-zero superdiagonals, and a bandwidth of `sa + sb + 1`.
 *  @param A            input matrix
 *  @param sa           number of superdiagonals of A
 *  @param B            input matrix
 *  @param sb           number of superdiagonals of B
 *  @param C            result
 */
static void multBand (final DoubleMatrix2D A, int sa,
                         final DoubleMatrix2D B, int sb,
                         DoubleMatrix2D C) {
      innerMultBand (A, sa, B, sb, C);
   }

   /**
    * Multiplication of the matrix @f$A@f$ by the scalar @f$h@f$. @f$A@f$
    * is a square banded upper triangular matrix. It has a non-zero
    * diagonal, `sa` superdiagonals, and thus a bandwidth of `sa + 1`. The
    * result of the multiplication is put back in @f$A@f$.
    *  @param A            input and output matrix
    *  @param sa           number of superdiagonals of A
    *  @param h            scalar
    */
   static void multBand (DoubleMatrix2D A, int sa, double h) {
      final int n = A.rows();
      double z;
      for (int i = 0; i < n; ++i) {
         int jmax = Math.min(i + sa, n - 1);
         for (int j = i; j <= jmax; ++j) {
            z = A.getQuick (i, j);
            A.setQuick (i, j, z*h);
         }
      }
   }

   /**
    * Addition of the matrices @f$gA + hB@f$ after multiplication with the
    * scalars @f$g@f$ and @f$h@f$. @f$A@f$ is a square banded upper
    * triangular matrix. It has a non-zero diagonal, `sa` superdiagonals,
    * and thus a bandwidth of `sa + 1`. Similarly for @f$B@f$ which has a
    * bandwidth of `sb + 1`. The result is put back in @f$A@f$.
    *  @param g            coefficient multiplying A
    *  @param A            input and output matrix
    *  @param sa           number of superdiagonals of A
    *  @param h            coefficient multiplying B
    *  @param B            input matrix
    *  @param sb           number of superdiagonals of B
    */
   static void addMultBand (double g, DoubleMatrix2D A, int sa,
                            double h, final DoubleMatrix2D B, int sb) {
      DoubleMatrix2D S;
      S = A.copy ();
      final int n = A.rows();
      double z;
      for (int i = 0; i < n; ++i) {
         int jmax = Math.max(i + sa, i + sb);
         jmax = Math.min(jmax, n - 1);
         for (int j = i; j <= jmax; ++j) {
            z = g*S.getQuick (i, j) + h*B.getQuick (i, j);
            A.setQuick (i, j, z);
         }
      }
   }

   /**
    * Copies the matrix @f$M@f$ into @f$R@f$.
    *  @param M            original matrix
    *  @param R            output matrix
    */
   public static void copy (double[][] M, double[][] R) {
       for (int i = 0; i < M.length; i++) {
         for(int j = 0; j < M[i].length; j++) {
            R[i][j] = M[i][j];
         }
      }
   }

   /**
    * Returns matrix @f$M@f$ as a string. It is displayed in matrix form,
    * with each row on a line.
    *  @return the content of @f$M@f$
    */
   public static String toString(double[][] M) {
      StringBuffer sb = new StringBuffer();

      sb.append("{" + PrintfFormat.NEWLINE);
      for (int i = 0; i < M.length; i++) {
         sb.append("   { ");
         for(int j = 0; j < M[i].length; j++) {
            sb.append(M[i][j] + " ");
            if (j < M.length - 1)
               sb.append(" ");
         }
         sb.append("}" + PrintfFormat.NEWLINE);
      }
      sb.append("}");

      return sb.toString();
   }

   /**
    * Creates a  String containing all the data of the `DMatrix`. The
    * result is displayed in matrix form, with each row on a line.
    *  @return the content of the `DMatrix`
    */
   public String toString() {
      return toString(mat);
   }

   /**
    * Returns the number of rows of the `DMatrix`.
    *  @return the number of rows
    */
   public int numRows() {
      return r;
   }

   /**
    * Returns the number of columns of the `DMatrix`.
    *  @return the number of columns
    */
   public int numColumns() {
      return c;
   }

   /**
    * Returns the matrix element in the specified row and column.
    *  @param row          the row of the selected element
    *  @param column       the column of the selected element
    *  @return the value of the element
    *
    *  @exception IndexOutOfBoundsException if the selected element would
    * be outside the `DMatrix`
    */
   public double get (int row, int column) {
      if (row >= r || column >= c)
         throw new IndexOutOfBoundsException();

      return mat[row][column];
   }

   /**
    * Sets the value of the element in the specified row and column.
    *  @param row          the row of the selected element
    *  @param column       the column of the selected element
    *  @param value        the new value of the element
    *  @exception IndexOutOfBoundsException if the selected element would
    * be outside the `DMatrix`
    */
   public void set (int row, int column, double value) {
      if (row >= r || column >= c)
         throw new IndexOutOfBoundsException();

      mat[row][column] = value;
   }

   /**
    * Returns the transposed matrix. The rows and columns are
    * interchanged.
    *  @return the transposed matrix
    */
   public DMatrix transpose() {
      DMatrix result = new DMatrix(c,r);

      for(int i = 0; i < r; i++)
         for(int j = 0; j < c; j++)
            result.mat[j][i] = mat[i][j];

      return result;
   }

}
