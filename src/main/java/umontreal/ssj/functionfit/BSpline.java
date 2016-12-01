/*
 * Class:        BSpline
 * Description:  B-spline
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
package umontreal.ssj.functionfit;

import umontreal.ssj.util.Misc;
import umontreal.ssj.util.RootFinder;
import umontreal.ssj.functions.MathFunction;
import umontreal.ssj.functions.MathFunctionWithIntegral;
import umontreal.ssj.functions.MathFunctionWithDerivative;
import umontreal.ssj.functions.MathFunctionWithFirstDerivative;
import umontreal.ssj.functions.MathFunctionUtil;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.DoubleFactory2D;
import java.util.Arrays;
import java.util.Random;
import java.io.*;

/**
 * Represents a B-spline with control points at @f$(X_i, Y_i)@f$. Let
 * @f$\mathbf{P_i}=(X_i, Y_i)@f$, for @f$i=0,…,n-1@f$, be a *control point*
 * and let @f$t_j@f$, for @f$j=0,…,m-1@f$ be a *knot*. A B-spline
 * @cite mDEB78a&thinsp; of degree @f$p=m-n-1@f$ is a parametric curve
 * defined as
 * @f[
 *   \mathbf{P(t)} = \sum_{i=0}^{n-1} N_{i, p}(t) \mathbf{P_i},\mbox{ for }t_p\le t\le t_{m-p-1}.
 * @f]
 * Here,
 * @f{align*}{
 *    N_{i, p}(t) 
 *    & 
 *   =
 *    \frac{t-t_i}{t_{i+p} - t_i}N_{i, p-1}(t) + \frac{t_{i+p+1} - t}{t_{i+p+1} - t_{i+1}}N_{i+1, p-1}(t)
 *    \\ 
 *    N_{i, 0}(t) 
 *    & 
 *   =
 *    \left\{
 *   \begin{array}{ll}
 *    1
 *    & 
 *   \mbox{ for }t_i\le t\le t_{i+1},
 *    \\ 
 *   0
 *   \mbox{ elsewhere}. 
 *   \end{array}
 *   \right.
 * @f}
 * This class provides methods to evaluate @f$\mathbf{P(t)}=(X(t), Y(t))@f$
 * at any value of @f$t@f$, for a B-spline of any degree @f$p\ge1@f$. Note
 * that the  #evaluate(double) method of this class can be slow, since it
 * uses a root finder to determine the value of @f$t^*@f$ for which
 * @f$X(t^*)=x@f$ before it computes @f$Y(t^*)@f$.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class BSpline implements MathFunction

,
MathFunctionWithIntegral, MathFunctionWithDerivative, MathFunctionWithFirstDerivative{
   // Formula taken from http://www.ibiblio.org/e-notes/Splines/Basis.htm
   // and http://en.wikipedia.org/wiki/De_Boor_algorithm
   private double[] x;     //x original
   private double[] y;     //y original

   private int degree;

   //variables sur lesquelles on travaille
   private double[] myX;
   private double[] myY;
   private double[] knots;

   /**
    * Constructs a new uniform B-spline of degree `degree` with control
    * points at (<tt>x[i]</tt>, <tt>y[i]</tt>). The knots of the resulting
    * B-spline are set uniformly from `x[0]` to `x[n-1]`.
    *  @param x            the values of @f$X@f$.
    *  @param y            the values of @f$Y@f$.
    *  @param degree       the degree of the B-spline.
    */
   public BSpline (final double[] x, final double[] y, final int degree) {
      if (x.length != y.length)
         throw new IllegalArgumentException("The arrays x and y must share the same length");
      this.degree = degree;
      this.x = x.clone();
      this.y = y.clone();
      init(x, y, null);
   }

   /**
    * Constructs a new uniform B-spline with control points at
    * (<tt>x[i]</tt>, <tt>y[i]</tt>), and knot vector given by the array
    * `knots`.
    *  @param x            the values of @f$X@f$.
    *  @param y            the values of @f$Y@f$.
    *  @param knots        the knots of the B-spline.
    */
   public BSpline (final double[] x, final double[] y, final double[] knots) {
      if (x.length != y.length)
         throw new IllegalArgumentException("The arrays x and y must share the same length");
      if (knots.length <= x.length + 1)
         throw new IllegalArgumentException("The number of knots must be at least n+2");

      this.x = x.clone();
      this.y = y.clone();
      this.knots = knots.clone();
      init(x, y, knots);
   }

   /**
    * Returns the @f$X_i@f$ coordinates for this spline.
    *  @return the @f$X_i@f$ coordinates.
    */
   public double[] getX() {
      return myX.clone ();
   }

   /**
    * Returns the @f$Y_i@f$ coordinates for this spline.
    *  @return the @f$Y_i@f$ coordinates.
    */
   public double[] getY() {
      return myY.clone ();
   }

   /**
    * Returns the knot maximal value.
    *  @return the @f$Y_i@f$ coordinates.
    */
   public double getMaxKnot() {
      return knots[knots.length-1];
   }

   /**
    * Returns the knot minimal value.
    *  @return the @f$Y_i@f$ coordinates.
    */
   public double getMinKnot() {
      return knots[0];
   }

   /**
    * Returns an array containing the knot vector @f$(t_0, t_{m-1})@f$.
    *  @return the knot vector.
    */
   public double[] getKnots() {
      return knots.clone ();
   }

   /**
    * Returns a B-spline curve of degree `degree` interpolating the
    * @f$(x_i, y_i)@f$ points @cite mDEB78a&thinsp;. This method uses the
    * uniformly spaced method for interpolating points with a B-spline
    * curve, and a uniformed clamped knot vector, as described in
    * [http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/](http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/).
    *  @param x            the values of @f$X@f$.
    *  @param y            the values of @f$Y@f$.
    *  @param degree       the degree of the B-spline.
    *  @return the B-spline curve.
    */
   public static BSpline createInterpBSpline (double[] x, double[] y,
                                              int degree) {
      if (x.length != y.length)
         throw new IllegalArgumentException("The arrays x and y must share the same length");
      if (x.length <= degree)
         throw new IllegalArgumentException("The arrays length must be greater than degree");

      int n = x.length-1;
      //compute t : parameters vector uniformly from 0 to 1
      double[] t = new double[x.length];
      for(int i =0; i<t.length; i++) {
         t[i] = (double)i/(t.length-1);
      }

      //compute U : clamped knots vector uniformly from 0 to 1
      double U[] = new double[x.length + degree + 1];
      int m = U.length-1;
      for(int i =0; i<=degree; i++)
         U[i] = 0;
      for(int i =1; i<x.length-degree; i++)
         U[i+degree] = (double)i/(x.length-degree);
      for(int i = U.length-1-degree; i<U.length; i++)
         U[i] = 1;


      //compute matrix N : made of BSpline coefficients
      double [][] N = new double[x.length][x.length];
      for(int i = 0; i<x.length; i++) {
            N[i] = computeN(U, degree, t[i], x.length);
      }

      //initialize D : initial points matrix
      double [][] D = new double[x.length][2];
      for(int i =0; i<x.length; i++) {
         D[i][0] = x[i];
         D[i][1] = y[i];
      }

      //solve the linear equation system using colt library
      DoubleMatrix2D coltN = DoubleFactory2D.dense.make(N);
      DoubleMatrix2D coltD = DoubleFactory2D.dense.make(D);
      DoubleMatrix2D coltP;
      coltP = Algebra.ZERO.solve(coltN, coltD);

      return new BSpline(coltP.viewColumn(0).toArray(), coltP.viewColumn(1).toArray(), U);
   }
   
   
   /**
    * Returns a B-spline curve of degree `degree` smoothing @f$(x_i,
    * y_i)@f$, for @f$i=0,…,n@f$ points. The precision depends on the
    * parameter @f$hp1@f$: @f$1 \le\mathtt{degree} \le hp1<n@f$, which
    * represents the number of control points used by the new B-spline
    * curve, minimizing the quadratic error
    * @f[
    *   L = \sum_{i=0}^n\left( \frac{Y_i - S_i(X_i)}{W_i}\right)^2.
    * @f]
    * This method uses the uniformly spaced method for interpolating
    * points with a B-spline curve and a uniformed clamped knot vector, as
    * described in
    * [http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/](http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/).
    *  @param x            the values of @f$X@f$.
    *  @param y            the values of @f$Y@f$.
    *  @param degree       the degree of the B-spline.
    *  @param hp1          the desired number of control points.
    *  @return the B-spline curve.
    */
   public static BSpline createApproxBSpline (double[] x, double[] y,
                                              int degree, int hp1) {
      if (x.length != y.length)
         throw new IllegalArgumentException("The arrays x and y must share the same length");
      if (x.length <= degree)
         throw new IllegalArgumentException("The arrays length must be greater than degree");

      int h = hp1 - 1;
      int n = x.length-1;
      
      //compute t : parameters vector uniformly from 0 to 1
      double[] t = new double[x.length];
      for(int i = 0; i < t.length; i++) {
         t[i] = (double)i / n;
      }

      //compute U : clamped knots vector uniformly from 0 to 1
      int m = h + degree + 1;
      double U[] = new double[m + 1];
      for(int i = 0; i <= degree; i++)
         U[i] = 0;
      for(int i = 1; i < hp1 - degree; i++)
         U[i+degree] = (double)i/(hp1 - degree);
      for(int i = m-degree; i<= m; i++)
         U[i] = 1;

      
      //compute matrix N : composed of BSpline coefficients
      double [][] N = new double[n+1][h+1];
      for(int i = 0; i < N.length; i++) {
         N[i] = computeN(U, degree, t[i], h+1);
      }

      //initialize D : initial points matrix
      double [][] D = new double[x.length][2];
      for(int i = 0; i < x.length; i++) {
         D[i][0] = x[i];
         D[i][1] = y[i];
      }

      //compute Q :
      double[][] tempQ = new double[x.length][2];
      for(int k = 1; k < n; k++) {
         tempQ[k][0] = D[k][0] - N[k][0]*D[0][0] - N[k][h]*D[D.length-1][0];
         tempQ[k][1] = D[k][1] - N[k][0]*D[0][1] - N[k][h]*D[D.length-1][1];
      }
      double[][] Q = new double[h-1][2];
      for(int i = 1; i < h; i++) {
         for(int k = 1; k < n; k++) {
            Q[i-1][0] += N[k][i]*tempQ[k][0];
            Q[i-1][1] += N[k][i]*tempQ[k][1];
         }
      }
      
      // compute N matrix for computation:
      double [][] N2 = new double[n-1][h-1];
      for(int i = 0; i < N2.length; i++) {
         for (int j = 0; j < h-1; j++)
            N2[i][j] = N[i+1][j+1];
      }

      //solve the linear equation system using colt library
      DoubleMatrix2D coltQ = DoubleFactory2D.dense.make(Q);
      DoubleMatrix2D coltN = DoubleFactory2D.dense.make(N2);
      DoubleMatrix2D coltM = Algebra.ZERO.mult(Algebra.ZERO.transpose(coltN), coltN);
      DoubleMatrix2D coltP = Algebra.ZERO.solve(coltM, coltQ);
      double[] pxTemp = coltP.viewColumn(0).toArray();
      double[] pyTemp = coltP.viewColumn(1).toArray();
      double[] px = new double[hp1];
      double[] py = new double[hp1];
      px[0] = D[0][0];
      py[0] = D[0][1];
      px[h] = D[D.length-1][0];
      py[h] = D[D.length-1][1];
      for(int i = 0; i < pxTemp.length; i++) {
         px[i+1] = pxTemp[i];
         py[i+1] = pyTemp[i];
      }

      return new BSpline(px, py, U);
      // return new BSpline(px, py, degree);
   }
   
   
   /**
    * Returns the derivative B-spline object of the current variable.
    * Using this function and the returned object, instead of the
    * `derivative` method, is strongly recommended if one wants to compute
    * many derivative values.
    *  @return the derivative B-spline of the current variable.
    */
   public BSpline derivativeBSpline() {
      double xTemp[] = new double[this.myX.length-1];
      double yTemp[] = new double[this.myY.length-1];
      for(int i = 0; i<xTemp.length; i++) {
         xTemp[i] = (myX[i+1]-myX[i])*degree/(knots[i+degree+1]-knots[i+1]);
         yTemp[i] = (myY[i+1]-myY[i])*degree/(knots[i+degree+1]-knots[i+1]);
      }

      double [] newKnots = new double[knots.length-2];
      for(int i = 0; i < newKnots.length; i++) {
         newKnots[i] = knots[i+1];
      }

      //tri pas optimise du tout
      double xTemp2[] = new double[this.myX.length-1];
      double yTemp2[] = new double[this.myY.length-1];
      for(int i = 0; i<xTemp.length; i++) {
         int k=0;
         for(int j = 0; j<xTemp.length; j++) {
            if(xTemp[i] > xTemp[j])
               k++;
         }
         while(xTemp2[k] != 0)
            k++;
         xTemp2[k] = xTemp[i];
         yTemp2[k] = yTemp[i];
      }

      return new BSpline(xTemp2, yTemp2, newKnots);
   }

   /**
    * Returns the @f$i@f$th derivative B-spline object of the current
    * variable; @f$i@f$ must be less than the degree of the original
    * B-spline. Using this function and the returned object, instead of
    * the `derivative` method, is strongly recommended if one wants to
    * compute many derivative values.
    *  @param i            the degree of the derivative.
    *  @return the ith derivative.
    */
   public BSpline derivativeBSpline (int i) {
      BSpline bs = this;
      while(i > 0) {
         i--;
         bs = bs.derivativeBSpline();
      }
      return bs;
   }

   public double evaluate(final double u) {
      final MathFunction xFunction = new MathFunction () {
         public double evaluate (double t) {
            return evalX(t) - u;
         }
      };
      // brentDekker may be unstable; using bisection instead
      // final double t = RootFinder.brentDekker (0, 1, xFunction, 1e-6);
      final double t = RootFinder.bisection (0, 1, xFunction, 1e-6);
      return evalY(t);
   }

   public double integral (double a, double b) {
      return MathFunctionUtil.simpsonIntegral (this, a, b, 500);
   }

   public double derivative(double u) {
      return derivativeBSpline().evaluate(u);
   }

   public double derivative(double u, int n) {
      return derivativeBSpline(n).evaluate(u);
   }

   private void init(double[] x, double[] y, double [] initialKnots) {
      if(initialKnots == null) {
         //Cree un vecteur de noeud uniforme entre 0 et 1
         knots = new double[x.length+degree+1];
         for(int i = degree; i < this.knots.length-degree; i++)
            this.knots[i]= (double)(i-degree)/(knots.length - (2.0*degree) -1);
         for(int i = this.knots.length-degree; i < this.knots.length; i++)
            this.knots[i]=this.knots[i-1];
         for(int i = degree; i > 0; i--)
            this.knots[i-1]=this.knots[i];

         // cree notre vecteur interne de Points de controle
         // ici, aucune modification a faire sur les tableaux originaux
         myX = x;
         myY = y;
      }
      else {
         degree = initialKnots.length - x.length -1;

      // on adapte le tableau des noeuds a notre algorithme
      // le tableau knot necessite d'avoir degree fois la meme valeur en debut et en fin de tableau
      // on adapte la taille des tableau X et Y en consequence afin de continuer a respecter la condition :
      // x.length + degree + 1 = this.knots.length
      // Cette modification n'influence pas le resultat et permet de fermer notre courbe

         //Compute the number of values that need to be added
         int iBorneInf = 1, iBorneSup = initialKnots.length-2;
         // while(initialKnots[iBorneInf] == initialKnots[0])
         while (areEqual(initialKnots[iBorneInf], initialKnots[0], 1e-10))
            iBorneInf++;
         if (iBorneInf <= degree)
            iBorneInf = degree-iBorneInf+1;
         else
            iBorneInf=0;//on a alors iBorneInf valeurs a rajouter en debut de tableau

         // while(initialKnots[iBorneSup] == initialKnots[initialKnots.length-1])
         while (areEqual(initialKnots[iBorneSup], initialKnots[initialKnots.length-1], 1e-10))
            iBorneSup--;
         if (iBorneSup >= initialKnots.length-1-degree)
            iBorneSup = degree+1-(initialKnots.length-1-iBorneSup);
         else
            iBorneSup = 0; //on a alors iBorneSup valeurs a rajouter en fin de tableau

         //add computed values
         this.knots = new double[initialKnots.length + iBorneInf + iBorneSup];
         myX = new double[x.length + iBorneInf + iBorneSup];
         myY = new double[y.length + iBorneInf + iBorneSup];
         for (int i = 0; i < iBorneInf; i++) {
            this.knots[i] = initialKnots[0];
            myX[i] = x[0];
            myY[i] = y[0];
         }
         for(int i = 0; i<initialKnots.length; i++)
            this.knots[iBorneInf + i] = initialKnots[i];
         for(int i = 0; i<x.length; i++) {
            myX[iBorneInf + i] = x[i];
            myY[iBorneInf + i] = y[i];
         }
         for(int i = 0; i<iBorneSup; i++) {
            this.knots[this.knots.length-1 - i] = initialKnots[initialKnots.length-1];
            myX[myX.length-1-i] = x[x.length-1];
            myY[myY.length-1-i] = y[y.length-1];
         }
      }
   }

   public double evalX(double u) {
      int k = Misc.getTimeInterval (knots, 0, knots.length - 1, u);
      if (k >= myX.length) // set to last control point
         k = myX.length-1;
      
      double[][] X = new double[degree+1][myX.length];

      for(int i = k-degree; i<=k; i++)
         X[0][i] = myX[i];

      for(int j=1; j<= degree; j++) {
         for(int i = k-degree+j; i <= k; i++) {
            double aij = (u - this.knots[i]) / (this.knots[i+1+degree-j] - this.knots[i]);
            X[j][i] = (1-aij) * X[j-1][i-1] + aij * X[j-1][i];
         }
      }
      return X[degree][k];
   }

   public double evalY(double u) {
      int k = Misc.getTimeInterval (knots, 0, knots.length - 1, u);
      if (k >= myY.length) // set to last control point
         k = myY.length-1;
      
      double[][] Y = new double[degree+1][myX.length];

      for(int i = k-degree; i<=k; i++)
         Y[0][i] = myY[i];
      for(int j=1; j<= degree; j++) {
         for(int i = k-degree+j; i <= k; i++) {
            double aij = (u - this.knots[i]) / (this.knots[i+1+degree-j] - this.knots[i]);
            Y[j][i] = (1-aij) * Y[j-1][i-1] + aij * Y[j-1][i];
         }
      }
      return Y[degree][k];
   }

   /**
    * Checks if two doubles {@code a} and {@code b} are equal with tolerance level.
    * 
    * @param a
    * @param b
    * @param tol absolute tolerance level
    * @return 
    */
   private static boolean areEqual(double a, double b, double tol) {
      return Math.abs(a - b) < tol;
   }
   
   private static double[] computeN(double[] U, int degree, double u, int np1) {
      double[] N = new double[np1];

      // special cases at bounds
      if (areEqual(u, U[0], 1e-10)) {
         N[0] = 1.0;
         return N;
      }
      else if (areEqual(u, U[U.length-1], 1e-10)) {
         N[N.length-1] = 1.0;
         return N;
      }

      // find the knot index k such that U[k]<= u < U[k+1]
      int k = Misc.getTimeInterval (U, 0, U.length - 1, u);

      N[k] = 1.0;
      for(int d = 1; d <= degree; d++) {
         N[k-d] = N[k-d+1] * (U[k+1] - u) / (U[k+1] - U[k-d+1]);
         for(int i = k-d+1; i<= k-1; i++)
            N[i] = (u - U[i]) / (U[i+d]-U[i]) * N[i] + ((U[i+d+1] - u)/(U[i+d+1] - U[i+1])) * N[i+1];
         N[k] = (u - U[k]) / (U[k+d] - U[k]) * N[k];
      }
      return N;
   }

}