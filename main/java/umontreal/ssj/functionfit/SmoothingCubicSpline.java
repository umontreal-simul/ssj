/*
 * Class:        SmoothingCubicSpline
 * Description:  smoothing cubic spline algorithm of Schoenberg
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
   import umontreal.ssj.functions.*;
   import umontreal.ssj.functions.Polynomial;
import umontreal.ssj.functions.MathFunctionWithFirstDerivative;
import umontreal.ssj.functions.MathFunctionWithDerivative;
import umontreal.ssj.functions.MathFunctionWithIntegral;
import umontreal.ssj.functions.SquareMathFunction;
import umontreal.ssj.functions.MathFunctionUtil;

/**
 * Represents a cubic spline with nodes at @f$(x_i, y_i)@f$ computed with the
 * smoothing cubic spline algorithm of Schoenberg @cite mDEB78a,
 * @cite mPOL93a&thinsp;. A smoothing cubic spline is made of @f$n+1@f$ cubic
 * polynomials. The @f$i@f$th polynomial of such a spline, for
 * @f$i=1,…,n-1@f$, is defined as @f$S_i(x)@f$ while the complete spline is
 * defined as
 * @f[
 *   S(x) = S_i(x), \qquad\mbox{ for }x\in[x_{i-1}, x_i].
 * @f]
 * For @f$x<x_0@f$ and @f$x>x_{n-1}@f$, the spline is not precisely defined,
 * but this class performs extrapolation by using @f$S_0@f$ and @f$S_n@f$
 * linear polynomials. The algorithm which calculates the smoothing spline is
 * a generalization of the algorithm for an interpolating spline. @f$S_i@f$
 * is linked to @f$S_{i+1}@f$ at @f$x_{i+1}@f$ and keeps continuity
 * properties for first and second derivatives at this point, therefore
 * @f$S_i(x_{i+1})=S_{i+1}(x_{i+1})@f$, @f$S’_i(x_{i+1})=S’_{i+1}(x_{i+1})@f$
 * and @f$S"_i(x_{i+1})=S"_{i+1}(x_{i+1})@f$.
 *
 * The spline is computed with a smoothing parameter @f$\rho\in[0, 1]@f$
 * which represents its accuracy with respect to the initial @f$(x_i, y_i)@f$
 * nodes. The smoothing spline minimizes
 * @f[
 *   L = \rho\sum_{i=0}^{n-1}{w_i}\left({y_i - S_i(x_i)}\right)^2 + (1-\rho)\int_{x_0}^{x_{n-1}}\left(S"(x)\right)^2dx
 * @f]
 * In fact, by setting @f$\rho= 1@f$, we obtain the interpolating spline;
 * and we obtain a linear function by setting @f$\rho= 0@f$. The weights
 * @f$w_i>0@f$, which default to 1, can be used to change the contribution of
 * each point in the error term. A large value @f$w_i@f$ will give a large
 * weight to the @f$i@f$th point, so the spline will pass closer to it. Here
 * is a small example that uses smoothing splines:
 *
 * @code
 *
 *    int n;
 *    double[] X = new double[n];
 *    double[] Y = new double[n];
 *    // here, fill arrays X and Y with n data points (x_i, y_i)
 *    // The points must be sorted with respect to x_i.
 *
 *    double rho = 0.1;
 *    SmoothingCubicSpline fit = new SmoothingCubicSpline(X, Y, rho);
 *
 *    int m = 40;
 *    double[] Xp = new double[m+1];       // Xp, Yp are spline points
 *    double[] Yp = new double[m+1];
 *    double h = (X[n-1] - X[0]) / m;      // step
 *
 *    for (int i = 0; i <= m; i++) {
 *       double z = X[0] + i * h;
 *       Xp[i] = z;
 *       Yp[i] = fit.evaluate(z);          // evaluate spline at z
 *    }
 *
 * @endcode
 *
 * <div class="SSJ-bigskip"></div>
 */
public class SmoothingCubicSpline implements MathFunction,
             MathFunctionWithFirstDerivative, MathFunctionWithDerivative,
             MathFunctionWithIntegral {

   private Polynomial[] splineVector;
   private double[] x, y, weight;
   private double rho;

   /**
    * Constructs a spline with nodes at @f$(x_i, y_i)@f$, with weights
    * @f$w_i@f$ and smoothing factor @f$\rho@f$ = `rho`. The @f$x_i@f$
    * *must* be sorted in increasing order.
    *  @param x            the @f$x_i@f$ coordinates.
    *  @param y            the @f$y_i@f$ coordinates.
    *  @param w            the weight for each point, must be @f$> 0@f$.
    *  @param rho          the smoothing parameter
    *  @exception IllegalArgumentException if `x`, `y` and `z` do not have
    * the same length, if rho has wrong value, or if the spline cannot be
    * calculated.
    */
   public SmoothingCubicSpline (double[] x, double[] y, double[] w,
                                double rho) {
      if (x.length != y.length)
         throw new IllegalArgumentException ("x.length != y.length");
      if (w != null && x.length != w.length)
         throw new IllegalArgumentException ("x.length != w.length");
      if (rho < 0 || rho > 1)
         throw new IllegalArgumentException ("rho not in [0, 1]");

      splineVector = new Polynomial[x.length+1];
      this.rho = rho;
      this.x = x.clone();
      this.y = y.clone();

      weight = new double[x.length];
      if (w == null) {
         for (int i = 0; i < weight.length; i++)
            weight[i] = 1.0;
      } else {
         for (int i = 0; i < weight.length; i++)
            weight[i] = w[i];
      }

      resolve();
   }

   /**
    * Constructs a spline with nodes at @f$(x_i, y_i)@f$, with weights
    * @f$= 1@f$ and smoothing factor @f$\rho@f$ = `rho`. The @f$x_i@f$
    * *must* be sorted in increasing order.
    *  @param x            the @f$x_i@f$ coordinates.
    *  @param y            the @f$y_i@f$ coordinates.
    *  @param rho          the smoothing parameter
    *  @exception IllegalArgumentException if `x` and `y` do not have the
    * same length, if rho has wrong value, or if the spline cannot be
    * calculated.
    */
   public SmoothingCubicSpline (double[] x, double[] y, double rho) {
      this (x, y, null, rho);
   }

   /**
    * Evaluates and returns the value of the spline at @f$z@f$.
    *  @param z            argument of the spline.
    *  @return value of spline.
    */
   public double evaluate (double z) {
      int i = getFitPolynomialIndex(z);
      if (i == 0)
         return splineVector[i].evaluate(z-x[0]);
      else
         return splineVector[i].evaluate(z-x[i-1]);
   }

   /**
    * Evaluates and returns the value of the integral of the spline from
    * @f$a@f$ to @f$b@f$.
    *  @param a            lower limit of integral.
    *  @param b            upper limit of integral.
    *  @return value of integral.
    */
   public double integral (double a, double b) {
      int iA = getFitPolynomialIndex(a);
      int iB = getFitPolynomialIndex(b);
      double retour;
      int i = 1;

      // Calcule l'integrale
      if (iA == iB) { // les deux valeurs sont sur le meme polynome
         retour = splineVector[iB].integral(a-x[iB], b-x[iB]);
      } else {
         if (iA == 0)
            retour = splineVector[iA].integral(a-x[iA], 0);
         else
            retour = splineVector[iA].integral(a-x[iA], x[iA+1]-x[iA]);
         for (i = iA+1; i<iB; i++) {
            retour += splineVector[i].integral(0, x[i+1]-x[i]);
         }
         retour += splineVector[iB].integral(0, b-x[iB]);
      }
      return retour;
   }

   /**
    * Evaluates and returns the value of the *first* derivative of the
    * spline at @f$z@f$.
    *  @param z            argument of the spline.
    *  @return value of first derivative.
    */
   public double derivative (double z) {
      int i = getFitPolynomialIndex(z);
      if (i == 0)
         return splineVector[i].derivative(z-x[0]);
      else
         return splineVector[i].derivative(z-x[i-1]);
   }

   /**
    * Evaluates and returns the value of the <em>n</em>-th derivative of
    * the spline at @f$z@f$.
    *  @param z            argument of the spline.
    *  @param n            order of the derivative.
    *  @return value of n-th derivative.
    */
   public double derivative (double z, int n)  {
      int i = getFitPolynomialIndex(z);
      if (i == 0)
         return splineVector[i].derivative(z - x[0], n);
      else
         return splineVector[i].derivative(z - x[i-1], n);
   }

   /**
    * Returns the @f$x_i@f$ coordinates for this spline.
    *  @return the @f$x_i@f$ coordinates.
    */
   public double[] getX() {
      return x.clone ();
   }

   /**
    * Returns the @f$y_i@f$ coordinates for this spline.
    *  @return the @f$y_i@f$ coordinates.
    */
   public double[] getY() {
      return y.clone ();
   }

   /**
    * Returns the weights of the points.
    *  @return the weights of the points.
    */
   public double[] getWeights() {
     return weight;
   }

   /**
    * Returns the smoothing factor used to construct the spline.
    *  @return the smoothing factor.
    */
   public double getRho() {
      return rho;
   }

   /**
    * Returns a table containing all fitting polynomials.
    *  @return Table containing the fitting polynomials.
    */
   public Polynomial[] getSplinePolynomials() {
      return splineVector.clone();
   }

   /**
    * Returns the index of @f$P@f$, the
    * @ref umontreal.ssj.functions.Polynomial instance used to evaluate
    * @f$x@f$, in an `ArrayList` table instance returned by
    * `getSplinePolynomials()`. This index @f$k@f$ gives also the interval
    * in table **X** which contains the value @f$x@f$ (i.e. such that
    * @f$x_k < x \leq x_{k+1}@f$).
    *  @return Index of the polynomial check with x in the Polynomial list
    * returned by  #getSplinePolynomials
    */
   public int getFitPolynomialIndex (double x) {
      // Algorithme de recherche binaire legerement modifie
      int j = this.x.length-1;
      if (x > this.x[j])
         return j+1;
      int tmp = 0;
      int i = 0;

      while (i+1 != j) {
         if (x > this.x[tmp]) {
            i = tmp;
            tmp = i+(j-i)/2;
         } else {
            j = tmp;
            tmp = i+(j-i)/2;
         }
         if (j == 0) // le point est < a x_0, on sort
            i--;
      }
      return i+1;
   }


   private void resolve () {
   /*
      taken from D.S.G Pollock's paper, "Smoothing with Cubic Splines",
      Queen Mary, University of London (1993)
      http://www.qmw.ac.uk/~ugte133/PAPERS/SPLINES.PDF
   */

      double[] h = new double[x.length];
      double[] r = new double[x.length];
      double[] u = new double[x.length];
      double[] v = new double[x.length];
      double[] w = new double[x.length];
      double[] q = new double[x.length+1];
      double[] sigma = new double[weight.length];

      for (int i = 0; i < weight.length; i++) {
         if (weight[i] <= 0.0)
            sigma[i] = 1.0e100;
         else
            sigma[i] = 1.0/Math.sqrt(weight[i]);
      }

      double mu;
      int n = x.length-1;

      if (rho <= 0)
         mu = 1.0e100;   // arbitrary large number to avoid 1/0
      else
         mu = 2 * (1 - rho)/(3 * rho);

      h[0] = x[1] - x[0];
      r[0] = 3/h[0];
      for (int i = 1; i < n; i++) {
         h[i] = x[i+1] - x[i];
         r[i] = 3/h[i];
         q[i] = 3 * (y[i+1] - y[i])/h[i] - 3 * (y[i] - y[i-1])/h[i - 1];
      }

      for (int i = 1; i < n; i++) {
         u[i] = r[i-1]*r[i-1] * sigma[i-1]
               + (r[i - 1] + r[i])*(r[i - 1] + r[i]) * sigma[i] +
                 r[i]*r[i] * sigma[i+1];
         u[i] = mu * u[i] + 2 * (x[i+1] - x[i-1]);
         v[i] = -(r[i - 1] + r[i]) * r[i] * sigma[i] - r[i] * (r[i] +
                 r[i+1]) * sigma[i+1];
         v[i] = mu * v[i] + h[i];
         w[i] = mu * r[i] * r[i+1] * sigma[i+1];
      }
      q = Quincunx(u, v, w, q);

      // extrapolation a gauche
      double[] params = new double[4];
      double dd;
      params[0] = y[0] - mu * r[0] * q[1] * sigma[0];
      dd = y[1] - mu * ((-r[0] - r[1]) * q[1] + r[1] * q[2]) * sigma[1];
      params[1] = (dd - params[0])/h[0] - q[1] * h[0]/3;
      splineVector[0] = new Polynomial(params);

      // premier polynome
      params[0] = y[0] - mu * r[0] * q[1] * sigma[0];
      dd = y[1] - mu * ((-r[0] - r[1]) * q[1] + r[1] * q[2]) * sigma[1];
      params[3] = q[1]/(3 * h[0]);
      params[2] = 0;
      params[1] = (dd - params[0])/h[0] - q[1] * h[0]/3;
      splineVector[1] = new Polynomial(params);

      // les polynomes suivants
      int j;
      for (j = 1; j < n; j++) {
         params[3] = (q[j + 1] - q[j])/(3 * h[j]);
         params[2] = q[j];
         params[1] = (q[j] + q[j - 1]) * h[j - 1] + splineVector[j].getCoefficient(1);
         params[0] = r[j - 1] * q[j - 1] + (-r[j-1] - r[j]) * q[j] + r[j] * q[j + 1];
         params[0] = y[j] - mu * params[0] * sigma[j];
         splineVector[j+1] = new Polynomial(params);
      }

      // extrapolation a droite
      j = n;
      params[3] = 0;
      params[2] = 0;
      params[1] = splineVector[j].derivative(x[x.length-1]-x[x.length-2]);
      params[0] = splineVector[j].evaluate(x[x.length-1]-x[x.length-2]);
      splineVector[n+1] = new Polynomial(params);
   }


   private static double[] Quincunx (double[] u, double[] v,
                                     double[] w, double[] q) {
      u[0] = 0;
      v[1] = v[1]/u[1];
      w[1] = w[1]/u[1];
      int j;

      for (j = 2; j < u.length-1; j++) {
         u[j] = u[j] - u[j - 2] * w[j-2]*w[j-2] - u[j - 1] * v[j-1]*v[j-1];
         v[j] = (v[j] - u[j - 1] * v[j-1] * w[j-1])/u[j];
         w[j] = w[j]/u[j];
      }

      // forward substitution
      q[1] = q[1] - v[0] * q[0];
      for (j = 2; j < u.length-1; j++) {
         q[j] = q[j] - v[j - 1] * q[j - 1] - w[j - 2] * q[j - 2];
      }
      for (j = 1; j < u.length-1; j++) {
         q[j] = q[j]/u[j];
      }

      // back substitution
      q[u.length-1] = 0;
      for (j = u.length-3; j > 0; j--) {
         q[j] = q[j] - v[j] * q[j + 1] - w[j] * q[j + 2];
      }
      return q;
   }
}