/*
 * Class:        MathFunctionUtil
 * Description:
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Éric Buist
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
package umontreal.ssj.functions;

import umontreal.ssj.util.Misc;

/**
 * Provides utility methods for computing derivatives and integrals of
 * functions.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class MathFunctionUtil {

/**
 * Step length in @f$x@f$ to compute derivatives. Default: @f$10^{-6}@f$.
 */
public static double H = 1e-6;

   private MathFunctionUtil() {}

   // For Gauss-Lobatto: nodes Cg and weights Wg
   private static final double[] Cg = { 0, 0.17267316464601142812, 0.5,
                                           0.82732683535398857188, 1 };
   private static final double[] Wg = { 0.05, 0.27222222222222222222,
                  0.35555555555555555555, 0.27222222222222222222, 0.05 };


   private static double[] fixBounds (MathFunction func, double a,
                                      double b, int numIntervals) {
   // For functions which are 0 on parts of [a, b], shorten the interval
   // [a, b] to the non-zero part of f(x). Returns the shortened interval.

       final double h = (b - a)/numIntervals;
       double x = b;
       while ((0 == func.evaluate (x)) && (x > a))
           x -= h;
       if (x < b)
           b = x + h;

       x = a;
       while ((0 == func.evaluate (x)) && (x < b))
           x += h;
       if (x > a)
          a = x - h;
       double[] D = {a, b};
       return D;
   }

/**
 * Default number of intervals for Simpson’s integral.
 */
public static int NUMINTERVALS = 1024;

   /**
    * Returns the first derivative of the function `func` evaluated at
    * `x`. If the given function implements
    * @ref MathFunctionWithFirstDerivative, this method calls
    * {@link MathFunctionWithFirstDerivative.derivative()
    * derivative(double)}. Otherwise, if the function implements
    * @ref MathFunctionWithDerivative, this method calls
    * {@link MathFunctionWithDerivative.derivative() derivative(double,
    * int)}. If the function does not implement any of these two
    * interfaces, the method uses
    * {@link #finiteCenteredDifferenceDerivative()
    * finiteCenteredDifferenceDerivative(MathFunction, double, double)} to
    * obtain an estimate of the derivative.
    *  @param func         the function to derivate.
    *  @param x            the evaluation point.
    *  @return the first derivative.
    */
   public static double derivative (MathFunction func, double x) {
      if (func instanceof MathFunctionWithFirstDerivative)
         return ((MathFunctionWithFirstDerivative)func).derivative (x);
      else if (func instanceof MathFunctionWithDerivative)
         return ((MathFunctionWithDerivative)func).derivative (x, 1);
      else
         return finiteCenteredDifferenceDerivative (func, x, H);
   }

   /**
    * Returns the @f$n@f$th derivative of function `func` evaluated at
    * `x`. If @f$n=0@f$, this returns @f$f(x)@f$. If @f$n=1@f$, this calls
    * {@link #derivative() derivative(MathFunction, double)} and returns
    * the resulting first derivative. Otherwise, if the function
    * implements  @ref MathFunctionWithDerivative, this method calls
    * {@link MathFunctionWithDerivative.derivative() derivative(double,
    * int)}. If the function does not implement this interface, the method
    * uses  {@link #finiteCenteredDifferenceDerivative()
    * finiteCenteredDifferenceDerivative(MathFunction, double, int,
    * double)} if @f$n@f$ is even, or
    * {@link #finiteDifferenceDerivative()
    * finiteDifferenceDerivative(MathFunction, double, int, double)} if
    * @f$n@f$ is odd, to obtain a numerical approximation of the
    * derivative.
    *  @param func         the function to derivate.
    *  @param x            the evaluation point.
    *  @param n            the order of the derivative.
    *  @return the @f$n@f$th derivative.
    */
   public static double derivative (MathFunction func, double x, int n) {
      if (n == 0)
         return func.evaluate (x);
      else if (n == 1)
         return derivative (func, x);
      else if (func instanceof MathFunctionWithDerivative)
         return ((MathFunctionWithDerivative)func).derivative (x, n);
      else if (n % 2 == 0)
         return finiteCenteredDifferenceDerivative (func, x, n, H);
      else
         return finiteDifferenceDerivative (func, x, n, H);
   }

   /**
    * Computes and returns an estimate of the @f$n@f$th derivative of the
    * function @f$f(x)@f$. This method estimates
    * @f[
    *   \frac{d^nf(x)}{dx^n},
    * @f]
    * the @f$n@f$th derivative of @f$f(x)@f$ evaluated at @f$x@f$. This
    * method first computes @f$f_i=f(x+i\epsilon)@f$, for @f$i=0,…,n@f$,
    * with @f$\epsilon=h^{1/n}@f$. The estimate is then given by
    * @f$\Delta^nf_0/h@f$, where @f$\Delta^nf_i=\Delta^{n-1}f_{i+1} -
    * \Delta^{n-1}f_i@f$, and @f$\Delta f_i = f_{i+1} - f_i@f$.
    *  @param func         the function to derivate.
    *  @param x            the evaluation point.
    *  @param n            the order of the derivative.
    *  @param h            the error.
    *  @return the estimate of the derivative.
    */
   public static double finiteDifferenceDerivative (
                 MathFunction func, double x, int n, double h) {
      if (n < 0)
         throw new IllegalArgumentException
         ("n must not be negative");
      if (n == 0)
         return func.evaluate (x);
      final double err = Math.pow (h, 1.0 / n);
      final double[] values = new double[n+1];
      for (int i = 0; i < values.length; i++)
         values[i] = func.evaluate (x + i*err);
      for (int j = 0; j < n; j++) {
         for (int i = 0; i < n - j; i++)
            values[i] = values[i + 1] - values[i];
      }
      return values[0] / h;
   }

   /**
    * Returns @f$(f(x + h) - f(x - h))/(2h)@f$, an estimate of the first
    * derivative of @f$f(x)@f$ using centered differences.
    *  @param func         the function to derivate.
    *  @param x            the evaluation point.
    *  @param h            the error.
    *  @return the estimate of the first derivative.
    */
   public static double finiteCenteredDifferenceDerivative (
                 MathFunction func, double x, double h) {
      final double fplus = func.evaluate (x + h);
      final double fminus = func.evaluate (x - h);
      return (fplus - fminus) / (2*h);
   }

   /**
    * Computes and returns an estimate of the @f$n@f$th derivative of the
    * function @f$f(x)@f$ using finite centered differences. If @f$n@f$ is
    * even, this method returns
    * {@link #finiteDifferenceDerivative(MathFunction,double,int,double)
    * finiteDifferenceDerivative(func, x - \epsilon*n/2, n, h)}, with
    * @f$h=\epsilon^n@f$.
    *  @param func         the function to derivate.
    *  @param x            the evaluation point.
    *  @param n            the order of the derivative.
    *  @param h            the error.
    *  @return the estimate of the derivative.
    */
   public static double finiteCenteredDifferenceDerivative (
                 MathFunction func, double x, int n, double h) {
      if (n < 0)
         throw new IllegalArgumentException
         ("n must not be negative");
      if (n == 0)
         return func.evaluate (x);
      if (n % 2 == 1)
         throw new IllegalArgumentException ("n must be even");
      final double err = Math.pow (h, 1.0 / n);
      return finiteDifferenceDerivative (func, x - n*err / 2, n, h);
   }

   /**
    * Removes any point `(NaN, y)` or `(x, NaN)` from `x` and `y`, and
    * returns a 2D array containing the filtered points. This method
    * filters each pair (<tt>x[i]</tt>, <tt>y[i]</tt>) containing at least
    * one NaN element. It constructs a 2D array containing the two
    * filtered arrays, whose size is smaller than or equal to `x.length`.
    *  @param x            the @f$X@f$ coordinates.
    *  @param y            the @f$Y@f$ coordinates.
    *  @return the filtered @f$X@f$ and @f$Y@f$ arrays.
    */
   public static double[][] removeNaNs (double[] x, double[] y) {
      if (x.length != y.length)
         throw new IllegalArgumentException();
      int numNaNs = 0;
      for (int i = 0; i < x.length; i++)
         if (Double.isNaN (x[i]) || Double.isNaN (y[i]))
            ++numNaNs;
      if (numNaNs == 0)
         return new double[][] { x, y };
      final double[] nx = new double[x.length - numNaNs];
      final double[] ny = new double[y.length - numNaNs];
      int j = 0;
      for (int i = 0; i < x.length; i++)
         if (!Double.isNaN (x[i]) && !Double.isNaN (y[i])) {
            nx[j] = x[i];
            ny[j++] = y[i];
         }
      return new double[][] { nx, ny };
   }

   /**
    * Returns the integral of the function `func` over @f$[a, b]@f$. If
    * the given function implements  @ref MathFunctionWithIntegral, this
    * returns
    * {@link umontreal.ssj.functions.MathFunctionWithIntegral.integral()
    * integral(double, double)}. Otherwise, this calls
    * {@link #simpsonIntegral() simpsonIntegral(MathFunction, double,
    * double, int)} with  #NUMINTERVALS intervals.
    *  @param func         the function to integrate.
    *  @param a            the lower bound.
    *  @param b            the upper bound.
    *  @return the value of the integral.
    */
   public static double integral (MathFunction func, double a, double b) {
      if (func instanceof MathFunctionWithIntegral)
         return ((MathFunctionWithIntegral)func).integral (a, b);
      else
         return simpsonIntegral (func, a, b, NUMINTERVALS);
   }

   /**
    * Computes and returns an approximation of the integral of `func` over
    * @f$[a, b]@f$, using the Simpson’s @f$1/3@f$ method with
    * `numIntervals` intervals. This method estimates
    * @f[
    *   \int_a^b f(x)dx,
    * @f]
    * where @f$f(x)@f$ is the function defined by `func` evaluated at
    * @f$x@f$, by dividing @f$[a, b]@f$ in @f$n=@f$&nbsp;`numIntervals`
    * intervals of length @f$h=(b - a)/n@f$. The integral is estimated by
    * @f[
    *   \frac{h}{3}(f(a)+4f(a+h)+2f(a+2h)+4f(a+3h)+\cdots+f(b))
    * @f]
    * This method assumes that @f$a\le b<\infty@f$, and @f$n@f$ is even.
    *  @param func         the function being integrated.
    *  @param a            the left bound
    *  @param b            the right bound.
    *  @param numIntervals the number of intervals.
    *  @return the approximate value of the integral.
    */
   public static double simpsonIntegral (MathFunction func, double a,
                                         double b, int numIntervals) {
      if (numIntervals % 2 != 0)
         throw new IllegalArgumentException
         ("numIntervals must be an even number");
      if (Double.isInfinite (a) || Double.isInfinite (b) ||
         Double.isNaN (a) || Double.isNaN (b))
         throw new IllegalArgumentException
             ("a and b must not be infinite or NaN");
      if (b < a)
         throw new IllegalArgumentException ("b < a");
      if (a == b)
         return 0;
      double[] D = fixBounds (func, a, b, numIntervals);
      a = D[0];
      b = D[1];
      final double h = (b - a) / numIntervals;
      final double h2 = 2*h;
      final int m = numIntervals / 2;
      double sum = 0;
      for (int i = 0; i < m - 1; i++) {
         final double x = a + h + h2*i;
         sum += 4*func.evaluate (x) + 2*func.evaluate (x + h);
      }
      sum += func.evaluate (a) + func.evaluate (b) + 4*func.evaluate (b - h);
      return sum * h / 3;
   }

   /**
    * Computes and returns a numerical approximation of the integral of
    * @f$f(x)@f$ over @f$[a, b]@f$, using Gauss-Lobatto adaptive
    * quadrature with 5 nodes, with tolerance `tol`. This method estimates
    * @f[
    *   \int_a^b f(x)dx,
    * @f]
    * where @f$f(x)@f$ is the function defined by `func`. Whenever the
    * estimated error is larger than `tol`, the interval @f$[a, b]@f$ will
    * be halved in two smaller intervals, and the method will recursively
    * call itself on the two smaller intervals until the estimated error
    * is smaller than `tol`.
    *  @param func         the function being integrated.
    *  @param a            the left bound
    *  @param b            the right bound.
    *  @param tol          error.
    *  @return the approximate value of the integral.
    */
   public static double gaussLobatto (MathFunction func, double a, double b,
                                      double tol) {
      if (b < a)
         throw new IllegalArgumentException ("b < a");
      if (Double.isInfinite (a) || Double.isInfinite (b) ||
          Double.isNaN (a) || Double.isNaN (b))
         throw new IllegalArgumentException ("a or b is infinite or NaN");
      if (a == b)
         return 0;
      double r0 = simpleGaussLob (func, a, b);
      final double h = (b - a)/2;
      double r1 = simpleGaussLob (func, a, a + h) +
                  simpleGaussLob (func, a + h, b);
      double maxi = Math.max(1.0, Math.abs(r1));
      if (Math.abs(r0 - r1) <= tol*maxi)
         return r1;
      return gaussLobatto (func, a, a + h, tol) +
             gaussLobatto (func, a + h, b, tol);
   }


   private static double simpleGaussLob (MathFunction func, double a, double b) {
      // Gauss-Lobatto integral over [a, b] with 5 nodes
      if (a == b)
         return 0;
      final double h = b - a;
      double sum = 0;
      for (int i = 0; i < 5; i++) {
         sum += Wg[i] * func.evaluate(a + h*Cg[i]);
      }
      return h*sum;
   }

   /**
    * Similar to method  {@link #gaussLobatto() gaussLobatto(MathFunction,
    * double, double, double)}, but also returns in `T[0]` the
    * subintervals of integration, and in `T[1]`, the partial values of
    * the integral over the corresponding subintervals. Thus `T[0][0]`
    * @f$= x_0 = a@f$ and `T[0][n]` @f$=x_n =b@f$; `T[1][i]` contains the
    * value of the integral over the subinterval @f$[x_{i-1}, x_i]@f$; we
    * also have `T[1][0]` @f$=0@f$. The sum over all `T[1][i]`, for
    * @f$i=1, …, n@f$ gives the value of the integral over @f$[a,b]@f$,
    * which is the value returned by this method. *WARNING:* The user
    * *must reserve* the 2 elements of the first dimension (<tt>T[0]</tt>
    * and <tt>T[1]</tt>) before calling this method.
    *  @param func         function being integrated
    *  @param a            left bound of interval
    *  @param b            right bound of interval
    *  @param tol          error
    *  @param T            @f$(x,y)@f$ = (values of partial
    *                      intervals,partial values of integral)
    *  @return value of the integral
    */
   public static double gaussLobatto (MathFunction func, double a, double b,
                                      double tol, double[][] T) {
      if (b < a)
         throw new IllegalArgumentException ("b < a");
      if (a == b) {
         T[0] = new double [1];
         T[1] = new double [1];
         T[0][0] = a;
         T[1][0] = 0;
         return 0;
      }

      int n = 1;         // initial capacity
      T[0] = new double [n];
      T[1] = new double [n];
      T[0][0] = a;
      T[1][0] = 0;
      int[] s = {0};    // actual number of intervals
      double res = innerGaussLob (func, a, b, tol, T, s);
      n = s[0] + 1;
      double[] temp = new double[n];
      System.arraycopy (T[0], 0, temp, 0, n);
      T[0] = temp;
      temp = new double[n];
      System.arraycopy (T[1], 0, temp, 0, n);
      T[1] = temp;
      return res;
   }


   private static double innerGaussLob (MathFunction func, double a, double b,
                                        double tol, double[][] T, int[] s) {
      double r0 = simpleGaussLob (func, a, b);
      final double h = (b - a) / 2;
      double r1 = simpleGaussLob (func, a, a + h) +
                  simpleGaussLob (func, a + h, b);
      if (Math.abs(r0 - r1) <= tol) {
         ++s[0];
         int len = s[0];
         if (len >= T[0].length) {
            double[] temp = new double[2 * len];
            System.arraycopy (T[0], 0, temp, 0, len);
            T[0] = temp;
            temp = new double[2 * len];
            System.arraycopy (T[1], 0, temp, 0, len);
            T[1] = temp;
         }
         T[0][len] = b;
         T[1][len] = r1;
         return r1;
      }

      return innerGaussLob (func, a, a + h, tol, T, s) +
             innerGaussLob (func, a + h, b, tol, T, s);
   }

}