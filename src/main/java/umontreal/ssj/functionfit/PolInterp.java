/*
 * Class:        PolInterp
 * Description:  polynomial that interpolates through a set of points
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

import java.io.Serializable;
import umontreal.ssj.functions.Polynomial;

import umontreal.ssj.util.PrintfFormat;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

/**
 * Represents a polynomial that interpolates through a set of points. More
 * specifically, let @f$(x_0,y_0), …, (x_n, y_n)@f$ be a set of points and
 * @f$p(x)@f$ the constructed polynomial of degree @f$n@f$. Then, for
 * @f$i=0,…,n@f$, @f$p(x_i)=y_i@f$.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class PolInterp extends Polynomial implements Serializable {
   private static final long serialVersionUID = -710451931485296501L;
   private static final Algebra alg = new Algebra ();
   private double[] x;
   private double[] y;

   /**
    * Constructs a new polynomial interpolating through the given points
    * `(x[0], y[0]), ..., (x[n], y[n])`. This constructs a polynomial of
    * degree `n` from `n+1` points.
    *  @param x            the @f$x@f$ coordinates of the points.
    *  @param y            the @f$y@f$ coordinates of the points.
    *  @exception NullPointerException if `x` or `y` are `null`.
    *  @exception IllegalArgumentException if the lengths of `x` and `y`
    * are different, or if less than two points are specified.
    */
   public PolInterp (double[] x, double[] y) {
      super (getCoefficients (x, y));
      this.x = x.clone ();
      this.y = y.clone ();
   }

   /**
    * Computes and returns the coefficients the polynomial interpolating
    * through the given points `(x[0], y[0]), ..., (x[n], y[n])`. This
    * polynomial has degree `n` and there are `n+1` coefficients.
    *  @param x            the @f$x@f$ coordinates of the points.
    *  @param y            the @f$y@f$ coordinates of the points.
    *  @return the coefficients the interpolating polynomial.
    */
   public static double[] getCoefficients (double[] x, double[] y) {
      if (x.length != y.length)
         throw new IllegalArgumentException (
               "x and y must have the same length");
      if (x.length <= 1)
         throw new IllegalArgumentException ("At least two points are needed");
      final DoubleMatrix2D u = new DenseDoubleMatrix2D (x.length, x.length);
      for (int i = 0; i < x.length; i++) {
         double v = 1;
         for (int j = 0; j < x.length; j++) {
            u.setQuick (i, j, v);
            v *= x[i];
         }
      }
      final DoubleMatrix2D yMat = new DenseDoubleMatrix2D (x.length, 1);
      yMat.viewColumn (0).assign (y);
      final DoubleMatrix2D bMat = alg.solve (u, yMat);
      return bMat.viewColumn (0).toArray ();
   }

   /**
    * Returns the @f$x@f$ coordinates of the interpolated points.
    *  @return the @f$x@f$ coordinates of the interpolated points.
    */
   public double[] getX() {
      return x.clone ();
   }

   /**
    * Returns the @f$y@f$ coordinates of the interpolated points.
    *  @return the @f$y@f$ coordinates of the interpolated points.
    */
   public double[] getY() {
      return y.clone ();
   }

   /**
    * Makes a string representation of a set of points.
    *  @param x            the @f$x@f$ coordinates of the points.
    *  @param y            the @f$y@f$ coordinates of the points.
    *  @return the string representing the points.
    */
   public static String toString (double[] x, double[] y) {
      final StringBuilder sb = new StringBuilder ("Points: ");
      for (int i = 0; i < x.length; i++) {
         if (i > 0)
            sb.append (", ");
         final String xstr = PrintfFormat.format (8, 3, 3, x[i]);
         final String ystr = PrintfFormat.format (8, 3, 3, y[i]);
         sb.append ('(').append (xstr).append (", ").append (ystr).append (')');
      }
      return sb.toString ();
   }

   /**
    * Calls  {@link #toString(double[],double[]) toString(double[],
    * double[])} with the associated points.
    *  @return a string containing the points.
    */
   @Override
   public String toString() {
      return toString (x, y);
   }


   public PolInterp clone() {
      final PolInterp p = (PolInterp) super.clone ();
      p.x = x.clone ();
      p.y = y.clone ();
      return p;
   }
}