/*
 * Class:        PiecewiseConstantFunction
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

import java.util.Arrays;

/**
 * Represents a piecewise-constant function.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class PiecewiseConstantFunction implements MathFunction {
   private double[] x;
   private double[] y;

/**
 * Constructs a new piecewise-constant function with @f$X@f$ and @f$Y@f$
 * coordinates given by `x` and `y`.
 *  @param x            the @f$X@f$ coordinates.
 *  @param y            the @f$Y@f$ coordinates.
 */
public PiecewiseConstantFunction (double[] x, double[] y) {
      if (x.length != y.length)
         throw new IllegalArgumentException();
      this.x = x.clone ();
      this.y = y.clone ();
   }

   /**
    * Returns the @f$X@f$ coordinates of the function.
    *  @return the @f$X@f$ coordinates of the function.
    */
   public double[] getX() {
      return x.clone ();
   }

   /**
    * Returns the @f$Y@f$ coordinates of the function.
    *  @return the @f$Y@f$ coordinates of the function.
    */
   public double[] getY() {
      return y.clone ();
   }


   public double evaluate (double x) {
      final int idx = Arrays.binarySearch (this.x, x);
      if (idx >= 0)
         return y[idx];
      final int insertionPoint = -(idx + 1);
      return y[insertionPoint];
   }
}