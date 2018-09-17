/*
 * Class:        MultivariateFunction
 * Description:  Represents a function of multiple variables.
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

/**
 * Represents a function of multiple variables. This interface specifies a
 * method `evaluate` that computes a @f$g(\mathbf{x})@f$ function, where
 * @f$\mathbf{x}=(x_0,…,x_{d-1})\in\mathbb{R}^d@f$. It also specifies a
 * method `evaluateGradient` for computing its gradient
 * @f$\nabla g(\mathbf{x})@f$.
 *
 * The dimension @f$d@f$ can be fixed or variable. When @f$d@f$ is fixed, the
 * methods specified by this interface always take the same number of
 * arguments. This is the case, for example, with a ratio of two variables.
 * When @f$d@f$ is variable, the implementation can compute the function for
 * a vector @f$\mathbf{x}@f$ of any length. This can happen for a product or
 * sum of variables.
 *
 * The methods of this interface take a variable number of arguments to
 * accomodate the common case of fixed dimension with more convenience; the
 * programmer can call the method without creating an array. For the generic
 * case, however, one can replace the arguments with an array.
 *
 * <div class="SSJ-bigskip"></div>
 */
public interface MultivariateFunction {

/**
 * Returns @f$d@f$, the dimension of the function computed by this
 * implementation. If the dimension is not fixed, this method must return a
 * negative value.
 *  @return the dimension.
 */
public int getDimension();

   /**
    * Computes the function @f$g(\mathbf{x})@f$ for the vector `x`. The
    * length of the given array must correspond to the dimension of this
    * function. The method must compute and return the result of the
    * function without modifying the elements in `x` since the array can
    * be reused for further computation.
    *  @param x            a vector @f$\mathbf{x}@f$.
    *  @return the value of @f$g(\mathbf{x})@f$.
    *
    *  @exception NullPointerException if `x` is `null`.
    *  @exception IllegalArgumentException if `x.length` does not
    * correspond to the dimension of this function.
    */
   public double evaluate (double... x);

   /**
    * Computes @f$\partial g(\mathbf{x})/\partial x_i@f$, the derivative
    * of @f$g(\mathbf{x})@f$ with respect to @f$x_i@f$. The length of the
    * given array must correspond to the dimension of this function. The
    * method must compute and return the result of the derivative without
    * modifying the elements in `x` since the array can be reused for
    * further computations, e.g., the gradient @f$\nabla g(\mathbf{x})@f$.
    *  @param i            the variable to derive with respect to.
    *  @param x            a vector @f$\mathbf{x}@f$.
    *  @return the value of the partial derivative.
    *
    *  @exception NullPointerException if `x` is `null`.
    *  @exception IllegalArgumentException if `x.length` does not
    * correspond to the dimension of this function.
    *  @exception IndexOutOfBoundsException if `i` is negative or greater
    * than or equal to the dimension of this function.
    */
   public double evaluateGradient (int i, double... x);

}