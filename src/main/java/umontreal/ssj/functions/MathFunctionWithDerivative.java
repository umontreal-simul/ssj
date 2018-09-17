/*
 * Class:        MathFunctionWithDerivative
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

/**
 * Represents a mathematical function whose @f$n@f$th derivative can be
 * computed using  #derivative(double,int).
 *
 * <div class="SSJ-bigskip"></div>
 */
public interface MathFunctionWithDerivative extends MathFunction {

/**
 * Computes (or estimates) the @f$n@f$th derivative of the function at point
 * `x`. For @f$n=0@f$, this returns the result of
 * umontreal.ssj.functions.MathFunction.evaluate(double).
 *  @param x            the point to evaluate the derivate to.
 *  @param n            the order of the derivative.
 *  @return the resulting derivative.
 *
 *  @exception IllegalArgumentException if `n` is negative or 0.
 */
public double derivative (double x, int n);
}