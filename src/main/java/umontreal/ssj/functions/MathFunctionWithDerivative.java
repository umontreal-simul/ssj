/*
 * Class:        MathFunctionWithDerivative
 * Description:
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Éric Buist
 * @since

 * SSJ is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, either version 3 of the License, or
 * any later version.

 * SSJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * A copy of the GNU General Public License is available at
   <a href="http://www.gnu.org/licenses">GPL licence site</a>.
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