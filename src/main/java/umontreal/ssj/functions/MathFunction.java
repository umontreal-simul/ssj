/*
 * Class:        MathFunction
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
 * This interface should be implemented by classes which represent univariate
 * mathematical functions. It is used to pass an arbitrary function of one
 * variable as argument to another function. For example, it is used in
 * @ref umontreal.ssj.util.RootFinder to find the zeros of a function.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public interface MathFunction {

   /**
    * Returns the value of the function evaluated at @f$x@f$.
    *  @param x            value at which the function is evaluated
    *  @return function evaluated at `x`
    */
   public double evaluate (double x);

}