/*
 * Class:        MathFunction
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