/*
 * Class:        IdentityMathFunction
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
 * Represents the identity function @f$f(x)=x@f$.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class IdentityMathFunction implements MathFunction

,
      MathFunctionWithFirstDerivative, MathFunctionWithDerivative,
      MathFunctionWithIntegral {
   public double evaluate (double x) {
      return x;
   }
   
   public double derivative (double x) {
      return 1;
   }

   public double derivative (double x, int n) {
      return n > 1 ? 0 : 1;
   }

   public double integral (double a, double b) {
      return (b*b - a*a) / 2;
   }
}