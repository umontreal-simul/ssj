/*
 * Class:        IdentityMathFunction
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