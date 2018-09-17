/*
 * Class:        ShiftedMathFunction
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
 * Represents a function computing @f$f(x) - \delta@f$ for a user-defined
 * function @f$f(x)@f$ and shift @f$\delta@f$.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class ShiftedMathFunction implements MathFunction

,
      MathFunctionWithFirstDerivative, MathFunctionWithDerivative,
      MathFunctionWithIntegral {
   MathFunction func;
   double delta;

/**
 * Constructs a new function shifting the function `func` by a shift `delta`.
 *  @param func         the function.
 *  @param delta        the shift.
 */
public ShiftedMathFunction (MathFunction func, double delta) {
      if (func == null)
         throw new NullPointerException ();
      this.func = func;
      this.delta = delta;
   }

   /**
    * Returns the function @f$f(x)@f$.
    *  @return the function.
    */
   public MathFunction getFunction () {
      return func;
   }

   /**
    * Returns the shift @f$\delta@f$ = `delta`.
    *  @return the shift.
    */
   public double getDelta () {
      return delta;
   }


   public double evaluate (double x) {
      return func.evaluate (x) - delta;
   }

   public double derivative (double x) {
      return MathFunctionUtil.derivative (func, x);
   }

   public double derivative (double x, int n) {
      return MathFunctionUtil.derivative (func, x, n);
   }

   public double integral (double a, double b) {
      return MathFunctionUtil.integral (func, a, b) - (b - a)*getDelta();
   }
}