/*
 * Class:        ShiftedMathFunction
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