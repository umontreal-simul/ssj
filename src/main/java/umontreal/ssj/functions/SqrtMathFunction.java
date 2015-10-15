/*
 * Class:        SqrtMathFunction
 * Description:  function computing the square root of another function
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
 * Represents a function computing the square root of another function
 * @f$f(x)@f$.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class SqrtMathFunction implements MathFunction {
   private MathFunction func;

/**
 * Computes and returns the square root of the function `func`.
 *  @param func         the function to compute square root for.
 */
public SqrtMathFunction (MathFunction func) {
      super ();
      if (func == null)
         throw new NullPointerException();
      this.func = func;
   }

   /**
    * Returns the function associated with this object.
    *  @return the associated function.
    */
   public MathFunction getFunction() {
      return func;
   }


   public double evaluate (double x) {
      return Math.sqrt (func.evaluate (x));
   }
}