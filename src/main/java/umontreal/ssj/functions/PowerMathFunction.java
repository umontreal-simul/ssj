/*
 * Class:        PowerMathFunction
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
 * Represents a function computing @f$(af(x) + b)^p@f$ for a user-defined
 * function @f$f(x)@f$ and power @f$p@f$.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class PowerMathFunction implements MathFunction

,
      MathFunctionWithFirstDerivative {
   private MathFunction func;
   private double a, b;
   private double power;

/**
 * Constructs a new power function for function `func` and power `power`. The
 * values of the constants are @f$a=1@f$ and @f$b=0@f$.
 *  @param func         the function @f$f(x)@f$.
 *  @param power        the power @f$p@f$.
 */
public PowerMathFunction (MathFunction func, double power) {
      this (func, 1, 0, power);
   }

   /**
    * Constructs a new power function for function `func`, power `power`,
    * and constants `a` and `b`.
    *  @param func         the function @f$f(x)@f$.
    *  @param power        the power @f$p@f$.
    *  @param a            the multiplicative constant.
    *  @param b            the additive constant.
    */
   public PowerMathFunction (MathFunction func, double a, double b, double power) {
      if (func == null)
         throw new NullPointerException ();
      this.func = func;
      this.a = a;
      this.b = b;
      this.power = power;
   }

   /**
    * Returns the function @f$f(x)@f$.
    *  @return the function.
    */
   public MathFunction getFunction () {
      return func;
   }

   /**
    * Returns the value of @f$a@f$.
    *  @return the value of @f$a@f$.
    */
   public double getA () {
      return a;
   }

   /**
    * Returns the value of @f$b@f$.
    *  @return the value of @f$b@f$.
    */
   public double getB () {
      return b;
   }

   /**
    * Returns the power @f$p@f$.
    *  @return the power.
    */
   public double getPower () {
      return power;
   }


   public double derivative (double x) {
      final double fder = MathFunctionUtil.derivative (func, x);
      return getA()*getPower()*Math.pow (getA() * func.evaluate (x) + getB(), getPower() - 1)*fder;
   }

   public double evaluate (double x) {
      final double v = func.evaluate (x);
      return Math.pow (a * v + b, power);
   }
}