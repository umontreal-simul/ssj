/*
 * Class:        AverageMathFunction
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
 * Represents a function computing the average of several functions. Let
 * @f$f_0(x), …, f_{n-1}(x)@f$ be a set of @f$n@f$ functions. This function
 * represents the average
 * @f[
 *   f(x)=\frac{1}{n}\sum_{i=0}^{n-1} f_i(x).
 * @f]
 * <div class="SSJ-bigskip"></div>
 */
public class AverageMathFunction implements MathFunction

,
   MathFunctionWithFirstDerivative, MathFunctionWithDerivative,
   MathFunctionWithIntegral {
   private MathFunction[] func;

/**
 * Constructs a function computing the average of the functions in the array
 * `func`.
 *  @param func         the array of functions to average.
 */
public AverageMathFunction (MathFunction... func) {
      if (func == null)
         throw new NullPointerException();
      this.func = func.clone ();
   }

   /**
    * Returns the functions being averaged.
    *  @return the averaged functions.
    */
   public MathFunction[] getFunctions() {
      return func.clone ();
   }


   public double evaluate (double x) {
      double sum = 0;
      for (final MathFunction fi : func)
         sum += fi.evaluate (x);
      return sum / func.length;
   }
   
   public double derivative (double x, int n) {
      if (n < 0)
         throw new IllegalArgumentException ("n must be greater than or equal to 0");
      if (n == 0)
         return evaluate (x);
      double sum = 0;
      for (final MathFunction fi : func)
         sum += MathFunctionUtil.derivative (fi, x, n);
      return sum / func.length;
   }

   public double derivative (double x) {
      double sum = 0;
      for (final MathFunction fi : func)
         sum += MathFunctionUtil.derivative (fi, x);
      return sum / func.length;
   }

   public double integral (double a, double b) {
      double sum = 0;
      for (final MathFunction fi : func)
         sum += MathFunctionUtil.integral (fi, a, b);
      return sum / func.length;
   }
}