/*
 * Class:        RatioFunction
 * Description:  Represents a function computing a ratio of two values.
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
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
package umontreal.ssj.util;

/**
 * Represents a function computing a ratio of two values.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class RatioFunction implements MultivariateFunction {
   private double zeroOverZero = Double.NaN;

   /**
    * Constructs a new ratio function.
    */
   public RatioFunction() {}

   /**
    * Constructs a new ratio function that returns `zeroOverZero` for the
    * special case of @f$0/0@f$. See the  #getZeroOverZeroValue method for
    * more information. The default value of `zeroOverZero` is
    * `Double.NaN`.
    *  @param zeroOverZero the value for @f$0/0@f$.
    */
   public RatioFunction (double zeroOverZero) {
      this.zeroOverZero = zeroOverZero;
   }

   /**
    * Returns the value returned by  #evaluate(double[]) in the case where
    * the @f$0/0@f$ function is calculated. The default value for
    * @f$0/0@f$ is `Double.NaN`.
    *
    * Generally, @f$0/0@f$ is undefined, and therefore associated with the
    * `Double.NaN` constant, meaning *not-a-number*. However, in certain
    * applications, it can be defined differently to accomodate some
    * special cases. For exemple, in a queueing system, if there are no
    * arrivals, no customers are served, lost, queued, etc. As a result,
    * many performance measures of interest turn out to be @f$0/0@f$.
    * Specifically, the loss probability, i.e., the ratio of lost
    * customers over the number of arrivals, should be 0 if there is no
    * arrival; in this case, @f$0/0@f$ means 0. On the other hand, the
    * service level, i.e., the fraction of customers waiting less than a
    * fixed threshold, could be fixed to 1 if there is no arrival.
    *  @return the value for @f$0/0@f$.
    */
   public double getZeroOverZeroValue() {
      return zeroOverZero;
   }

   /**
    * Sets the value returned by  #evaluate(double[]) for the undefined
    * function @f$0/0@f$ to `zeroOverZero`. See  #getZeroOverZeroValue for
    * more information.
    *  @param zeroOverZero the new value for @f$0/0@f$.
    */
   public void setZeroOverZeroValue (double zeroOverZero) {
      this.zeroOverZero = zeroOverZero;
   }


   public int getDimension() {
      return 2;
   }

   public double evaluate (double... x) {
      if (x.length != 2)
         throw new IllegalArgumentException
            ("Invalid length of x");
      if (x[0] == 0 && x[1] == 0)
         return zeroOverZero;
      return x[0]/x[1];
   }

   public double evaluateGradient (int i, double... x) {
      if (x.length != 2)
         throw new IllegalArgumentException
            ("Invalid length of x");
      switch (i) {
      case 0: return 1.0/x[1];
      case 1: return -x[0]/(x[1]*x[1]);
      default: throw new IndexOutOfBoundsException
         ("Invalid value of i: " + i); 
      }
   }
}