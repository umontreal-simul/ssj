/*
 * Class:        ExponentialDistFromMean
 * Description:  exponential distribution
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
package umontreal.ssj.probdist;

/**
 * Extends the  @ref ExponentialDist class with a constructor accepting as
 * argument the mean @f$1/\lambda@f$ instead of the rate @f$\lambda@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class ExponentialDistFromMean extends ExponentialDist {

   /**
    * Constructs a new exponential distribution with mean `mean`.
    *  @param mean         the required mean.
    */
   public ExponentialDistFromMean (double mean) {
      super (1.0 / mean);
   }

   /**
    * Calls  umontreal.ssj.probdist.ExponentialDist.setLambda(double) with
    * argument `1/mean` to change the mean of this distribution.
    *  @param mean         the new mean.
    */
   public void setMean (double mean) {
      setLambda (1.0 / mean);
   }
}