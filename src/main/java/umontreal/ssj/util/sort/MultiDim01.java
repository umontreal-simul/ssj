/*
 * Class:        MultiDim01
 * Description:  Represents an object which can be compared in many dimensions.
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
package umontreal.ssj.util.sort;

/**
 * This interface represents a point or array of @f$d@f$ dimensions in a unit
 * hypercube @f$[0, 1)^d@f$. The value of the @f$j@f$th dimension can be
 * accessed with the method  {@link #getCoordinate() getCoordinate(j)}.
 *
 * <div class="SSJ-bigskip"></div>
 */
public interface MultiDim01 {

   /**
    * This method returns the number dimensions of this point.
    */
   public int dimension();

   /**
    * Returns the @f$d@f$ coordinates of this point.
    */
   public double[] getPoint();

   /**
    * Returns the value of @f$j@f$th coordinate (or dimension). This value
    * should be in the interval @f$[0, 1)@f$. If @f$j@f$ is outside the
    * range @f$\{0,…,d-1\}@f$, this method should throw an
    * IllegalArgumentException.
    */
   public double getCoordinate (int j);

}
