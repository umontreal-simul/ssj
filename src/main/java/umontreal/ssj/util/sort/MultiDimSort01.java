/*
 * Class:        MultiDimSort01
 * Description:  Represents a sort on multidimensional arrays of MultiDimComparable.
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
 import java.util.Comparator;
 import java.util.Arrays;

/**
 * This interface extends  @ref MultiDimSort<T> to implement multivariate
 * sorting algorithms that sort points of @f$d@f$ dimensions in the unit
 * hypercube @f$[0, 1)^d@f$.
 *
 * <div class="SSJ-bigskip"></div>
 */
public interface MultiDimSort01 <T extends MultiDim01> extends 
   MultiDimSort<T> {
}
