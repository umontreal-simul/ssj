/*
 * Class:        CloneableRandomStream
 * Description:  
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
package umontreal.ssj.rng;

/**
 *  @ref CloneableRandomStream extends  @ref RandomStream and  Cloneable. All
 * classes that implements this interface are able to produce cloned objects.
 *
 * The cloned object is entirely independent of the older odject. Moreover
 * the cloned object has all the same properties as the older one. All his
 * seeds are duplicated, and therefore both generators will produce the same
 * random number sequence.
 *
 * <div class="SSJ-bigskip"></div>
 */
public interface CloneableRandomStream extends RandomStream, Cloneable {

   /**
    * Clones the current object and returns its copy.
    *  @return A deep copy of the current object
    */
   public CloneableRandomStream clone();
 
}