/*
 * Class:        CloneableRandomStream
 * Description:  
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
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