/*
 * Class:        UnuranException
 * Description:  unchecked exception for errors inside the UNURAN package
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
package umontreal.ssj.randvar;

/**
 * This type of unchecked exception is thrown when an error occurs *inside*
 * the UNURAN package. Usually, such an exception will come from the native
 * side.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_unuran
 */
public class UnuranException extends RuntimeException {

   /**
    * Constructs a new generic UNURAN exception.
    */
   public UnuranException() {
      super();
   }

   /**
    * Constructs a UNURAN exception with the error message `message`
    *  @param message      error message describing the problem that
    *                      occurred
    */
   public UnuranException (String message) {
      super (message);
   }
}