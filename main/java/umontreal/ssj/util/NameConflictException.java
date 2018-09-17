/*
 * Class:        NameConflictException
 * Description:  An exception 
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
package umontreal.ssj.util;

/**
 * This exception is thrown by a  @ref ClassFinder when two or more fully
 * qualified class names can be associated with a simple class name.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class NameConflictException extends Exception {
   private static final long serialVersionUID = -5124156035520217708L;
   private ClassFinder finder;
   private String name;

   /**
    * Constructs a new name conflict exception.
    */
   public NameConflictException() {}

   /**
    * Constructs a new name conflict exception with message `message`.
    *  @param message      the error message.
    */
   public NameConflictException (String message) {
      super (message);
   }

   /**
    * Constructs a new name conflict exception with class finder `finder`,
    * simple name `name`, and message `message`.
    *  @param finder       the class finder in which the name conflict
    *                      occurred.
    *  @param name         the simple conflicting name.
    *  @param message      the message describint the conflict.
    */
   public NameConflictException (ClassFinder finder, String name,
                                 String message) {
      super (message);
      this.finder = finder;
      this.name = name;
   }

   /**
    * Returns the class finder associated with this exception.
    *  @return the associated class finder.
    */
   public ClassFinder getClassFinder() {
      return finder;
   }

   /**
    * Returns the simple name associated with this exception.
    *  @return the associated simple name.
    */
   public String getName() {
      return name;
   }
}