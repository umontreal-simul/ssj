/*
 * Class:        NameConflictException
 * Description:  An exception 
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