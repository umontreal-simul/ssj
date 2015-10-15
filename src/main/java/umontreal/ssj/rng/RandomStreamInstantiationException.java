/*
 * Class:        RandomStreamInstantiationException
 * Description:  thrown when a random stream factory cannot instantiate a stream
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
 * This exception is thrown when a random stream factory cannot instantiate a
 * stream on a call to its  umontreal.ssj.rng.RandomStreamFactory.newInstance
 * method.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class RandomStreamInstantiationException extends RuntimeException {
   private RandomStreamFactory factory;

/**
 * Constructs a new random stream instantiation exception with no message, no
 * cause, and thrown by the given `factory`.
 *  @param factory      the random stream factory which thrown the exception.
 */
public RandomStreamInstantiationException (RandomStreamFactory factory) {
      super();
      this.factory = factory;
   }

   /**
    * Constructs a new random stream instantiation exception with the
    * given `message`, no cause, and concerning `factory`.
    *  @param factory      the random stream factory concerned by the
    *                      exception.
    *  @param message      the error message describing the exception.
    */
   public RandomStreamInstantiationException (RandomStreamFactory factory,
                                         String message) {
      super (message);
      this.factory = factory;
   }

   /**
    * Constructs a new random stream instantiation exception with no
    * message, the given `cause`, and concerning `factory`.
    *  @param factory      the random stream factory concerned by the
    *                      exception.
    *  @param cause        the cause of the exception.
    */
   public RandomStreamInstantiationException (RandomStreamFactory factory,
                                         Throwable cause) {
      super (cause);
      this.factory = factory;
   }

   /**
    * Constructs a new random stream instantiation exception with the
    * given `message`, the supplied `cause`, and concerning `factory`.
    *  @param factory      the random stream factory concerned by the
    *                      exception.
    *  @param message      the error message describing the exception.
    *  @param cause        the cause of the exception.
    */
   public RandomStreamInstantiationException (RandomStreamFactory factory,
                                         String message, Throwable cause) {
      super (message, cause);
      this.factory = factory;
   }

   /**
    * Returns the random stream factory concerned by this exception.
    *  @return the random stream factory concerned by this exception.
    */
   public RandomStreamFactory getRandomStreamFactory() {
      return factory;
   }

   /**
    * Returns a short description of the exception. If
    * #getRandomStreamFactory returns `null`, this calls `super.toString`.
    * Otherwise, the result is the concatenation of: a) the name of the
    * actual class of the exception;<br>b) the string `": For random
    * stream factory "`;<br>c) the result of  #getRandomStreamFactory
    * `.toString()`;<br>d) if  java.lang.Throwable.getMessage is
    * non-<tt>null</tt>, `", "` followed by the result of
    * java.lang.Throwable.getMessage.
    *  @return a string representation of the exception.
    */
   public String toString() {
      if (factory == null)
         return super.toString();

      StringBuffer sb = new StringBuffer (getClass().getName());
      sb.append (": For random stream factory ");
      sb.append (factory.toString());
      String msg = getMessage();
      if (msg != null)
         sb.append (", ").append (msg);
      return sb.toString();
   }
}