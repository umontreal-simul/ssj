/*
 * Class:        BasicRandomStreamFactory
 * Description:  basic random stream factory
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

import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.rng.MRG32k3a;

/**
 * Represents a basic random stream factory that can constructs new instances
 * of a given  @ref RandomStream implementation via the  #newInstance method.
 * The class name of the implementation to be used must be passed to the
 * constructor as a `String`, which must be the name of a nullary constructor
 * of a  @ref RandomStream object (i.e., a constructor that has no
 * parameters). The streams are constructed by the factory by reflection from
 * this `String`.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class BasicRandomStreamFactory implements RandomStreamFactory {
   private Class rsClass;

/**
 * Constructs a new basic random stream factory with random stream class
 * `rsClass`. The supplied class object must represent an implementation of
 * @ref RandomStream and must provide a nullary constructor. For example, to
 * construct a factory producing  @ref MRG32k3a random streams, this
 * constructor must be called with `MRG33k3a.class`.
 *  @param rsClass      the random stream class being used.
 *  @exception NullPointerException if `rsClass` is `null`.
 *  @exception IllegalArgumentException if `rsClass` does not represent an
 * implementation of  @ref RandomStream, or does not provide a nullary
 * constructor.
 */
public BasicRandomStreamFactory (Class rsClass) {
      checkRandomStreamClass (rsClass);
      this.rsClass = rsClass;
   }

   /**
    * Returns the random stream class associated with this object.
    *  @return the associated random stream class.
    */
   public Class getRandomStreamClass() {
      return rsClass;
   }

   /**
    * Sets the associated random stream class to `rsClass`. The supplied
    * class object must represent an implementation of  @ref RandomStream
    * and must provide a nullary constructor.
    *  @param rsClass      the random stream class being used.
    *  @exception NullPointerException if `rsClass` is `null`.
    *  @exception IllegalArgumentException if `rsClass` does not represent
    * an implementation of  @ref RandomStream, or does not provide a
    * nullary constructor.
    */
   public void setRandomStreamClass (Class rsClass) {
      checkRandomStreamClass (rsClass);
      this.rsClass = rsClass;
   }


   private void checkRandomStreamClass (final Class rsClass) {
      if (!RandomStream.class.isAssignableFrom (rsClass))
         throw new IllegalArgumentException
            ("The random class must implement the RandomStream interface");
      try {
         rsClass.getConstructor (new Class[0]);
      }
      catch (NoSuchMethodException nme) {
         throw new IllegalArgumentException
            ("The random stream class " + rsClass.getName() + " does not have a " +
             "nullary public constructor.");
      }
   }

   public RandomStream newInstance() {
      try {
         return (RandomStream)rsClass.newInstance();
      }
      catch (IllegalAccessException iae) {
         throw new RandomStreamInstantiationException
            (this, "Cannot access constructor for random stream class " + rsClass.getName(),
             iae);
      }
      catch (InstantiationException ie) {
         throw new RandomStreamInstantiationException
            (this, "Cannot instantiate random stream class "
             + rsClass.getName(), ie);
      }
      catch (Exception e) {
         throw new RandomStreamInstantiationException
            (this, "Exception while calling the nullary constructor for random stream class "
             + rsClass.getName(), e);
      }
   }

   public String toString() {
      return "Basic random stream factory constructing streams of class " + rsClass.getName();
   }
}