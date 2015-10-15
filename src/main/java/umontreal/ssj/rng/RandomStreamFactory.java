/*
 * Class:        RandomStreamFactory
 * Description:  random stream factory that can construct instances of
                 a given type of random stream
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

import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.rng.MRG32k3a;

/**
 * Represents a random stream factory capable of constructing instances of a
 * given type of random stream by invoking the  #newInstance method each time
 * a new random stream is needed, instead of invoking directly the specific
 * constructor of the desired type. Hence, if several random streams of a
 * given type (class) must be constructed at different places in a large
 * simulation program, and if we decide to change the type of stream in the
 * future, there is no need to change the code at those different places.
 * With the random stream factory, the class-specific code for constructing
 * these streams appears at a single place, where the factory is constructed.
 *
 * The class  @ref BasicRandomStreamFactory provides an implementation of
 * this interface.
 *
 * <div class="SSJ-bigskip"></div>
 */
public interface RandomStreamFactory {

/**
 * Constructs and returns a new random stream. If the instantiation of the
 * random stream fails, this method throws a
 * @ref RandomStreamInstantiationException.
 *  @return the newly-constructed random stream.
 *
 *  @exception RandomStreamInstantiationException if the new random stream
 * cannot be instantiated.
 */
public RandomStream newInstance();

}