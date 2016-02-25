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