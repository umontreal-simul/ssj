/*
 * Class:        ArrayOfObservationListener
 * Description:  Array of observation listener
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Éric Buist 
 * @since        2007
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
package umontreal.ssj.stat.list;
import umontreal.ssj.stat.StatProbe;

/**
 * Represents an object that can listen to observations broadcast by lists of
 * statistical probes.
 *
 * <div class="SSJ-bigskip"></div>
 */
public interface ArrayOfObservationListener {

/**
 * Receives the new array of observations `x` broadcast by the list of
 * statistical probes `listOfProbes`.
 *  @param listOfProbes the list of statistical probes broadcasting the
 *                      observation.
 *  @param x            the array of observations being broadcast.
 */
public void newArrayOfObservations (ListOfStatProbes<?> listOfProbes,
                                       double[] x);

}