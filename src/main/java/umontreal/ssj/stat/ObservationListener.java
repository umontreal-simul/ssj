/*
 * Class:        ObservationListener
 * Description:  Observation listener
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
package umontreal.ssj.stat;

/**
 * Represents an object that can listen to observations broadcast by
 * statistical probes.
 *
 * <div class="SSJ-bigskip"></div>
 */
public interface ObservationListener {

/**
 * Receives the new observation `x` broadcast by `probe`.
 *  @param probe        the statistical probe broadcasting the observation.
 *  @param x            the observation being broadcast.
 */
public void newObservation (StatProbe probe, double x);

}