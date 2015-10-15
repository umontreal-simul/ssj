/*
 * Class:        ObservationListener
 * Description:  Observation listener
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