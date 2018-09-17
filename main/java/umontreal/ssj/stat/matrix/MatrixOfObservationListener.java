/*
 * Class:        MatrixOfObservationListener
 * Description:  Observation listener
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Éric Buist
 * @since        2006

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
package umontreal.ssj.stat.matrix;
import umontreal.ssj.stat.StatProbe;
import cern.colt.matrix.DoubleMatrix2D;

/**
 * Represents an object that can listen to observations broadcast by matrices
 * of statistical probes.
 *
 * <div class="SSJ-bigskip"></div>
 */
public interface MatrixOfObservationListener {

/**
 * Receives the new matrix of observations `x` broadcast by the matrix of
 * statistical probes `matrixOfProbes`.
 *  @param matrixOfProbes the matrix of statistical probes broadcasting the
 *                        observation.
 *  @param x            the matrix of observations being broadcast.
 */
public void newMatrixOfObservations (MatrixOfStatProbes<?> matrixOfProbes,
                                        DoubleMatrix2D x);

}