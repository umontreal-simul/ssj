/*
 * Class:        DataWriter
 * Description:  Data writer interface
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       David Munger 
 * @since        August 2009

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
package umontreal.ssj.util.io;
import java.io.IOException;

/**
 * Data writer interface.
 *
 * <div class="SSJ-bigskip"></div>
 */
public interface DataWriter {

   /**
    * @name Writing atomic data
    * @{
    */

   /**
    * Writes an atomic string field. Writes an anonymous field if `label`
    * is `null`.
    */
   public void write (String label, String s) throws IOException;

   /**
    * Writes an atomic 32-bit integer (big endian). Writes an anonymous
    * field if `label` is `null`.
    */
   public void write (String label, int a) throws IOException;

   /**
    * Writes an atomic 32-bit float (big endian). Writes an anonymous
    * field if `label` is `null`.
    */
   public void write (String label, float a) throws IOException;

   /**
    * Writes an atomic 64-bit double (big endian). Writes an anonymous
    * field if `label` is `null`.
    */
   public void write (String label, double a) throws IOException;

   /**
    * @}
    */

   /**
    * @name Writing one-dimensional arrays
    * @{
    */

   /**
    * Writes a one-dimensional array of strings. Writes an anonymous field
    * if `label` is `null`.
    */
   public void write (String label, String[] a) throws IOException;

   /**
    * Writes the first `n` elements of a one-dimensional array of strings.
    * Writes an anonymous field if `label` is `null`.
    */
   public void write (String label, String[] a, int n) throws IOException;

   /**
    * Writes a one-dimensional array of 32-bit integers (big endian).
    * Writes an anonymous field if `label` is `null`.
    */
   public void write (String label, int[] a) throws IOException;

   /**
    * Writes the first `n` elements of a one-dimensional array of 32-bit
    * integers (big endian). Writes an anonymous field if `label` is
    * `null`.
    */
   public void write (String label, int[] a, int n) throws IOException;

   /**
    * Writes a one-dimensional array of 32-bit floats (big endian). Writes
    * an anonymous field if `label` is `null`.
    */
   public void write (String label, float[] a) throws IOException;

   /**
    * Writes the first `n` elements of a one-dimensional array of 32-bit
    * floats (big endian). Writes an anonymous field if `label` is `null`.
    */
   public void write (String label, float[] a, int n) throws IOException;

   /**
    * Writes a one-dimensional array of 64-bit doubles (big endian).
    * Writes an anonymous field if `label` is `null`.
    */
   public void write (String label, double[] a) throws IOException;

   /**
    * Writes the first `n` elements of a one-dimensional array of 64-bit
    * doubles (big endian). Writes an anonymous field if `label` is
    * `null`.
    */
   public void write (String label, double[] a, int n) throws IOException;

   /**
    * @}
    */

   /**
    * @name Writing two-dimensional arrays
    * @{
    */

   /**
    * Writes a two-dimensional array of strings. Writes an anonymous field
    * if `label` is `null`.
    */
   public void write (String label, String[][] a) throws IOException;

   /**
    * Writes a two-dimensional array of 32-bit integers (big endian).
    * Writes an anonymous field if `label` is `null`.
    */
   public void write (String label, int[][] a) throws IOException;

   /**
    * Writes a two-dimensional array of 32-bit floats (big endian). Writes
    * an anonymous field if `label` is `null`.
    */
   public void write (String label, float[][] a) throws IOException;

   /**
    * Writes a two-dimensional array of 64-bit doubles (big endian).
    * Writes an anonymous field if `label` is `null`.
    */
   public void write (String label, double[][] a) throws IOException;

   /**
    * @}
    */

   /**
    * @name Other methods
    * @{
    */

   /**
    * Flushes any pending data and closes the output stream.
    */
   public void close() throws IOException;

}

/**
 * @}
 */