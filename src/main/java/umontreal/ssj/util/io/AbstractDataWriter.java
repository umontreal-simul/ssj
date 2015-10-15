/*
 * Class:        AbstractDataWriter
 * Description:  
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
 * This abstract class implements shared functionality for data writers.
 *
 * <div class="SSJ-bigskip"></div>
 */
public abstract class AbstractDataWriter implements DataWriter {

   /**
    * @name Writing one-dimensional arrays
    * @{
    */

   /**
    * Writes a one-dimensional array of strings. If `label` is `null`,
    * writes an anonymous field.
    */
   public void write (String label, String[] a) throws IOException {
      write(label, a, a.length);
   }

   /**
    * Writes a one-dimensional array of 32-bit integers (big endian). If
    * `label` is `null`, writes an anonymous field.
    */
   public void write (String label, int[] a) throws IOException {
      write(label, a, a.length);
   }

   /**
    * Writes a one-dimensional array of 32-bit floats (big endian). If
    * `label` is `null`, writes an anonymous field.
    */
   public void write (String label, float[] a) throws IOException {
      write(label, a, a.length);
   }

   /**
    * Writes a one-dimensional array of 64-bit doubles (big endian). If
    * `label` is `null`, writes an anonymous field.
    */
   public void write (String label, double[] a) throws IOException {
      write(label, a, a.length);
   }

}

/**
 * @}
 */