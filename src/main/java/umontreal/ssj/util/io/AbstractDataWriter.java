/*
 * Class:        AbstractDataWriter
 * Description:  
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       David Munger 
 * @since        August 2009
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