/*
 * Class:        CachedDataWriter
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
import java.io.*;
import java.util.LinkedList;

/**
 * This abstract class implements shared functionality for data writers that
 * store all fields in memory before outputing them with
 * umontreal.ssj.util.io.DataWriter.close.
 *
 * <div class="SSJ-bigskip"></div>
 */
public abstract class CachedDataWriter extends AbstractDataWriter {

   // don't use a map because ordering is important
   protected LinkedList<DataField> fields;

   protected LinkedList<DataField> getFields() {
      return fields;
   }

   /**
    * Class constructor.
    */
   public CachedDataWriter() {
      this.fields = new LinkedList<DataField>();
   }

   /**
    * @name Writing atomic data
    * @{
    */

   /**
    * Writes an atomic string field. Writes an anonymous field if `label`
    * is `null`.
    */
   public void write (String label, String s) throws IOException {
      fields.add(new DataField(label, s));
   }

   /**
    * Writes an atomic 32-bit integer (big endian). Writes an anonymous
    * field if `label` is `null`.
    */
   public void write (String label, int a) throws IOException {
      fields.add(new DataField(label, a));
   }

   /**
    * Writes an atomic 32-bit float (big endian). Writes an anonymous
    * field if `label` is `null`.
    */
   public void write (String label, float a) throws IOException {
      fields.add(new DataField(label, a));
   }

   /**
    * Writes an atomic 64-bit double (big endian). Writes an anonymous
    * field if `label` is `null`.
    */
   public void write (String label, double a) throws IOException {
      fields.add(new DataField(label, a));
   }

   /**
    * @}
    */

   /**
    * @name Writing one-dimensional arrays
    * @{
    */

   /**
    * Writes the first `n` elements of a one-dimensional array of strings.
    * Writes an anonymous field if `label` is `null`.
    */
   public void write (String label, String[] a, int n) throws IOException {
      fields.add(new DataField(label, a.clone(), n));
   }

   /**
    * Writes the first `n` elements of a one-dimensional array of 32-bit
    * integers (big endian). Writes an anonymous field if `label` is
    * `null`.
    */
   public void write (String label, int[] a, int n) throws IOException {
      fields.add(new DataField(label, a.clone(), n));
   }

   /**
    * Writes the first `n` elements of a one-dimensional array of 32-bit
    * floats (big endian). Writes an anonymous field if `label` is `null`.
    */
   public void write (String label, float[] a, int n) throws IOException {
      fields.add(new DataField(label, a.clone(), n));
   }

   /**
    * Writes the first `n` elements of a one-dimensional array of 64-bit
    * doubles (big endian). Writes an anonymous field if `label` is
    * `null`.
    */
   public void write (String label, double[] a, int n) throws IOException {
      fields.add(new DataField(label, a.clone(), n));
   }

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
   public void write (String label, String[][] a) throws IOException {
      fields.add(new DataField(label, a.clone()));
   }

   /**
    * Writes a two-dimensional array of 32-bit integers (big endian).
    * Writes an anonymous field if `label` is `null`.
    */
   public void write (String label, int[][] a) throws IOException {
      fields.add(new DataField(label, a.clone()));
   }

   /**
    * Writes a two-dimensional array of 32-bit floats (big endian). Writes
    * an anonymous field if `label` is `null`.
    */
   public void write (String label, float[][] a) throws IOException {
      fields.add(new DataField(label, a.clone()));
   }

   /**
    * Writes a two-dimensional array of 64-bit doubles (big endian).
    * Writes an anonymous field if `label` is `null`.
    */
   public void write (String label, double[][] a) throws IOException {
      fields.add(new DataField(label, a.clone()));
   }

}

/**
 * @}
 */