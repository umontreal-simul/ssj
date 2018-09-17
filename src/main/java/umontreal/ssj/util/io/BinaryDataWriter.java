/*
 * Class:        BinaryDataWriter
 * Description:  Binary data writer
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

/**
 * Binary data writer.
 *
 * Stores a sequence of fields in binary file, which can be either atoms or
 * arrays, each of which having the following format:
 *
 * <ul><li>
 * Field label:
 *
 * <ul><li>
 * Pipe character (<tt>|</tt>)
 * </li>
 * <li>
 * Label length (32-bit integer, big endian)
 * </li>
 * <li>
 * Label string (array of bytes of the specified length)
 * </li>
 * </ul>
 * </li>
 * <li>
 * Field type (byte):
 *
 * <ul><li>
 * `i` (32-bit integer)
 * </li>
 * <li>
 * `f` (32-bit float)
 * </li>
 * <li>
 * `d` (64-bit double)
 * </li>
 * <li>
 * `S` (string)
 * </li>
 * </ul>
 * </li>
 * <li>
 * Number of dimensions (8-bit integer)
 * </li>
 * <li>
 * Dimensions (array of 32-bit integers, big endian)
 * </li>
 * <li>
 * Field data (in the specified format, big endian)
 * </li>
 * </ul>
 *
 * In the case of an atomic field, the number of dimensions is set to zero.
 *
 * A string field is stored in the following format:
 *
 * <ul><li>
 * String length (32-bit integer)
 * </li>
 * <li>
 * Array of bytes of the specified length
 * </li>
 * </ul>
 *
 * Also supports anonymous fields (fields with an empty label).
 *
 * Arrays up to two dimensions are supported.
 *
 * Modules for reading data exported with this class are available in Java (
 * @ref BinaryDataReader ), Matlab and Python (numpy).
 *
 *  Provide links for the import modules.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class BinaryDataWriter extends AbstractDataWriter {
   protected DataOutputStream out;
   
   /**
    * Utility method to write string data.
    *
    */    
   protected void writeStringData(String s) throws IOException {
      if (s != null) {
         out.writeInt(s.length());
         out.writeBytes(s);
      }
      else {
         out.writeInt(0);
      }
   }
   
   /**
    * Starts a new field by writing its label.
    *
    * @param label   name of the field (can be {@code null})
    *
    */
   protected void writeLabel(String label) throws IOException {
      out.writeByte(TYPECHAR_LABEL);
      writeStringData(label);
   }

   /**
    * @name Fields
    * @{
    */

   /**
    * *Field-type* symbol indicating a label (it more accurately a field
    * separator symbol).
    */
   public final static byte TYPECHAR_LABEL   = '|';

   /**
    * *Field-type* symbol indicating `String` data.
    */
   public final static byte TYPECHAR_STRING  = 'S';

   /**
    * *Field-type* symbol indicating `int` data.
    */
   public final static byte TYPECHAR_INTEGER = 'i';

   /**
    * *Field-type* symbol indicating `float` data.
    */
   public final static byte TYPECHAR_FLOAT   = 'f';

   /**
    * *Field-type* symbol indicating `double` data.
    */
   public final static byte TYPECHAR_DOUBLE  = 'd';

   /**
    * @}
    */

   /**
    * Data will be output to the file with the specified name.
    *  @param filename     name of the file to be created or appended to
    *  @param append       an existing file with the specified name will
    *                      be appended to if `true` or truncated if
    *                      `false`
    */
   public BinaryDataWriter (String filename, boolean append)
         throws IOException {
      this.out = new DataOutputStream(new FileOutputStream(filename, append));
   }

   /**
    * Data will be output to the specified file.
    *  @param file         file to be created or appended to
    *  @param append       an existing file with the specified name will
    *                      be appended to if `true` or truncated if
    *                      `false`
    */
   public BinaryDataWriter (File file, boolean append) throws IOException {
      this.out = new DataOutputStream(new FileOutputStream(file, append));
   }

   /**
    * Truncates any existing file with the specified name.
    *  @param filename     name of the file to be created
    */
   public BinaryDataWriter (String filename) throws IOException {
      this.out = new DataOutputStream(new FileOutputStream(filename));
   }

   /**
    * Truncates any existing file with the specified name.
    *  @param file         file to be created
    */
   public BinaryDataWriter (File file) throws IOException {
      this.out = new DataOutputStream(new FileOutputStream(file));
   }

   /**
    * Constructor.
    *  @param outputStream output stream to write to
    */
   public BinaryDataWriter (OutputStream outputStream) throws IOException {
      this.out = new DataOutputStream(outputStream);
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
      writeLabel(label);
      out.writeByte(TYPECHAR_STRING);
      out.writeByte(0);
      writeStringData(s);
   }

   /**
    * Writes an atomic 32-bit integer (big endian). Writes an anonymous
    * field if `label` is `null`.
    */
   public void write (String label, int a) throws IOException {
      writeLabel(label);
      out.writeByte(TYPECHAR_INTEGER);
      out.writeByte(0);
      out.writeInt(a);
   }

   /**
    * Writes an atomic 32-bit float (big endian). Writes an anonymous
    * field if `label` is `null`.
    */
   public void write (String label, float a) throws IOException {
      writeLabel(label);
      out.writeByte(TYPECHAR_FLOAT);
      out.writeByte(0);
      out.writeFloat(a);
   }

   /**
    * Writes an atomic 64-bit double (big endian). Writes an anonymous
    * field if `label` is `null`.
    */
   public void write (String label, double a) throws IOException {
      writeLabel(label);
      out.writeByte(TYPECHAR_DOUBLE);
      out.writeByte(0);
      out.writeDouble(a);
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
      writeLabel(label);
      out.writeByte(TYPECHAR_STRING);
      out.writeByte(1);
      out.writeInt(n);
      for (int i = 0; i < n; i++)
         writeStringData(a[i]);
   }

   /**
    * Writes the first `n` elements of a one-dimensional array of 32-bit
    * integers (big endian). Writes an anonymous field if `label` is
    * `null`.
    */
   public void write (String label, int[] a, int n) throws IOException {
      writeLabel(label);
      out.writeByte(TYPECHAR_INTEGER);
      out.writeByte(1);
      out.writeInt(n);
      for (int i = 0; i < n; i++)
         out.writeInt(a[i]);
   }

   /**
    * Writes the first `n` elements of a one-dimensional array of 32-bit
    * floats (big endian). Writes an anonymous field if `label` is `null`.
    */
   public void write (String label, float[] a, int n) throws IOException {
      writeLabel(label);
      out.writeByte(TYPECHAR_FLOAT);
      out.writeByte(1);
      out.writeInt(n);
      for (int i = 0; i < n; i++)
         out.writeFloat(a[i]);
   }

   /**
    * Writes the first `n` elements of a one-dimensional array of 64-bit
    * doubles (big endian). Writes an anonymous field if `label` is
    * `null`.
    */
   public void write (String label, double[] a, int n) throws IOException {
      writeLabel(label);
      out.writeByte(TYPECHAR_DOUBLE);
      out.writeByte(1);
      out.writeInt(n);
      for (int i = 0; i < n; i++)
         out.writeDouble(a[i]);
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
      writeLabel(label);
      out.writeByte(TYPECHAR_STRING);
      out.writeByte(2);
      out.writeInt(a.length);
      out.writeInt(a[0].length);
      for (int i = 0; i < a.length; i++)
         for (int j = 0; j < a[i].length; j++)
            writeStringData(a[i][j]);
   }

   /**
    * Writes a two-dimensional array of 32-bit integers (big endian).
    * Writes an anonymous field if `label` is `null`.
    */
   public void write (String label, int[][] a) throws IOException {
      writeLabel(label);
      out.writeByte(TYPECHAR_INTEGER);
      out.writeByte(2);
      out.writeInt(a.length);
      out.writeInt(a[0].length);
      for (int i = 0; i < a.length; i++)
         for (int j = 0; j < a[i].length; j++)
            out.writeInt(a[i][j]);
   }

   /**
    * Writes a two-dimensional array of 32-bit floats (big endian). Writes
    * an anonymous field if `label` is `null`.
    */
   public void write (String label, float[][] a) throws IOException {
      writeLabel(label);
      out.writeByte(TYPECHAR_FLOAT);
      out.writeByte(2);
      out.writeInt(a.length);
      out.writeInt(a[0].length);
      for (int i = 0; i < a.length; i++)
         for (int j = 0; j < a[i].length; j++)
            out.writeFloat(a[i][j]);
   }

   /**
    * Writes a two-dimensional array of 64-bit doubles (big endian).
    * Writes an anonymous field if `label` is `null`.
    */
   public void write (String label, double[][] a) throws IOException {
      writeLabel(label);
      out.writeByte(TYPECHAR_DOUBLE);
      out.writeByte(2);
      out.writeInt(a.length);
      out.writeInt(a[0].length);
      for (int i = 0; i < a.length; i++)
         for (int j = 0; j < a[i].length; j++)
            out.writeDouble(a[i][j]);
   }

   /**
    * @}
    */

   /**
    * @name Other methods
    * @{
    */

   /**
    * Flushes any pending data and closes the file.
    */
   public void close() throws IOException {
      out.close();
   }

}

/**
 * @}
 */