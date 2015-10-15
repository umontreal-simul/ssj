/*
 * Class:        DataReader
 * Description:  Data reader interface
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
import java.util.Map;

/**
 * Data reader interface.
 *
 * <div class="SSJ-bigskip"></div>
 */
public interface DataReader {

   /**
    * @name Reading atomic data
    * @{
    */

   /**
    * Reads the first field labeled as `label` and returns its `String`
    * value.
    */
   public String readString (String label) throws IOException;

   /**
    * Reads the first field labeled as `label` and returns its `int`
    * value.
    */
   public int readInt (String label) throws IOException;

   /**
    * Reads the first field labeled as `label` and returns its `float`
    * value.
    */
   public float readFloat (String label) throws IOException;

   /**
    * Reads the first field labeled as `label` and returns its `double`
    * value.
    */
   public double readDouble (String label) throws IOException;

   /**
    * @}
    */

   /**
    * @name Reading one-dimensional arrays
    * @{
    */

   /**
    * Reads the first field labeled as `label` and returns its value as a
    * one-dimensional array of <tt>String</tt>’s.
    */
   public String[] readStringArray (String label) throws IOException;

   /**
    * Reads the first field labeled as `label` and returns its value as a
    * one-dimensional array of <tt>int</tt>’s.
    */
   public int[] readIntArray (String label) throws IOException;

   /**
    * Reads the first field labeled as `label` and returns its value as a
    * one-dimensional array of <tt>float</tt>’s.
    */
   public float[] readFloatArray (String label) throws IOException;

   /**
    * Reads the first field labeled as `label` and returns its value as a
    * one-dimensional array of <tt>double</tt>’s.
    */
   public double[] readDoubleArray (String label) throws IOException;

   /**
    * @}
    */

   /**
    * @name Reading two-dimensional arrays
    * @{
    */

   /**
    * Reads the first field labeled as `label` and returns its value as a
    * two-dimensional array of <tt>String</tt>’s.
    */
   public String[][] readStringArray2D (String label) throws IOException;

   /**
    * Reads the first field labeled as `label` and returns its value as a
    * two-dimensional array of <tt>int</tt>’s.
    */
   public int[][] readIntArray2D (String label) throws IOException;

   /**
    * Reads the first field labeled as `label` and returns its value as a
    * two-dimensional array of <tt>float</tt>’s.
    */
   public float[][] readFloatArray2D (String label) throws IOException;

   /**
    * Reads the first field labeled as `label` and returns its value as a
    * two-dimensional array of <tt>double</tt>’s.
    */
   public double[][] readDoubleArray2D (String label) throws IOException;

   /**
    * @}
    */

   /**
    * @name Reading fields of unknown type
    * @{
    */

   /**
    * Reads all remaining fields in the file and returns a hashmap indexed
    * by field labels. Anonymous fields are mapped to
    * <code>"_data01_"</code>, <code>"_data02_"</code>, …
    */
   public Map<String, DataField> readAllNextFields() throws IOException;

   /**
    * Reads all fields in the file and returns a hashmap indexed by field
    * labels. Anonymous fields are mapped to <code>"_data01_"</code>,
    * <code>"_data02_"</code>, …
    */
   public Map<String, DataField> readAllFields() throws IOException;

   /**
    * Reads the next available field.
    *  @return a newly created DataField instance or `null` if not found
    */
   public DataField readNextField() throws IOException;

   /**
    * Reads the first field labeled as `label`.
    *  @return a newly created DataField instance or `null` if not found
    */
   public DataField readField (String label) throws IOException;

   /**
    * @}
    */

   /**
    * @name Other methods
    * @{
    */

   /**
    * Closes the input stream.
    */
   public void close() throws IOException;

   /**
    * Resets the reader to its initial state, i.e. goes back to the
    * beginning of the data stream, if possible.
    */
   public void reset() throws IOException;

   /**
    * Returns `true` if there remains data to be read.
    */
   public boolean dataPending() throws IOException;

}

/**
 * @}
 */