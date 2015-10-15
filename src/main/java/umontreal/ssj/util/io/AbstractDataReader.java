/*
 * Class:        AbstractDataReader
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
import java.util.Map;
import java.util.HashMap;

/**
 * This abstract class implements shared functionality for data readers.
 *
 * <div class="SSJ-bigskip"></div>
 */
public abstract class AbstractDataReader implements DataReader {

   /**
    * @name Reading atomic data
    * @{
    */

   /**
    * Reads first field labeled as `label` and returns its `String` value.
    */
   public String readString (String label) throws IOException {
      return readField(label).asString();
   }

   /**
    * Reads first field labeled as `label` and returns its `int` value.
    */
   public int readInt (String label) throws IOException {
      return readField(label).asInt();
   }

   /**
    * Reads first field labeled as `label` and returns its `float` value.
    */
   public float readFloat (String label) throws IOException {
      return readField(label).asFloat();
   }

   /**
    * Reads first field labeled as `label` and returns its `double` value.
    */
   public double readDouble (String label) throws IOException {
      return readField(label).asDouble();
   }

   /**
    * @}
    */

   /**
    * @name Reading one-dimensional arrays
    * @{
    */

   /**
    * Reads first field labeled as `label` and returns its value as a
    * one-dimensional array of <tt>String</tt>’s.
    */
   public String[] readStringArray (String label) throws IOException {
      return readField(label).asStringArray();
   }

   /**
    * Reads first field labeled as `label` and returns its value as a
    * one-dimensional array of <tt>int</tt>’s.
    */
   public int[] readIntArray (String label) throws IOException {
      return readField(label).asIntArray();
   }

   /**
    * Reads first field labeled as `label` and returns its value as a
    * one-dimensional array of <tt>float</tt>’s.
    */
   public float[] readFloatArray (String label) throws IOException {
      return readField(label).asFloatArray();
   }

   /**
    * Reads first field labeled as `label` and returns its value as a
    * one-dimensional array of <tt>double</tt>’s.
    */
   public double[] readDoubleArray (String label) throws IOException {
      return readField(label).asDoubleArray();
   }

   /**
    * @}
    */

   /**
    * @name Reading two-dimensional arrays
    * @{
    */

   /**
    * Reads first field labeled as `label` and returns its value as a
    * two-dimensional array of <tt>String</tt>’s.
    */
   public String[][] readStringArray2D (String label) throws IOException {
      return readField(label).asStringArray2D();
   }

   /**
    * Reads first field labeled as `label` and returns its value as a
    * two-dimensional array of <tt>int</tt>’s.
    */
   public int[][] readIntArray2D (String label) throws IOException {
      return readField(label).asIntArray2D();
   }

   /**
    * Reads first field labeled as `label` and returns its value as a
    * two-dimensional array of <tt>float</tt>’s.
    */
   public float[][] readFloatArray2D (String label) throws IOException {
      return readField(label).asFloatArray2D();
   }

   /**
    * Reads first field labeled as `label` and returns its value as a
    * two-dimensional array of <tt>double</tt>’s.
    */
   public double[][] readDoubleArray2D (String label) throws IOException {
      return readField(label).asDoubleArray2D();
   }

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
   public Map<String, DataField> readAllNextFields() throws IOException {

      HashMap<String,DataField> fields = new HashMap<String,DataField>();
      
      int iAnonymous = 0;
      
      while (dataPending()) {

         DataField data = readNextField();

         String key = data.getLabel();
         if (key == null)
            key = String.format("_data%02d_", ++iAnonymous);
         fields.put(key, data);

      }
      
      return fields;
   }

   /**
    * Reads all fields in the file and returns a hashmap indexed by field
    * labels. Anonymous fields are mapped to <code>"_data01_"</code>,
    * <code>"_data02_"</code>, …
    */
   public Map<String, DataField> readAllFields() throws IOException {
      reset();
      return readAllNextFields();
   }

}

/**
 * @}
 */