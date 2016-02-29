/*
 * Class:        DataField
 * Description:  Represents a data field
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
import java.lang.reflect.Array;

/**
 * This class represents a data field from a file read by an instance of a
 * class implementing  @ref DataReader.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class DataField {
   protected String label;
   protected Object data;
   protected int effectiveLength;

   /**
    * Constructor. Creates a field named `label` of value `data`.
    *  @param label        name of the field
    *  @param data         value of the field
    */
   public DataField (String label, Object data) {
      this(label, data, -1);
   }

   /**
    * Constructor. Creates a field named `label` of value `data`.
    * `effectiveLength` is the number of significant elements contained in
    * `data` if it is an array.
    *  @param label        name of the field
    *  @param data         value of the field
    *  @param effectiveLength number of significant elements contained in
    *                         `data`
    */
   public DataField (String label, Object data, int effectiveLength) {
      this.label = label;
      this.data = data;
      this.effectiveLength = effectiveLength;
   }

   /**
    * @name Information on the field
    * @{
    */

   /**
    * Returns the field label (or name).
    */
   public String getLabel() {
      return label;
   }

   /**
    * Returns the type of the field.
    */
   public Class getType() {
      return data.getClass();
   }

   /**
    * Returns `true` if the field value is atomic data.
    */
   public boolean isAtomic() {
      return !isArray();
   }

   /**
    * Returns `true` if the field contains an array.
    */
   public boolean isArray() {
      return data.getClass().isArray();
   }

   /**
    * Returns `true` if the field contains a two-dimensional array.
    */
   public boolean isArray2D() {
      return isArray() && Array.get(data, 0).getClass().isArray();
   }

   /**
    * Returns the length of the array contained by the field, or `-1` if
    * it is not an array.
    */
   public int getArrayLength() {
      if (!isArray()) return -1;
      if (effectiveLength < 0) return Array.getLength(data);
      return effectiveLength;
   }

   /**
    * Returns `true` if the field value is an atomic `String`.
    */
   public boolean isString() {
      return (data instanceof String);
   }

   /**
    * Returns `true` if the field value is an atomic `int`.
    */
   public boolean isInt() {
      return (data instanceof Integer);
   }

   /**
    * Returns `true` if the field value is an atomic `float`.
    */
   public boolean isFloat() {
      return (data instanceof Float);
   }

   /**
    * Returns `true` if the field value is an atomic `double`.
    */
   public boolean isDouble() {
      return (data instanceof Double);
   }

   /**
    * @}
    */

   /**
    * @name Obtaining the value as atomic data
    * @{
    */

   /**
    * Returns the value as `String`, or `null` if it is not of type
    * `String`. See  #isString.
    */
   public String asString() {
      return (data instanceof String) ? (String)data : null;
   }

   /**
    * Returns the value as `int` or `0` if it is not of type `int` See
    * #isInt.
    */
   public int asInt() {
      return (data instanceof Integer) ? ((Integer)data).intValue() : 0;
   }

   /**
    * Returns the value as `float` or `0` if it is not of type `float` See
    * #isFloat.
    */
   public float asFloat() {
      return (data instanceof Float) ? ((Float)data).floatValue() : 0;
   }

   /**
    * Returns the value as `double` or `0` if it is not of type `double`
    * See  #isDouble.
    */
   public double asDouble() {
      return (data instanceof Double) ? ((Double)data).doubleValue() : 0;
   }

   /**
    * @}
    */

   /**
    * @name Obtaining the value as a one-dimensional array
    * @{
    */

   /**
    * Returns the value as one-dimensional `String` array or `null` if it
    * is not of type `String[]`.
    */
   public String[] asStringArray() {
      return (data instanceof String[]) ? (String[])data : null;
   }

   /**
    * Returns the value as one-dimensional `int` array or `null` if it is
    * not of type `int[]`.
    */
   public int[] asIntArray() {
      return (data instanceof int[]) ? (int[])data : null;
   }

   /**
    * Returns the value as one-dimensional `float` array or `null` if it
    * is not of type `float[]`.
    */
   public float[] asFloatArray() {
      return (data instanceof float[]) ? (float[])data : null;
   }

   /**
    * Returns the value as one-dimensional `double` array or `null` if it
    * is not of type `double[]`.
    */
   public double[] asDoubleArray() {
      return (data instanceof double[]) ? (double[])data : null;
   }

   /**
    * @}
    */

   /**
    * @name Obtaining the value as a two-dimensional array
    * @{
    */

   /**
    * Returns the value as two-dimensional `String` array or `null` if it
    * is not of type `String[][]`.
    */
   public String[][] asStringArray2D() {
      return (data instanceof String[][]) ? (String[][])data : null;
   }

   /**
    * Returns the value as two-dimensional `int` array or `null` if it is
    * not of type `int[][]`.
    */
   public int[][] asIntArray2D() {
      return (data instanceof int[][]) ? (int[][])data : null;
   }

   /**
    * Returns the value as two-dimensional `float` array or `null` if it
    * is not of type `float[][]`.
    */
   public float[][] asFloatArray2D() {
      return (data instanceof float[][]) ? (float[][])data : null;
   }

   /**
    * Returns the value as two-dimensional `double` array or `null` if it
    * is not of type `double[][]`.
    */
   public double[][] asDoubleArray2D() {
      return (data instanceof double[][]) ? (double[][])data : null;
   }

   /**
    * @}
    */

   /**
    * @name Obtaining the value as an `Object`
    * @{
    */

   /**
    * Returns the value of the field as an `Object`.
    */
   public Object asObject() {
      return data;
   }

}

/**
 * @}
 */