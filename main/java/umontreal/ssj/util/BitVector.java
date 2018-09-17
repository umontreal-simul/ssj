/*
 * Class:        BitVector
 * Description:  implements vectors of bits and their operations
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since
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
package umontreal.ssj.util;

import java.io.Serializable;

/**
 * This class implements vectors of bits and the operations needed to use
 * them. The vectors can be of arbitrary length. The operations provided are
 * all the binary operations available to the `int` and `long` primitive
 * types in Java.
 *
 * All bit operations are present in two forms: a normal form and a `self`
 * form. The normal form returns a newly created object containing the
 * result, while the `self` form puts the result in the calling object
 * (<tt>this</tt>). The return value of the `self` form is the calling object
 * itself. This is done to allow easier manipulation of the results, making
 * it possible to chain operations.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class BitVector implements Serializable, Cloneable {

   static final long serialVersionUID = -3448233092524725148L;

   private int[] v;       //the bits data
   private int length;    //number of data bits (in bits, not in bytes)

   private final static int all_1 = -1;  //integer with all bits set to 1
   private final static int one_1 = 1;   //integer with only his last bit set to 1
   /*
     Note sur le format interne du vecteur de bits :
     On fait toujours en sorte que les bits redondants (ceux qui apparaissent
     quand length % 32 != 0) soient mis a 0. Ceci permet de faire des 
     operations entre des vecteurs de longeurs differentes en posant que
     les bits manquants sur le plus petit des deux vecteurs ont la valeur 0.
   */

   /**
    * Creates a new `BitVector` of length `length` with all its bits set
    * to 0.
    *  @param length       the length of the `BitVector`
    */
   public BitVector (int length) {
      this.length = length;
      v = new int[(length + 31) / 32];
      for(int i = 0; i < v.length; i++)
         v[i] = 0;
   }

   /**
    * Creates a new `BitVector` of length `length` using the data in
    * `vect`. Component `vect[0]` makes the 32 lowest order bits, with
    * `vect[1]` being the 32 next lowest order bits, and so on. The normal
    * bit order is then used to fill the 32 bits (the first bit is the
    * lowest order bit and the last bit is largest order bit). Note that
    * the sign bit is used as the largest order bit.
    *  @param vect         the bits data
    *  @param length       the length of the vector
    *  @exception IllegalArgumentException when the length of `vect` is
    * not compatible with the `length` provided
    */
   public BitVector (int[] vect, int length) {
      if (((length + 31)/ 32) != vect.length)
         throw new IllegalArgumentException
         ("The int[] length must be equal to the (length + 31) / 32");

      this.length = length;
      v = new int[vect.length];
      for(int i = 0; i < vect.length; i++)
         v[i] = vect[i];

      //the extra bits must be set at zero
      v[v.length - 1] &= (all_1 >>> (31 - (length - 1) % 32));
   }

   /**
    * Creates a new `BitVector` using the data in `vect`. The length of
    * the `BitVector` is always equals to 32 times the length of `vect`.
    *  @param vect         the bits data
    */
   public BitVector (int[] vect) {
      this(vect, vect.length * 32);
   }

   /**
    * Creates a copy of the `BitVector that`.
    *  @param that         the `BitVector` to copy
    */
   public BitVector (BitVector that) {
      this.length = that.length;
      this.v = new int[that.v.length];
      for(int i = 0; i < that.v.length; i++)
         this.v[i] = that.v[i];
   }

   /**
    * Creates a copy of the `BitVector`.
    *  @return a deep copy of the `BitVector`
    */
   public Object clone() {
      try{
         BitVector c = (BitVector)super.clone();
         c.v = (int[])v.clone();
         return c;
      }catch(CloneNotSupportedException e) {
         IllegalStateException ne = new IllegalStateException();
         ne.initCause(e);
         throw ne;
      }      
   }

   /**
    * Verifies if two <tt>BitVector</tt>’s have the same length and the
    * same data.
    *  @param that         the other `BitVector` to compare to
    *  @return if the two <tt>BitVector</tt>’s are identiqual
    */
   public boolean equals (BitVector that) {
      if(this.length != that.length)
         return false;
      for(int i = 0; i < this.v.length; i++)
         if(this.v[i] != that.v[i])
            return false;
      return true;
   }

   /**
    * Returns the length of the `BitVector`.
    *  @return the length of the `BitVector`
    */
   public int size() {
      return length;
   }

   /**
    * Resizes the `BitVector` so that its length is equal to `size`. If
    * the `BitVector` is enlarged, then the newly added bits are given the
    * value 1 if `filling` is set to `true` and 0 otherwise.
    *  @param size         the new size of the `BitVector`
    *  @param filling      the state of the new bits
    */
   public void enlarge (int size, boolean filling) {
      if(size < 0)
         throw new NegativeArraySizeException
         ("The BitVector must have a non-negative size");
      if(filling && (length % 32) != 0)
         v[v.length - 1] ^= all_1 << (length % 32);
      if((size + 31) / 32 != v.length) {
         int[] new_v = new int[(size + 31) / 32];
         int i;
         for(i = 0; i < new_v.length && i < v.length; i++)
            new_v[i] = v[i];
         for(; i < new_v.length; i++)
            new_v[i] = filling ? all_1 : 0;
         v = new_v;
      }
      length = size;

      //the extra bits must be set at zero
      v[v.length - 1] &= (all_1 >>> (31 - (length - 1) % 32));
   }

   /**
    * Resizes the `BitVector` so that its length is equal to `size`. Any
    * new bit added is set to 0.
    *  @param size         the new size of the `BitVector`
    */
   public void enlarge (int size) {
      enlarge(size, false);
   }

   /**
    * Gives the value of the bit in position `pos`. If the value is 1,
    * returns `true`; otherwise, returns `false`.
    *  @param pos          the position of the checked bit
    *  @return the value of the bit as a boolean
    *
    *  @exception ArrayIndexOutOfBoundsException if `pos` is outside the
    * range of the `BitVector`
    */
   public boolean getBool (int pos) {
      if(pos < 0 || pos >= length)
         throw new ArrayIndexOutOfBoundsException(pos);
      return (v[pos >>> 5] & (one_1 << (pos & 31))) != 0;
   }

   /**
    * Sets the value of the bit in position `pos`. If `value` is equal to
    * `true`, sets it to 1; otherwise, sets it to 0.
    *  @param pos          the position of the bit to modify
    *  @param value        the new value of the bit as a boolean
    *  @exception ArrayIndexOutOfBoundsException if `pos` is outside the
    * range of the `BitVector`
    */
   public void setBool (int pos, boolean value) {
      if(pos > length || pos < 0)
         throw new ArrayIndexOutOfBoundsException(pos);
      if(value)
         v[pos / 32] |= (one_1 << (pos % 32));
      else
         v[pos / 32] &= (one_1 << (pos % 32)) ^ all_1;
   }

   /**
    * Returns an `int` containing all the bits in the interval
    * @f$[\mathtt{pos} \times32, \mathtt{pos} \times32 + 31]@f$.
    *  @param pos          the selected position
    *  @return the int at the specified position
    *
    *  @exception ArrayIndexOutOfBoundsException if `pos` is outside the
    * range of the `BitVector`
    */
   public int getInt (int pos) {
      if(pos >= v.length || pos < 0)
         throw new ArrayIndexOutOfBoundsException(pos);
      return v[pos];
   }

   /**
    * Returns a string containing all the bits of the `BitVector`,
    * starting with the highest order bit and finishing with the lowest
    * order bit. The bits are grouped by groups of 8 bits for ease of
    * reading.
    *  @return all the bits of the `BitVector`
    */
   public String toString() {
      StringBuffer sb = new StringBuffer();

      for(int i = length - 1; i > 0; i--)
         sb.append(getBool(i) ? "1" : "0").append(i%8 == 0 ? " " : "");
      sb.append(getBool(0) ? "1" : "0");

      return sb.toString();
   }

   /**
    * Returns a `BitVector` which is the result of the `not` operator on
    * the current `BitVector`. The `not` operator is equivalent to the
    * <code>~</code> operator in Java and thus swap all bits (bits
    * previously set to 0 become 1 and bits previously set to 1 become 0).
    *  @return the effect of the `not` operator
    */
   public BitVector not() {
      BitVector bv = new BitVector(length);
      for(int i = 0; i < v.length; i++)
         bv.v[i] = v[i] ^ all_1;

      //the extra bits must be set at zero
      bv.v[v.length - 1] &= (all_1 >>> (31 - (length - 1) % 32));

      return bv;
   }

   /**
    * Applies the `not` operator on the current `BitVector` and returns
    * it.
    *  @return the `BitVector` itself
    */
   public BitVector selfNot() {
      for(int i = 0; i < v.length; i++)
         v[i] = v[i] ^ all_1;

      //the extra bits must be set at zero
      v[v.length - 1] &= (all_1 >>> (31 - (length - 1) % 32));

      return this;
   }

   /**
    * Returns a `BitVector` which is the result of the `xor` operator
    * applied on `this` and `that`. The `xor` operator is equivalent to
    * the <code>^</code> operator in Java. All bits which were set to 0 in
    * one of the vector and to 1 in the other vector are set to 1. The
    * others are set to 0. This is equivalent to the addition in modulo 2
    * arithmetic.
    *  @param that         the second operand to the `xor` operator
    *  @return the result of the `xor` operation
    */
   public BitVector xor (BitVector that) {
      //we always consider that this is longer than that
      if(that.length > this.length)
         return that.xor(this);

      BitVector bv = new BitVector(this.length);

      int max = this.v.length;
      int min = that.v.length;

      for(int i = 0; i < min; i++)
         bv.v[i] = this.v[i] ^ that.v[i];
      for(int i = min; i < max; i++)
         bv.v[i] = this.v[i];

      return bv;
   }

   /**
    * Applies the `xor` operator on `this` with `that`. Stores the result
    * in `this` and returns it.
    *  @param that         the second operand to the `xor` operator
    *  @return `this`
    */
   public BitVector selfXor (BitVector that) {
      //we assume that this is large enough to contain that
      if(this.length < that.length)
         this.enlarge(that.length);

      int min = that.v.length;

      for(int i = 0; i < min; i++)
         this.v[i] ^= that.v[i];

      return this;
   }

   /**
    * Returns a `BitVector` which is the result of the `and` operator with
    * both the `this` and `that` <tt>BitVector</tt>’s. The `and` operator
    * is equivalent to the <code>&</code> operator in Java. Only bits
    * which are set to 1 in both `this` and `that` are set to 1 in the
    * result, all the others are set to 0.
    *  @param that         the second operand to the `and` operator
    *  @return the result of the `and` operation
    */
   public BitVector and (BitVector that) {
      //we always consider this as longer than that
      if(that.length > this.length)
         return that.and(this);

      BitVector bv = new BitVector(this.length);

      int max = this.v.length;
      int min = that.v.length;

      for(int i = 0; i < min; i++)
         bv.v[i] = this.v[i] & that.v[i];
      if(that.length % 32 != 0)
         bv.v[min - 1] |= this.v[min - 1] & (all_1 << (that.length % 32));
      for(int i = min; i < max; i++)
         bv.v[i] = 0;

      return bv;
   }

   /**
    * Applies the `and` operator on `this` with `that`. Stores the result
    * in `this` and returns it.
    *  @param that         the second operand to the `and` operator
    *  @return `this`
    */
   public BitVector selfAnd (BitVector that) {
      //we assume that this is large enough to contain that
      if(this.length < that.length)
         this.enlarge(that.length, true);

      int min = that.v.length;

      for(int i = 0; i < min - 1; i++)
         this.v[i] &= that.v[i];
      this.v[min - 1] &= (that.v[min - 1] | (all_1 << (that.length % 32)));

      return this;
   }

   /**
    * Returns a `BitVector` which is the result of the `or` operator with
    * both the `this` and `that` <tt>BitVector</tt>’s. The `or` operator
    * is equivalent to the `|` operator in Java. Only bits which are set
    * to 0 in both `this` and `that` are set to to 0 in the result, all
    * the others are set to 1.
    *  @param that         the second operand to the `or` operator
    *  @return the result of the `or` operation
    */
   public BitVector or (BitVector that) {
      //we always consider this is longer than that
      if(that.length > this.length)
         return that.or(this);

      BitVector bv = new BitVector(this.length);

      int max = this.v.length;
      int min = that.v.length;

      for(int i = 0; i < min; i++)
         bv.v[i] = this.v[i] | that.v[i];
      for(int i = min; i < max; i++)
         bv.v[i] = 0;

      return bv;
   }

   /**
    * Applies the `or` operator on `this` with `that`. Stores the result
    * in `this` and returns it.
    *  @param that         the second operand to the `or` operator
    *  @return `this`
    */
   public BitVector selfOr (BitVector that) {
      //we assume that this is large enough to contain that
      if(this.length < that.length)
         this.enlarge(that.length);

      int min = that.v.length;

      for(int i = 0; i < min; i++)
         this.v[i] |= that.v[i];

      return this;
   }

   /**
    * Returns a `BitVector` equal to the original with all the bits
    * shifted `j` positions to the right if `j` is positive, and shifted
    * `j` positions to the left if `j` is negative. The new bits that
    * appears to the left or to the right are set to 0. If `j` is
    * positive, this operation is equivalent to the `>>>` operator in
    * Java, otherwise, it is equivalent to the `<<` operator.
    *  @param j            the size of the shift
    *  @return the shifted `BitVector`
    */
   public BitVector shift (int j) {
      BitVector bv = new BitVector(length);

      if(j == 0)
         return bv;
      else if(j > 0) {
         int a = j / 32;
         int b = j % 32;
         int c = 32 - b;

         int i;
         for(i = 0; i+a < v.length; i++)
            bv.v[i] = v[i+a] >>> b;
         // la retenue
         for(i = 0; i+a+1 < v.length; i++)
            bv.v[i] ^= v[i+a+1] << c;
      } else // if(j < 0)
      {
         j = -j;
         int a = j / 32;
         int b = j % 32;
         int c = 32 - b;

         int i;
         for(i = a; i < v.length; i++)
            bv.v[i] ^= v[i-a] << b;
         // la retenue
         for(i = a+1; i < v.length; i++)
            bv.v[i] ^= v[i-a-1] >>> c;
      }

      return bv;
   }

   /**
    * Shift all the bits of the current `BitVector` `j` positions to the
    * right if `j` is positive, and `j` positions to the left if `j` is
    * negative. The new bits that appears to the left or to the rigth are
    * set to 0. Returns `this`.
    *  @param j            the size of the shift
    *  @return `this`
    */
   public BitVector selfShift (int j) {
      if(j == 0)
         return this;
      else if(j > length || j < -length) {
         for(int i = 0; i < v.length; i++)
            v[i] = 0;
      } else if(j > 0) {
         int a = j / 32;
         int b = j % 32;
         int c = 32 - b;

         int i;
         for(i = 0; i+a+1 < v.length; i++) {
            v[i] = v[i+a] >>> b;
            // la retenue
            v[i] ^= v[i+a+1] << c;
         }
         v[i] = v[i+a] >>> b;
         for(i += 1; i < v.length; i++)
            v[i] = 0;
      } else // if(j < 0)
      {
         j = -j;
         int a = j / 32;
         int b = j % 32;
         int c = 32 - b;

         int i;
         for(i = v.length - 1; i > a; i--) {
            v[i] = v[i-a] << b;
            // la retenue
            v[i] ^= v[i-a-1] >>> c;
         }
         v[i] = v[i-a] << b;
         for(i -= 1; i >= 0; i--)
            v[i] = 0;
      }

      return this;
   }

   /**
    * Returns the scalar product of two <tt>BitVector</tt>’s modulo 2. It
    * returns `true` if there is an odd number of bits with a value of 1
    * in the result of the `and` operator applied on `this` and `that`,
    * and returns `false` otherwise.
    *  @param that         the other `BitVector` with which to do the
    *                      scalar product
    *  @return the scalar product
    */
   public boolean scalarProduct (BitVector that) {
      //we must take that is not longer than this
      if(that.v.length > this.v.length)
         return that.scalarProduct(this);

      boolean result = false;
      int prod;

      for(int i = 0; i < that.v.length; i++) {
         prod = this.v[i] & that.v[i];
         while(prod != 0) {
            // a chaque iteration, on enleve le 1 le plus a droite
            prod &= prod - 1;
            result = !result;
         }
      }

      return result;
   }

}