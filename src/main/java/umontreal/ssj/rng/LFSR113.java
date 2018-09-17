/*
 * Class:        LFSR113
 * Description:  32-bit composite linear feedback shift register proposed by L'Ecuyer
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
package umontreal.ssj.rng;

import java.io.Serializable;

/**
 * Extends  @ref RandomStreamBase using a composite linear feedback shift
 * register (LFSR) (or Tausworthe) RNG as defined in @cite rLEC96a,
 * @cite rTEZ91b&thinsp;. This generator is the `LFSR113` proposed by
 * @cite rLEC99a&thinsp;. It has four 32-bit components combined by a bitwise
 * xor. Its period length is  @f$\rho\approx2^{113}@f$. The values of
 * @f$V@f$, @f$W@f$ and @f$Z@f$ are @f$2^{35}@f$, @f$2^{55}@f$ and
 * @f$2^{90}@f$ respectively (see  @ref RandomStream for their definition).
 * The seed of the RNG, and the state of a stream at any given step, are
 * four-dimensional vectors of 32-bit integers. The default initial seed of
 * the RNG is <sup title="In previous versions, it was @f$(12345, 12345,
 * 12345, 12345)@f$.">[1]</sup> (987654321, 987654321, 987654321, 987654321).
 * The `nextValue` method returns numbers with 32 bits of precision.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class LFSR113 extends RandomStreamBase {

   private static final long serialVersionUID = 70510L;
   // La date de modification a l'envers, lire 10/05/2007

   // generator constant: make sure that double values 0 and 1 never occur
   private static final double NORM = 1.0 / 0x100000001L;   // 2^32 + 1


   // state variables:
   private int z0;
   private int z1;
   private int z2;
   private int z3;

   //stream and substream variables :
   private int[] stream;
   private int[] substream;
   private static int[] curr_stream = {987654321, 987654321, 987654321, 987654321};

   /**
    * Constructs a new stream.
    */
   public LFSR113() {
      name = null;

      stream = new int[4];
      substream = new int[4];

      for(int i = 0; i < 4; i++)
         stream[i] = curr_stream[i];

      resetStartStream();


      // Les operations qui suivent permettent de faire sauter en avant
      // de 2^90 iterations chacunes des composantes du generateur.
      // L'etat interne apres le saut est cependant legerement different
      // de celui apres 2^90 iterations puisqu'il ignore l'etat dans
      // lequel se retrouvent les premiers bits de chaque composantes,
      // puisqu'ils sont ignores dans la recurrence. L'etat redevient
      // identique a ce que l'on aurait avec des iterations normales
      // apres un appel a nextValue().

      int z, b;

      z = curr_stream[0] & -2;
      b = (z <<  6) ^ z;
      z = (z) ^ (z << 2) ^ (z << 3) ^ (z << 10) ^ (z << 13) ^
         (z << 16) ^ (z << 19) ^ (z << 22) ^ (z << 25) ^
         (z << 27) ^ (z << 28) ^
         (b >>> 3) ^ (b >>> 4) ^ (b >>> 6) ^ (b >>> 9) ^ (b >>> 12) ^
         (b >>> 15) ^ (b >>> 18) ^ (b >>> 21);
      curr_stream[0] = z;


      z = curr_stream[1] & -8;
      b = (z <<  2) ^ z;
      z = (b >>> 13) ^ (z << 16);
      curr_stream[1] = z;


      z = curr_stream[2] & -16;
      b = (z <<  13) ^ z;
      z = (z << 2) ^ (z << 4) ^ (z << 10) ^ (z << 12) ^ (z << 13) ^
         (z << 17) ^ (z << 25) ^
         (b >>> 3) ^ (b >>> 11) ^ (b >>> 15) ^ (b >>> 16) ^ (b >>> 24);
      curr_stream[2] = z;


      z = curr_stream[3] & -128;
      b = (z <<  3) ^ z;
      z = (z << 9) ^ (z << 10) ^ (z << 11) ^ (z << 14) ^ (z << 16) ^
         (z << 18) ^ (z << 23) ^ (z << 24) ^
         (b >>> 1) ^ (b >>> 2) ^ (b >>> 7) ^ (b >>> 9) ^ (b >>> 11) ^
         (b >>> 14) ^ (b >>> 15) ^ (b >>> 16) ^ (b >>> 23) ^ (b >>> 24);
      curr_stream[3] = z;
   }

   /**
    * Constructs a new stream with the identifier `name`.
    *  @param name         name of the stream
    */
   public LFSR113 (String name) {
      this();
      this.name = name;
   }

   /**
    * Sets the initial seed for the class `LFSR113` to the four integers
    * of the vector `seed[0..3]`. This will be the initial state of the
    * next created stream. The default seed for the first stream is <sup
    * title="In previous versions, it was @f$(12345, 12345, 12345,
    * 12345)@f$.">[2]</sup> (987654321, 987654321, 987654321, 987654321).
    * The first, second, third and fourth integers of `seed` must be
    * either negative, or greater than or equal to 2, 8, 16 and 128
    * respectively.
    *  @param seed         array of 4 elements representing the seed
    */
   public static void setPackageSeed (int[] seed) {
      checkSeed (seed);
      for(int i = 0; i < 4; i++)
         curr_stream[i] = seed[i];
   }

   private static void checkSeed  (int[] seed) {
      if (seed.length < 4)
         throw new IllegalArgumentException("Seed must contain 4 values");
      if ((seed[0] >= 0 && seed[0] < 2)  ||
          (seed[1] >= 0 && seed[1] < 8)  ||
          (seed[2] >= 0 && seed[2] < 16) ||
          (seed[3] >= 0 && seed[3] < 128))
         throw new IllegalArgumentException
         ("The seed elements must be either negative or greater than 1, 7, 15 and 127 respectively");
   }

   /**
    * This method is discouraged for normal use. Initializes the stream at
    * the beginning of a stream with the initial seed `seed[0..3]`. The
    * seed must satisfy the same conditions as in `setPackageSeed`. This
    * method only affects the specified stream; the others are not
    * modified, so the beginning of the streams will not be spaced @f$Z@f$
    * values apart. For this reason, this method should only be used in
    * very exceptional cases; proper use of the `reset...` methods and of
    * the stream constructor is preferable.
    *  @param seed         array of 4 elements representing the seed
    */
   public void setSeed (int[] seed) {
      checkSeed  (seed);
      for(int i = 0; i < 4; i++)
         stream[i] = seed[i];
      resetStartStream();
   }

   /**
    * Returns the current state of the stream, represented as an array of
    * four integers.
    *  @return the current state of the stream
    */
   public int[] getState() {
      return new int[]{z0, z1, z2, z3};
   }

   /**
    * Clones the current generator and return its copy.
    *  @return A deep copy of the current generator
    */
   public LFSR113 clone() {
      LFSR113 retour = null;
      retour = (LFSR113)super.clone();
      retour.stream = new int[4];
      retour.substream = new int[4];
      for (int i = 0; i<4; i++) {
         retour.substream[i] = substream[i];
         retour.stream[i] = stream[i];
      }
      return retour;
   }


   public void resetStartStream() {
      for(int i = 0; i < 4; i++)
         substream[i] = stream[i];
      resetStartSubstream();
   }

   public void resetStartSubstream() {
      z0 = substream[0];
      z1 = substream[1];
      z2 = substream[2];
      z3 = substream[3];
   }
/*
   // La version de Mario: beaucoup plus lent que l'ancienne version: on garde
   // l'ancienne version.
   public void resetNextSubstream() {
      byte [] c0 = new byte[] {0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
                            0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1};
      byte [] c1 = new byte[] {1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                            0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0};
      byte [] c2 = new byte[] {0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1,
                            0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0};
      byte [] c3 = new byte[] {1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0,
                            0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      int x0 = 0;
      int x1 = 0;
      int x2 = 0;
      int x3 = 0;
      resetStartSubstream();
      for (int i = 0; i < 31; ++i) {
         if (c0[i] == 1) x0 ^= z0;
         if (c1[i] == 1) x1 ^= z1;
         if (c2[i] == 1) x2 ^= z2;
	 if (c3[i] == 1) x3 ^= z3;

         int b;
         b  = (((z0 <<   6) ^ z0) >>> 13);
         z0 = (((z0 &   -2) << 18) ^ b);
         b  = (((z1 <<   2) ^ z1) >>> 27);
         z1 = (((z1 &   -8) <<  2) ^ b);
         b  = (((z2 <<  13) ^ z2) >>> 21);
         z2 = (((z2 &  -16) <<  7) ^ b);
         b  = (((z3 <<   3) ^ z3) >>> 12);
         z3 = (((z3 & -128) << 13) ^ b);
      }
      substream[0] = x0;
      substream[1] = x1;
      substream[2] = x2;
      substream[3] = x3;
      resetStartSubstream();
   }
*/

   public void resetNextSubstream() {
      // Les operations qui suivent permettent de faire sauter en avant
      // de 2^55 iterations chacunes des composantes du generateur.
      // L'etat interne apres le saut est cependant legerement different
      // de celui apres 2^55 iterations puisqu'il ignore l'etat dans
      // lequel se retrouvent les premiers bits de chaque composantes,
      // puisqu'ils sont ignores dans la recurrence. L'etat redevient
      // identique a ce que l'on aurait avec des iterations normales
      // apres un appel a nextValue().

      int z, b;

      z = substream[0] & -2;
      b = (z <<  6) ^ z;
      z = (z) ^ (z << 3) ^ (z << 4) ^ (z << 6) ^ (z << 7) ^
         (z << 8) ^ (z << 10) ^ (z << 11) ^ (z << 13) ^ (z << 14) ^
         (z << 16) ^ (z << 17) ^ (z << 18) ^ (z << 22) ^
         (z << 24) ^ (z << 25) ^ (z << 26) ^ (z << 28) ^ (z << 30);
      z ^= (b >>> 1) ^ (b >>> 3) ^ (b >>> 5) ^ (b >>> 6) ^
         (b >>> 7) ^ (b >>> 9) ^ (b >>> 13) ^ (b >>> 14) ^
         (b >>> 15) ^ (b >>> 17) ^ (b >>> 18) ^ (b >>> 20) ^
         (b >>> 21) ^ (b >>> 23) ^ (b >>> 24) ^ (b >>> 25) ^
         (b >>> 26) ^ (b >>> 27) ^ (b >>> 30);
      substream[0] = z;


      z = substream[1] & -8;
      b = z ^ (z << 1);
      b ^= (b << 2);
      b ^= (b << 4);
      b ^= (b << 8);

      b <<= 8;
      b ^= (z << 22) ^ (z << 25) ^ (z << 27);
      if((z & 0x80000000) != 0) b ^= 0xABFFF000;
      if((z & 0x40000000) != 0) b ^= 0x55FFF800;
      z = b ^ (z >>> 7) ^ (z >>> 20) ^ (z >>> 21);
      substream[1] = z;


      z = substream[2] & -16;
      b = (z <<  13) ^ z;
      z = (b >>> 3) ^ (b >>> 17) ^
         (z << 10) ^ (z << 11) ^ (z << 25);
      substream[2] = z;


      z = substream[3] & -128;
      b = (z <<  3) ^ z;
      z = (z << 14) ^ (z << 16) ^ (z << 20) ^
         (b >>> 5) ^ (b >>> 9) ^ (b >>> 11);
      substream[3] = z;

      resetStartSubstream();
   }


   public String toString()  {
      if (name == null)
         return "The state of the LFSR113 is: { " +
                z0 + ", " + z1 + ", " + z2 + ", " + z3 + " }";
      else
         return "The state of " + name + " is: { " +
                z0 + ", " + z1 + ", " + z2 + ", " + z3 + " }";
   }

   private long nextNumber() {
      int b;
      b  = (((z0 <<   6) ^ z0) >>> 13);
      z0 = (((z0 &   -2) << 18) ^ b);
      b  = (((z1 <<   2) ^ z1) >>> 27);
      z1 = (((z1 &   -8) <<  2) ^ b);
      b  = (((z2 <<  13) ^ z2) >>> 21);
      z2 = (((z2 &  -16) <<  7) ^ b);
      b  = (((z3 <<   3) ^ z3) >>> 12);
      z3 = (((z3 & -128) << 13) ^ b);

      long r = (z0 ^ z1 ^ z2 ^ z3);

      if (r <= 0)
         r += 0x100000000L;      //2^32

      return r;
   }

   protected double nextValue() {
      // Make sure that double values 0 and 1 never occur
      return nextNumber() * NORM;
   }

   public int nextInt (int i, int j) {
      if (i > j)
         throw new IllegalArgumentException(i + " is larger than " + j + ".");
      long d = j-i+1L;
      long q = 0x100000000L / d;
      long r = 0x100000000L % d;
      long res;

      do {
         res = nextNumber();
      } while (res >= 0x100000000L - r);

      return (int) (res / q) + i;
   }

}