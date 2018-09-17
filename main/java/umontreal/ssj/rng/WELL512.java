/*
 * Class:        WELL512
 * Description:  a Well Equidistributed Long period Linear Random Number
                 Generator with a state size of 512 bits
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
 * This class implements the  @ref RandomStream interface via inheritance
 * from  @ref RandomStreamBase. The backbone generator is a Well
 * Equidistributed Long period Linear Random Number Generator (WELL),
 * proposed by F. Panneton in @cite rPAN06b, @cite rPAN04t&thinsp;, and which
 * has a state size of 512 bits and a period length of  @f$\rho\approx@f$
 * @f$2^{512}@f$. The values of @f$V@f$, @f$W@f$ and @f$Z@f$ are
 * @f$2^{150}@f$, @f$2^{200}@f$ and @f$2^{350}@f$ respectively (see
 * @ref RandomStream for their definition). The seed of the RNG, and the
 * state of a stream at any given step, is a 16-dimensional vector of 32-bit
 * integers.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class WELL512 extends RandomStreamBase {

   private static final long serialVersionUID = 120307L;
   //La date de modification a l'envers, lire 07/03/2012

   private static final double NORM = (1.0 / 0x100000001L);

   private static final int W = 32;
   private static final int R = 16;
   private static final int P = 0;
   private static final int M1 = 13;
   private static final int M2 = 9;
   private static final int M3 = 5;
   private static final int MASK = 0xF;    // = 15

   //state variables
   private int state_i;
   private int[] state;

   //stream and substream variables :
   private int[] stream;
   private int[] substream;
   private static int[] curr_stream = new int[]
                              {0xA341BF9A, 0xAFE4901B, 0x6B10DE18, 0x05FE1420,
                               0xE48B1A9C, 0x590AE15E, 0xC5EB82A7, 0x37EAB2F9,
                               0x90E1C6EA, 0x3AE63902, 0x735DC91C, 0x902E3A8C,
                               0x6CB28A5D, 0x8474E7D1, 0x843E01A3, 0x5A7370EF};

   // P(z) = {0xa7600001, 0xe0f4f3e2, 0xcb30e185, 0x7d6b79a9,
   //         0xf3d46237, 0x13a524cb, 0x38e3c2d2, 0xa1381bcb,
   //         0xf7ab5f06, 0x04a72cda, 0x4e302521, 0xaca072f1,
   //         0x4dd96181, 0x24aa25c9, 0x3c417e7,  0x0, 0x1}

   // Ce tableau represente les 512 coefficients du polynome suivant
   // (z^(2^200) mod P(z)) mod 2
   // P(z) est le polynome caracteristique du generateur.
   private static final int [] pw = new int[]
                          {0x280009a9, 0x31e221d0, 0xa00c0296, 0x763d492b,
                           0x63875b75, 0xef2acc3a, 0x1400839f, 0x5e0c8526,
                            0x514e11b, 0x56b398e4, 0x9436c8b9, 0xa6d8130b,
                           0xc0a48a78, 0x26ad57d0, 0xa3a0c62a, 0x3ff16c9b};

   // Ce tableau represente les 512 coefficients du polynome suivant
   // (z^(2^350) mod P(z)) mod 2
   // P(z) est le polynome caracteristique du generateur.
   private static final int [] pz = new int[]
                          {0xcd68f2fe, 0x183e969a, 0x760449ae, 0xaa0ce54e,
                           0xfb5363af, 0x79deea9b, 0xef66c516, 0x103543cb,
                           0x244d1a97, 0x7570bc91, 0x31203fc7, 0x455ea2ca,
                           0xd77d327d, 0xd8c6a83c, 0xc51b05e7, 0x300c1501};


   private void advanceSeed(int[] seed, int [] p) {
      int b;
      int [] x = new int[R];

      for (int i = 0; i < R; i++) {
         state[i] = seed[i];
      }
      state_i = 0;

      for (int j = 0; j < R; ++j) {
         b = p[j];
         for (int k = 0; k < W; ++k) {
            if ((b & 1) == 1) {
               for (int i = 0; i < R; i++) {
                  x[i] ^= state[(state_i + i) & MASK];
               }
            }
            b >>= 1;

            nextValue();
         }
      }

      for (int i = 0; i < R; i++) {
         seed[i] = x[i];
      }
   }


   private static void verifySeed(int[] seed) {
      if (seed.length < R)
         throw new IllegalArgumentException("Seed must contain " + R +
                                            "values");
      for (int i = 0; i < R; i++)
         if (seed[i] != 0)
            return;
      throw new IllegalArgumentException
      ("At least one of the element of the seed must not be 0.");
   }

   private WELL512(int i) {
      //unit vector (to build the state transition matrice)
      state = new int[R];
      for(int j = 0; j < R; j++)
         state[j] = 0;
      state[i / W] = 1 << (i % W);
      state_i = 0;
   }

   /**
    * Constructs a new stream.
    */
   public WELL512() {
      state = new int[R];
      stream = new int[R];
      substream = new int[R];

      for(int i = 0; i < R; i++)
         stream[i] = curr_stream[i];

      advanceSeed(curr_stream, pz);
      resetStartStream();
   }

   /**
    * Constructs a new stream with the identifier `name` (used in the
    * `toString` method).
    *  @param name         name of the stream
    */
   public WELL512 (String name) {
      this();
      this.name = name;
   }

   /**
    * Sets the initial seed of the class `WELL512` to the 16 integers of
    * the vector `seed[0..15]`. This will be the initial seed of the class
    * of the next created stream. At least one of the integers must be
    * non-zero.
    *  @param seed         array of 16 elements representing the seed
    */
   public static void setPackageSeed (int seed[]) {
      verifySeed(seed);
      for(int i = 0; i < R; i++)
         curr_stream[i] = seed[i];
   }

   /**
    * This method is discouraged for normal use. Initializes the stream at
    * the beginning of a stream with the initial seed `seed[0..15]`. The
    * seed must satisfy the same conditions as in `setPackageSeed`. This
    * method only affects the specified stream; the others are not
    * modified. Hence after calling this method, the beginning of the
    * streams will no longer be spaced @f$Z@f$ values apart. For this
    * reason, this method should only be used in very exceptional cases;
    * proper use of the `reset...` methods and of the stream constructor
    * is preferable.
    *  @param seed         array of 16 elements representing the seed
    */
   public void setSeed (int seed[]) {
      verifySeed(seed);
      for(int i = 0; i < R; i++)
         stream[i] = seed[i];
      resetStartStream();
   }

   /**
    * Returns the current state of the stream, represented as an array of
    * 16 integers.
    *  @return the current state of the stream
    */
   public int[] getState() {
      int[] result = new int[R];
      for(int i = 0; i < R; i++)
         result[i] = state[(state_i + i) & MASK];
      return result;
   }

   /**
    * Clones the current generator and return its copy.
    *  @return A deep copy of the current generator
    */
   public WELL512 clone() {
      WELL512 retour = null;

      retour = (WELL512)super.clone();
      retour.state = new int[R];
      retour.substream = new int[R];
      retour.stream = new int[R];

      for (int i = 0; i<R; i++) {
         retour.state[i] = state[i];
         retour.substream[i] = substream[i];
         retour.stream[i] = stream[i];
      }
      return retour;
   }

   public void resetStartStream() {
      for(int i = 0; i < R; i++)
         substream[i] = stream[i];
      resetStartSubstream();
   }

   public void resetStartSubstream() {
      state_i = 0;
      for(int i = 0; i < R; i++)
         state[i] = substream[i];
   }

   public void resetNextSubstream() {
      advanceSeed(substream, pw);
      resetStartSubstream();
   }

   public String toString()  {
      StringBuffer sb = new StringBuffer();

      if(name == null)
         sb.append("The state of this WELL512 is : {");
      else
         sb.append("The state of " + name + " is : {");
      for(int i = 0; i < R - 1; i++)
         sb.append(state[(state_i + i) & MASK] + ", ");
      sb.append(state[(state_i + R - 1) & MASK] + "}");

      return sb.toString();
   }

   protected double nextValue() {
      int z0, z1, z2;
      z0 = state[(state_i + 15) & MASK];
      z1 = (state[state_i] ^ (state[state_i] << 16)) ^
           (state[(state_i+M1) & MASK] ^ (state[(state_i+M1) & MASK] << 15));
      z2 = (state[(state_i+M2) & MASK] ^
           (state[(state_i+M2) & MASK] >>> 11));
      state[state_i] = z1 ^ z2;
      state[(state_i + 15) & MASK] = (z0 ^ (z0 << 2)) ^ (z1 ^ (z1 << 18)) ^
                       (z2 << 28) ^ (state[state_i] ^
                           ((state[state_i] << 5) & 0xDA442D24));
      state_i = (state_i + 15) & MASK;

      long result = state[state_i];

      return (double)(result > 0 ? result : (result + 0x100000000L)) * NORM;
   }

}