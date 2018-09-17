/*
 * Class:        WELL1024
 * Description:  a Well Equidistributed Long period Linear Random Number
                 Generator with a state size of 1024 bits
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
 * Implements the  @ref RandomStream interface via inheritance from
 * @ref RandomStreamBase. The backbone generator is a Well Equidistributed
 * Long period Linear Random Number Generator (WELL), proposed by F. Panneton
 * in @cite rPAN06b, @cite rPAN04t&thinsp;, and which has a state size of
 * 1024 bits and a period length of  @f$\rho\approx@f$ @f$2^{1024}@f$. The
 * values of @f$V@f$, @f$W@f$ and @f$Z@f$ are @f$2^{300}@f$, @f$2^{400}@f$
 * and @f$2^{700}@f$ respectively (see  @ref RandomStream for their
 * definition). The seed of the RNG, and the state of a stream at any given
 * step, is a 16-dimensional vector of 32-bit integers. The output of
 * `nextValue` has 32 bits of precision.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class WELL1024 extends RandomStreamBase {

   private static final long serialVersionUID = 120307L;
   //La date de modification a l'envers, lire 12/03/07

   private static final double NORM = 2.32830643653869628906e-10;

   private static final int MASK = 0x1F;

   private static final int W = 32;
   private static final int R = 32;
   private static final int P = 0;
   private static final int M1 = 3;
   private static final int M2 = 24;
   private static final int M3 = 10;

   private static final int A1 = 0xDB79FB31;
   private static final int B1 = 0x0FDAA2C1;
   private static final int B2 = 0x27C5A58D;
   private static final int C1 = 0x71E993FE;

   private int state_i = 0;
   private int[] state;

   //stream and substream variables
   private int[] stream;
   private int[] substream;
   private static int[] curr_stream =
                    new int[] {0xDE410B75, 0x904FA5C7, 0x8BD4701E, 0x011EA361,
                               0x6EB189E0, 0x7A2B0CE1, 0xE02631CA, 0x72EBA132,
                               0x5189DA0F, 0x3EB72A2C, 0x51ABE513, 0x6D9EA57C,
                               0x4D690BF1, 0x84217FCA, 0x7290DE1A, 0x429F5A48,
                               0x6EC42EF3, 0x960AB315, 0x72C3A743, 0x48E13BF1,
                               0x8917EAC8, 0x284AE026, 0x357BF240, 0x913B51AC,
                               0x136AF195, 0x361ABC18, 0x731AB725, 0x63D3A7C9,
                               0xE5F32A18, 0x91A8E164, 0x04EA61B5, 0xC72A6091};


   // P(z) = {0x00000001, 0x00000000, 0x00000000, 0x00000000,
   //         0x028a0008, 0x02288020, 0x2baaa20a, 0x0209aa00,
   //         0x3f871248, 0x80172a7b, 0xee101d14, 0xef2221f3,
   //         0xb5bf7be1, 0xab57e80c, 0xfa24ee53, 0x37dab9aa,
   //         0xd353180b, 0xf1c5d9ed, 0xd6465866, 0x7a048625,
   //         0x892b7ef6, 0x2ca9170f, 0xa8a3f324, 0x36be065f,
   //         0x57aee2ab, 0xb20f4dd9, 0xa0eaa2ee, 0xa678c37a,
   //         0x5792d2ae, 0xac449456, 0x51549f89, 0x0, 0x1}

   // Ce tableau represente les 1024 coefficients du polynome suivant
   // (z^(2^400) mod P(z)) mod 2
   // P(z) est le polynome caracteristique du generateur.
   private static final int [] pw = new int[]
                        { 0xe44294e, 0xef237eff, 0x5e8b6bfb, 0xa724e67a,
                         0x59994cfd, 0x6f7c3de1, 0x6735d50d, 0x4bfe199a,
                         0x39c28e61, 0xfd075266, 0x96cc6d1f, 0x5dc1a685,
                         0xd67fa444, 0xccc01b86,  0x8ff861c, 0xce113725,
                         0x66707603, 0x38abb0fd,  0x7681f64, 0x104535c5,
                         0xce4ae5f4, 0x50e37105, 0xd0c5f77f, 0x74c1ebf6,
                         0x2ccf1505, 0xd1f21b86, 0x9a6c402e, 0xea34a31c,
                         0x65e13d13, 0xde8f2f05, 0x89db804f, 0x8dc387f2};

   // Ce tableau represente les 1024 coefficients du polynome suivant
   // (z^(2^700) mod P(z)) mod 2
   // P(z) est le polynome caracteristique du generateur.
   private static final int [] pz = new int[]
                        {0x7cab7da4, 0xef28b275, 0x18ffa66a, 0x2aa41e52,
                         0x15b6bd86, 0x560d0d76, 0xcdeda011, 0x96231727,
                         0xeec6a7f2, 0x99fd2be6, 0x92afa886, 0xcca777f0,
                         0x972eff38, 0xa29f8e49, 0x22b4b9b6, 0x1089c898,
                         0x6d569b25, 0x879044c2, 0x5e41b523, 0x33f19dd6,
                         0x7c005fc5, 0x7f9a1907, 0x39bf9eed, 0x4bd86a74,
                          0xe1e47e3, 0x96ead7ac, 0xc834f9ee, 0xd9ff4a4f,
                         0x717f044c, 0xfd0e15e6,  0x6c18ef3, 0xbfdd2942};



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
                                            " values");
      for(int i = 0; i < R; i++)
         if (seed[i] != 0)
            return;
      throw new IllegalArgumentException
      ("At least one of the element of the seed must not be 0.");
   }

   private WELL1024(int i) {
      //unit vector (to build the state transition matrices)
      state = new int[R];
      for(int j = 0; j < R; j++)
         state[j] = 0;
      state[i / W] = 1 << (i % W);
      state_i = 0;
   }

   /**
    * Constructs a new stream.
    */
   public WELL1024() {
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
   public WELL1024 (String name) {
      this();
      this.name = name;
   }

   /**
    * Sets the initial seed of this class to the 32 integers of array
    * `seed[0..31]`. This will be the initial seed of the class and of the
    * next created stream. At least one of the integers must be non-zero.
    *  @param seed         array of 32 elements representing the seed
    */
   public static void setPackageSeed (int seed[]) {
      verifySeed (seed);
      for(int i = 0 ; i < R; i++)
         curr_stream[i] = seed[i];
   }

   /**
    * This method is discouraged for normal use. Initializes the stream at
    * the beginning of a stream with the initial seed `seed[0..31]`. The
    * seed must satisfy the same conditions as in `setPackageSeed`. This
    * method only affects the specified stream; the others are not
    * modified. Hence after calling this method, the beginning of the
    * streams will no longer be spaced @f$Z@f$ values apart. For this
    * reason, this method should only be used in very exceptional cases;
    * proper use of the `reset...` methods and of the stream constructor
    * is preferable.
    *  @param seed         array of 32 elements representing the seed
    */
   public void setSeed (int seed[]) {
      verifySeed (seed);
      for(int i = 0 ; i <  R; i ++)
         stream[i] = seed[i];
      resetStartStream();
   }

   /**
    * Returns the current state of the stream, represented as an array of
    * 32 integers.
    *  @return the current state of the stream
    */
   public int[] getState() {
      int[] result = new int[R];
      for(int i = 0 ; i < R; i ++)
         result[i] = state[(state_i + i) & MASK];
      return result;
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
         sb.append("The state of this WELL1024 is : {");
      else
         sb.append("The state of " + name + " is : {");
      for(int i = 0; i < R - 1; i++)
         sb.append(state[(state_i + i) & MASK] + ", ");
      sb.append(state[(state_i + R - 1) & MASK] + "}");

      return sb.toString();
   }

   protected double nextValue() {
      int z0, z1, z2;

      z0    = state[(state_i + 31) & MASK];
      z1    = state[state_i] ^ (state[(state_i + M1) & MASK] ^
                                (state[(state_i + M1) & MASK] >>> 8));
      z2    = (state[(state_i + M2) & MASK] ^
               (state[(state_i + M2) & MASK] << 19)) ^
              (state[(state_i + M3) & MASK] ^
               (state[(state_i + M3) & MASK] << 14));
      state[state_i] = z1 ^ z2;
      state[(state_i + 31) & MASK] = (z0 ^ (z0 << 11)) ^
                                     (z1 ^ (z1 << 7)) ^ (z2 ^ (z2 << 13));
      state_i = (state_i + 31) & MASK;

      long result = state[state_i];

      return ((double) (result > 0 ? result : result + 0x100000000L) * NORM);
   }

/**
 * Clones the current generator and return its copy.
 *  @return A deep copy of the current generator
 */
public WELL1024 clone() {
      WELL1024 retour = null;
      retour = (WELL1024)super.clone();
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

}