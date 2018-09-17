/*
 * Class:        WELL607
 * Description:  a Well Equidistributed Long period Linear Random Number
                 Generator with a state size of 607 bits
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
 * proposed by F. Panneton in @cite rPAN06b, @cite rPAN04t&thinsp;. The
 * implemented generator is the `WELL607`, which has a state size of 607 bits
 * and a period length of  @f$\rho\approx@f$ @f$2^{607}@f$. The values of
 * @f$V@f$, @f$W@f$ and @f$Z@f$ are @f$2^{150}@f$, @f$2^{250}@f$ and
 * @f$2^{400}@f$ respectively (see  @ref RandomStream for their definition).
 * The seed of the RNG, and the state of a stream at any given step, is a
 * 19-dimensional vector of 32-bit integers. The output of `nextValue` has 32
 * bits of precision.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class WELL607 extends WELL607base {

   private static final long serialVersionUID = 70510L;
   //La date de modification a l'envers, lire 10/05/2007

   private static int[] curr_stream = {0xD6AFB71C, 0x82ADB18E, 0x326E714E,
                                       0xB1EE42B6, 0xF1A834ED, 0x04AE5721,
                                       0xC5EA2843, 0xFA04116B, 0x6ACE14EF,
                                       0xCD5781A0, 0x6B1F731C, 0x7E3B8E3D,
                                       0x8B34DE2A, 0x74EC15F5, 0x84EBC216,
                                       0x83EA2C61, 0xE4A83B1E, 0xA5D82CB9,
                                       0x9E1A6C89};

   // P(z) = {0x987b2631, 0x2e33283d, 0x6a398474, 0xe9d24da1,
   //         0x31235359, 0x6a2baf48, 0x7f97efd4, 0x468280f4,
   //         0x7d9d9424, 0xa3238f8e, 0xe3edb4ef, 0x0e0a25f7,
   //         0x92c4dff5, 0x55d0b8da, 0x7b982dec, 0xa06c078f,
   //         0x38b65c31, 0xc8c3788d, 0x8000b200}

   // Ce tableau represente les 512 coefficients du polynome suivant
   // (z^(2^250) mod P(z)) mod 2
   // P(z) est le polynome caracteristique du generateur.
   static final int [] pw = new int[]
                     {0x83167621, 0x6b5515c8, 0x61a62bd2, 0xbceaa78f,
                      0xac04b304, 0x28a75ea4, 0xa9104058, 0x595ea53b,
                      0x35687e95, 0x7f8eca9b, 0x30beffb8, 0xc61e6111,
                      0x284ee30e, 0x4e9cd901, 0x659633ba, 0x344cc69e,
                      0xd6052ac1, 0x5d508b69,  0x62cf130};

   // Ce tableau represente les 512 coefficients du polynome suivant
   // (z^(2^400) mod P(z)) mod 2
   // P(z) est le polynome caracteristique du generateur.
   static final int [] pz = new int[]
                     {0x70b2bdee, 0x595828f1, 0x85a17885, 0x5100c7b2,
                      0xd3333da2, 0xb42857de, 0xf8a7a4a7, 0xabad2a33,
                       0xa2580cf, 0xf94c465e, 0x7df951d5, 0x35467053,
                       0xb3c9a4e,  0x6a33977,  0x443910e, 0xc25aec3d,
                      0xeb72e8c5,  0x8873b01, 0x7da57636};

   /**
    * Constructs a new stream.
    */
   public WELL607() {
      state = new int[BUFFER_SIZE];
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
   public WELL607 (String name) {
      this();
      this.name = name;
   }

   /**
    * Sets the initial seed of the class `WELL607` to the 19 integers of
    * the vector `seed[0..18]`. This will be the initial seed of the next
    * created stream. At least one of the integers must not be zero and if
    * this integer is the last one, it must not be equal to `0x80000000`.
    *  @param seed         array of 19 elements representing the seed
    */
   public static void setPackageSeed (int seed[]) {
      verifySeed(seed);
      for(int i = 0; i < R; i++)
         curr_stream[i] = seed[i];
   }

   /**
    * This method is discouraged for normal use. Initializes the stream at
    * the beginning of a stream with the initial seed `seed[0..18]`. The
    * seed must satisfy the same conditions as in `setPackageSeed`. This
    * method only affects the specified stream; the others are not
    * modified. Hence after calling this method, the beginning of the
    * streams will no longer be spaced @f$Z@f$ values apart. For this
    * reason, this method should only be used in very exceptional cases;
    * proper use of the `reset...` methods and of the stream constructor
    * is preferable.
    *  @param seed         array of 19 elements representing the seed
    */
   public void setSeed (int seed[]) {
      verifySeed(seed);
      for(int i = 0; i < R; i++)
         stream[i] = seed[i];
      resetStartStream();
   }

   /**
    * Returns the current state of the stream, represented as an array of
    * 19 integers.
    *  @return the current state of the stream
    */
   public int[] getState() {
      return super.getState();
   }


   public void resetStartStream() {
      for(int i = 0; i < R; i++)
         substream[i] = stream[i];
      resetStartSubstream();
   }

   public void resetStartSubstream() {
      state_i= 0;
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
         sb.append("The state of this WELL607 is : {");
      else
         sb.append("The state of " + name + " is : {");
      sb.append(super.stringState());
      return sb.toString();
   }

   protected double nextValue() {
      long result = nextInt();
      if(result <= 0)
         result += 0x100000000L;
      return result * NORM;
   }

/**
 * Clones the current generator and return its copy.
 *  @return A deep copy of the current generator
 */
public WELL607 clone() {
      WELL607 retour = null;
      retour = (WELL607)super.clone();
      retour.state = new int[BUFFER_SIZE];
      retour.substream = new int[R];
      retour.stream = new int[R];

      for (int i = 0; i<R; i++) {
         retour.substream[i] = substream[i];
         retour.stream[i] = stream[i];
      }
      for (int i = 0; i<BUFFER_SIZE; i++) {
         retour.state[i] = state[i];
      }

      return retour;
   }

}