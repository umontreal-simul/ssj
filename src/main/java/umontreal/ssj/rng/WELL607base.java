package umontreal.ssj.rng;

import umontreal.ssj.util.BitVector;
import umontreal.ssj.util.BitMatrix;
import java.io.Serializable;
import java.io.*;


/*
 * Class:        WELL607base
 * Description:
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


abstract class WELL607base extends RandomStreamBase
{

   private static final long serialVersionUID = 120307L;
   // La date de modification a l'envers, lire 07/03/2012

   // constants
   static final double NORM = (1.0 / 0x100000001L);
   // static final double NORM = 2.32830643653869628906e-10;

   static final int R = 19;       // useful length of the state
   static final int BUFFER_SIZE = 32; // length of the state
   static final int NUM_BITS = 608;
   static final int MASK_STATE = 0x0000001F; // = 31
   static final int W = 32;
   static final int P = 1;
   static final int MASKU = (0xffffffff >>> (W - P));
   static final int MASKL = (~MASKU);
   static final int M1 = 16;
   static final int M2 = 15;
   static final int M3 = 14;
   static final int R1 = 18;      // R - 1
   static final int R2 = 17;      // R - 2

   // state variables
   int state_i;
   int[] state;

   // stream and substream
   int[] stream;
   int[] substream;

   // length of the jumps
   static final int w = 250;
   static final int v = 150;

   // advance the state by a transition matrice
   protected void advanceSeed(int[] seed, int [] p)
   {
      int b;
      int[] x = new int[R];

      for (int i = 0; i < R; i++) {
         state[i] = seed[i];
      }
      state_i = 0;

      for (int j = 0; j < R; ++j) {
         b = p[j];
         for (int k = 0; k < W; ++k) {
            if ((b & 1) == 1) {
               for (int i = 0; i < R; i++) {
                  x[i] ^= state[(state_i + i) & MASK_STATE];
               }
            }
            b >>= 1;
            nextInt ();
         }
      }

      for (int i = 0; i < R; i++) {
         seed[i] = x[i];
      }
   }


   static void verifySeed (int seed[])
   {
      if (seed.length < R)
         throw new IllegalArgumentException ("Seed must contain " + R +
                                             " values");

      boolean goodSeed = false;
      for (int i = 0; !goodSeed && i < R; i++)
         if (seed[i] != 0)
            goodSeed = true;
      if (!goodSeed)
         if (seed[R - 1] == 0x80000000)
            throw new IllegalArgumentException
            ("At least one of the element of the seed must not be 0. " +
             "If this element is the last one, it mustn't be equal " +
             "to 0x80000000 (" + 0x80000000 + ").");
   }


   int[] getState ()
   {
      int[] result = new int[R];
      for (int i = 0; i < R; i++)
         result[i] = state[(state_i + i) & MASK_STATE];
      return result;
   }


   // just like formatState, but not public
   String stringState ()
   {
      StringBuffer sb = new StringBuffer ();
      for (int i = 0; i < R - 1; i++)
         sb.append (state[(state_i + i) & MASK_STATE] + ", ");
      sb.append (state[(state_i + R - 1) & MASK_STATE] + "}");
      return sb.toString ();
   }


   int nextInt ()
   {
      int z0, z1, z2;

      z0 = (state[(state_i + R1) & MASK_STATE] & MASKL) |
           (state[(state_i + R2) & MASK_STATE] & MASKU);
      z1 = (state[state_i] ^ (state[state_i] >>> 19)) ^
           (state[(state_i + M1) & MASK_STATE] ^
            (state[(state_i + M1) & MASK_STATE] >>> 11));
      z2 = (state[(state_i + M2) & MASK_STATE] ^
            (state[(state_i + M2) & MASK_STATE] << (14))) ^
           state[(state_i + M3) & MASK_STATE];
      state[state_i] = z1 ^ z2;
      state[(state_i - 1) & MASK_STATE] = (z0 ^ (z0 >>> 18)) ^
                                          z1 ^ (state[state_i] ^ (state[state_i] << 5));

      state_i--;
      state_i &= MASK_STATE;
      return state[state_i];
   }


   public WELL607base clone ()
   {
      WELL607base retour = null;
      retour = (WELL607base) super.clone ();
      return retour;
   }
}
