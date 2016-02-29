/*
 * Class:        MRG31k3p
 * Description:  combined multiple recursive generator proposed by L'Ecuyer and Touzin
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


import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.rng.RandomStreamBase;
import umontreal.ssj.util.ArithmeticMod;
import java.io.Serializable;

/**
 * Extends the abstract class  @ref RandomStreamBase, thus implementing the
 * @ref RandomStream interface indirectly. The backbone generator is the
 * combined multiple recursive generator (CMRG) `MRG31k3p` proposed by
 * L’Ecuyer and Touzin @cite rLEC00b&thinsp;, implemented in 32-bit integer
 * arithmetic. This RNG has a period length of  @f$\rho\approx@f$
 * @f$2^{185}@f$. The values of @f$V@f$, @f$W@f$ and @f$Z@f$ are
 * @f$2^{62}@f$, @f$2^{72}@f$ and @f$2^{134}@f$ respectively. (See
 * @ref RandomStream for their definition.) The seed and the state of a
 * stream at any given step are six-dimensional vectors of 32-bit integers.
 * The default initial seed is @f$(12345, 12345, 12345, 12345, 12345,
 * 12345)@f$. The method `nextValue` provides 31 bits of precision.
 *
 * The difference between the RNG of class  @ref MRG32k3a and this one is
 * that this one has all its coefficients of the form @f$a = \pm2^q
 * \pm2^r@f$. This permits a faster implementation than for arbitrary
 * coefficients. On a 32-bit computer, `MRG31k3p` is about twice as fast as
 * `MRG32k3a`. On the other hand, the latter does a little better in the
 * spectral test and has been more extensively tested.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class MRG31k3p extends RandomStreamBase {

   private static final long serialVersionUID = 70510L;
   //La date de modification a l'envers, lire 10/05/2007

   //generator constants :
   private static final int M1 = 2147483647;    //2^31 - 1
   private static final int M2 = 2147462579;    //2^31 - 21069
   private static final int MASK12 = 511;       //2^9 - 1
   private static final int MASK13 = 16777215;  //2^24 - 1
   private static final int MASK2 = 65535;      //2^16 - 1
   private static final int MULT2 = 21069;
   private static final double NORM = 4.656612873077392578125e-10;
   //private static final double NORM = 1.0 / (M1 + 1.0);

   //state variables
   private int x11;
   private int x12;
   private int x13;
   private int x21;
   private int x22;
   private int x23;

   //stream and substream variables
   private int[] stream;
   private int[] substream;
   private static int[] curr_stream = {12345, 12345, 12345,
                                       12345, 12345, 12345};

   //streams constants :
   private static final int[][] A1p0 =
      {{0, 4194304, 129},
       {1, 0, 0},
       {0, 1, 0}};
   private static final int[][] A2p0 =
      {{32768, 0, 32769},
       {1, 0, 0},
       {0, 1, 0}};

   private static final int[][] A1p72 =
      {{1516919229, 758510237, 499121365},
       {1884998244, 1516919229, 335398200},
       {601897748, 1884998244, 358115744}};
   private static final int[][] A2p72 =
      {{1228857673, 1496414766, 954677935},
       {1133297478, 1407477216, 1496414766},
       {2002613992, 1639496704, 1407477216}};

   private static final int[][] A1p134 =
      {{1702500920, 1849582496, 1656874625},
       {828554832, 1702500920, 1512419905},
       {1143731069, 828554832, 102237247}};
   private static final int[][] A2p134 =
      {{796789021, 1464208080, 607337906},
       {1241679051, 1431130166, 1464208080},
       {1401213391, 1178684362, 1431130166}};


   //multiply the first half of v by A with a modulo of m1
   //and the second half by B with a modulo of m2
   private static void multMatVect(int[] v, int[][] A, int m1,
                                   int[][] B, int m2) {
      int[] vv = new int[3];
      for(int i = 0; i < 3; i++)
         vv[i] = v[i];
      ArithmeticMod.matVecModM(A, vv, vv, m1);
      for(int i = 0; i < 3; i++)
         v[i] = vv[i];

      for(int i = 0; i < 3; i++)
         vv[i] = v[i + 3];
      ArithmeticMod.matVecModM(B, vv, vv, m2);
      for(int i = 0; i < 3; i++)
         v[i + 3] = vv[i];

   }

   /**
    * Constructs a new stream, initialized at its beginning. Its seed is
    * @f$Z = 2^{134}@f$ steps away from the previous seed.
    */
   public MRG31k3p() {
      name = null;

      prec53 = false;
      anti = false;

      stream = new int[6];
      substream = new int[6];
      for(int i = 0; i < 6; i++)
         stream[i] = curr_stream[i];

      resetStartStream();

      multMatVect(curr_stream, A1p134, M1, A2p134, M2);
   }

   /**
    * Constructs a new stream with the identifier `name` (used when
    * formatting the stream state).
    *  @param name         name of the stream
    */
   public MRG31k3p (String name) {
      this();
      this.name = name;
   }

   /**
    * Sets the initial seed for the class `MRG31k3p` to the six integers
    * of the vector `seed[0..5]`. This will be the initial state (or seed)
    * of the next created stream. By default, if this method is not
    * called, the first stream is created with the seed @f$(12345, 12345,
    * 12345, 12345, 12345, 12345)@f$. If it is called, the first 3 values
    * of the seed must all be less than @f$m_1 = 2147483647@f$, and not
    * all 0; and the last 3 values must all be less than @f$m_2 =
    * 2147462579@f$, and not all 0.
    *  @param seed         array of 6 elements representing the seed
    */
   public static void setPackageSeed (int seed[]) {
      if (seed.length < 6)
         throw new IllegalArgumentException ("Seed must contain 6 values");
      if (seed[0] == 0 && seed[1] == 0 && seed[2] == 0)
         throw new IllegalArgumentException ("The first 3 values must not be 0");
      if (seed[5] == 0 && seed[3] == 0 && seed[4] == 0)
         throw new IllegalArgumentException ("The last 3 values must not be 0");
      if (seed[0] >= M1 || seed[1] >= M1 || seed[2] >= M1)
         throw new IllegalArgumentException ("The first 3 values must be less than " + M1);

      if (seed[5] >= M2 || seed[3] >= M2 || seed[4] >= M2)
         throw new IllegalArgumentException ("The last 3 values must be less than " + M2);
      for (int i = 0; i < 6;  ++i)
         curr_stream[i] = seed[i];
   }

   /**
    * *Use of this method is strongly discouraged*. Initializes the stream
    * at the beginning of a stream with the initial seed `seed[0..5]`.
    * This vector must satisfy the same conditions as in `setPackageSeed`.
    * This method only affects the specified stream, all the others are
    * not modified, so the beginning of the streams are no longer spaced
    * @f$Z@f$ values apart. For this reason, this method should be used
    * only in very exceptional situations; proper use of `reset...` and of
    * the stream constructor is preferable.
    *  @param seed         array of 6 integers representing the new seed
    */
   public void setSeed (int seed[]) {
      if (seed.length < 6)
         throw new IllegalArgumentException ("Seed must contain 6 values");
      if (seed[0] == 0 && seed[1] == 0 && seed[2] == 0)
         throw new IllegalArgumentException ("The first 3 values must not be 0");
      if (seed[3] == 0 && seed[4] == 0 && seed[5] == 0)
         throw new IllegalArgumentException ("The last 3 values must not be 0");
      if (seed[0] >= M1 || seed[1] >= M1 || seed[2] >= M1)
         throw new IllegalArgumentException ("The first 3 values must be less than " + M1);
      if (seed[3] >= M2 || seed[4] >= M2 || seed[5] >= M2)
         throw new IllegalArgumentException ("The last 3 values must be less than " + M2);
      for (int i = 0; i < 6;  ++i)
         stream[i] = seed[i];
      resetStartStream();
   }


   public void resetStartStream() {
      for(int i = 0; i < 6; i++)
         substream[i] = stream[i];
      resetStartSubstream();
   }

   public void resetStartSubstream() {
      x11 = substream[0];
      x12 = substream[1];
      x13 = substream[2];
      x21 = substream[3];
      x22 = substream[4];
      x23 = substream[5];
   }

   public void resetNextSubstream() {
      multMatVect(substream, A1p72, M1, A2p72, M2);
      resetStartSubstream();
   }

/**
 * Returns the current state @f$C_g@f$ of this stream. This is a vector of 6
 * integers represented. This method is convenient if we want to save the
 * state for subsequent use.
 *  @return the current state of the generator
 */
public int[] getState() {
      return new int[]{x11, x12, x13, x21, x22, x23};
   }

   /**
    * Clones the current generator and return its copy.
    *  @return A deep copy of the current generator
    */
   public MRG31k3p clone() {
      MRG31k3p retour = null;
      retour = (MRG31k3p)super.clone();
      retour.substream = new int[6];
      retour.stream = new int[6];
      for (int i = 0; i<6; i++) {
         retour.substream[i] = substream[i];
         retour.stream[i] = stream[i];
      }
      return retour;
   }

   public String toString() {
      if(name == null)
         return "The state of the MRG31k3p is: " +
                x11 + ", " + x12 + ", " + x13 + ";  " +
                x21 + ", " + x22 + ", " + x23;
      else
         return "The state of " + name + " is: " +
                x11 + ", " + x12 + ", " + x13 + ";  " +
                x21 + ", " + x22 + ", " + x23;
   }

   protected double nextValue()  {
      int y1, y2;

      //first component
      y1 = ((x12 & MASK12) << 22) + (x12 >>> 9)
           + ((x13 & MASK13) << 7) + (x13 >>> 24);
      if(y1 < 0 || y1 >= M1)     //must also check overflow
         y1 -= M1;
      y1 += x13;
      if(y1 < 0 || y1 >= M1)
         y1 -= M1;

      x13 = x12;
      x12 = x11;
      x11 = y1;

      //second component
      y1 = ((x21 & MASK2) << 15) + (MULT2 * (x21 >>> 16));
      if(y1 < 0 || y1 >= M2)
         y1 -= M2;
      y2 = ((x23 & MASK2) << 15) + (MULT2 * (x23 >>> 16));
      if(y2 < 0 || y2 >= M2)
         y2 -= M2;
      y2 += x23;
      if(y2 < 0 || y2 >= M2)
         y2 -= M2;
      y2 += y1;
      if(y2 < 0 || y2 >= M2)
         y2 -= M2;

      x23 = x22;
      x22 = x21;
      x21 = y2;

      //Must never return either 0 or 1
      if(x11 <= x21)
         return (x11 - x21 + M1) * NORM;
      else
         return (x11 - x21) * NORM;
   }

}