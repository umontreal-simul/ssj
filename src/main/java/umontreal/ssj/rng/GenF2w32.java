/*
 * Class:        GenF2w32
 * Description:  generator proposed by Panneton
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

import umontreal.ssj.util.BitVector;
import umontreal.ssj.util.BitMatrix;
import java.io.*;

/**
 * Implements the  @ref RandomStream interface via inheritance from
 * @ref RandomStreamBase. The backbone generator is a Linear Congruential
 * Generator (LCG) in the finite field @f$\mathbb F_{2^w}@f$ instead of
 * @f$\mathbb F_2@f$. The implemented generator is the <tt>GenF2w2_32</tt>
 * proposed by Panneton @cite rPAN04a, @cite rPAN04t&thinsp;. Its state is 25
 * 32-bit words and it has a period length of @f$2^{800} - 1@f$. The values
 * of @f$V@f$, @f$W@f$ and @f$Z@f$ are @f$2^{200}@f$, @f$2^{300}@f$ and
 * @f$2^{500}@f$ respectively (see  @ref RandomStream for their definition).
 * The seed of the RNG, and the state of a stream at any given step, is a
 * 25-dimensional vector of 32-bits integers. Its `nextValue` method returns
 * numbers with 32 bits of precision.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class GenF2w32 extends RandomStreamBase {

   private static final long serialVersionUID = 120307L;
   //La date de modification a l'envers, lire 07/03/2012

   //tas de constantes
   private final static double NORM = 1.0 / (0x100000001L);

   private final static int Z0 = 0x80000000;

   private final static int W = 32;

   private final static int R = 25;
   private final static int T = 7;

   private final static int W1 = 11, W2 = 11, W3 = 10;
   private final static int Wsplit1 = W - W1;
   private final static int Wsplit2 = Wsplit1 - W2;
   private final static int MASK1 = 0xFFFFFFFF << Wsplit1;
   private final static int MASK2 = (0xFFFFFFFF << Wsplit2) ^ MASK1;
   private final static int MASK3 = 0xFFFFFFFF ^ MASK1 ^ MASK2;

   private static int[] BrmT1, BrmT2, BrmT3;
   private static int[] Br1, Br2, Br3;


   private final static int BrmT = 0xE6A68D20;
   private final static int Br = 0x287AB842;

   private final static int modQ = 0xFA4F9B3F;

   //stream variables
   private final static int w = 300;
   private final static int v = 200;
   //private static F2wPoly.F2wPolyElem jumpW;
   //private static F2wPoly.F2wPolyElem jumpZ;
   private static BitMatrix Apw;
   private static BitMatrix Apz;

   //private static F2wPoly polyBase;

   /*
   private static F2wPoly.F2wPolyElem curr_stream;

   private F2wPoly.F2wPolyElem stream;
   private F2wPoly.F2wPolyElem substream;
   */
   private int[] stream;
   private int[] substream;
   private static int[] curr_stream;

   //if the generator was initialised
   private static boolean initialised = false;

   private int[] state;
   private int state_i;


   static
   {
      initTables();

      /*
      F2wPoly.F2w f2wBase = new F2wPoly.F2w(modQ);

      polyBase = new F2wPoly
                 (R, new F2wPoly.F2w.F2wElem[]{f2wBase.createElem(Br),
                                               f2wBase.createElem(BrmT)},
                  new int[]{0,T}, f2wBase);


      curr_stream = polyBase.createElem
                    (new int[]{0x95F24DAB, 0x0B685215, 0xE76CCAE7, 0xAF3EC239,
                               0x715FAD23, 0x24A590AD, 0x69E4B5EF, 0xBF456141,
                               0x96BC1B7B, 0xA7BDF825, 0xC1DE75B7, 0x8858A9C9,
                               0x2DA87693, 0xB657F9DD, 0xFFDC8A8F, 0x8121DA71,
                               0x8B823ECB, 0x885D05F5, 0x4E20CD47, 0x5A9AD5D9,
                               0x512C0C03, 0xEA857CCD, 0x4CC1D30F, 0x8891A8A1,
                               0xA6B7AADB});

      curr_stream = new int[]{0x95F24DAB, 0x0B685215, 0xE76CCAE7, 0xAF3EC239,
                              0x715FAD23, 0x24A590AD, 0x69E4B5EF, 0xBF456141,
                              0x96BC1B7B, 0xA7BDF825, 0xC1DE75B7, 0x8858A9C9,
                              0x2DA87693, 0xB657F9DD, 0xFFDC8A8F, 0x8121DA71,
                              0x8B823ECB, 0x885D05F5, 0x4E20CD47, 0x5A9AD5D9,
                              0x512C0C03, 0xEA857CCD, 0x4CC1D30F, 0x8891A8A1,
                              0xA6B7AADB};

      // Il n'est pas utile de fournir les resultats des exponentiation
      // puisqu'elles ne prennent qu'un temps negligeable par rapport
      // a la creation des tables de pre-calculs.


      jumpW = polyBase.createElem(new int[]
         {0, 1960343840, 1516339037, -333505122, -1976104464, -482584330,
          -1279385510, -1806644808, -1348555459, -588715441, 2019804456,
          533392567, 1830622053, -348374534, 335303887, 193005475, 2020690292,
          1810924850, -1106874017, 988574120, -581662300, 1744525859,
          697086516, 1542714170, -916619518});

      jumpZ = polyBase.createElem(new int[]
         {0, -946485540, 1238370284, -1390769524, -1573619753, -524893823,
          -873333569, 600407275, 631551291, 1484250168, 81475604, 721883080,
          767883194, 1119863295, -948128525, -2061810380, 631461516,
          351243833, 417985878, 828823193, 1965856614, -673630318,
          -2078590119, -1319538924, 838929999});


      jumpW = polyBase.createZ().exponentiateBase2(w);
      jumpZ = jumpW.exponentiateBase2(v);



      //transfert des vecteurs dans des matrices

      BitVector[] bv = new BitVector[800];
      int[] vect = new int[R];

      for(int i = 0; i < 800; i++)
         {
            for(int j = 0; j < vect.length; j++)
               vect[j] = 0;
            vect[i / 32] = 1 << (i % 32);
            polyBase.createElem(vect).multiply(jumpW).copyTo(vect);
            bv[i] = new BitVector(vect, 800);
         }
      Apw = (new BitMatrix(bv)).transpose();


      for(int i = 0; i < 800; i++)
         {
            for(int j = 0; j < vect.length; j++)
               vect[j] = 0;
            vect[i / 32] = 1 << (i % 32);
            polyBase.createElem(vect).multiply(jumpZ).copyTo(vect);
            bv[i] = new BitVector(vect, 800);
         }
      Apz = (new BitMatrix(bv)).transpose();
      */
   }



   static private void initialisation() {
      //initialise all of the state variables

      curr_stream = new int[]{0x95F24DAB, 0x0B685215, 0xE76CCAE7, 0xAF3EC239,
                              0x715FAD23, 0x24A590AD, 0x69E4B5EF, 0xBF456141,
                              0x96BC1B7B, 0xA7BDF825, 0xC1DE75B7, 0x8858A9C9,
                              0x2DA87693, 0xB657F9DD, 0xFFDC8A8F, 0x8121DA71,
                              0x8B823ECB, 0x885D05F5, 0x4E20CD47, 0x5A9AD5D9,
                              0x512C0C03, 0xEA857CCD, 0x4CC1D30F, 0x8891A8A1,
                              0xA6B7AADB};

      //reading the state transition matrices

      try {
         InputStream is = GenF2w32.class.getClassLoader().
                          getResourceAsStream ("umontreal/ssj/rng/GenF2w32.dat");
         ObjectInputStream ois = new ObjectInputStream(is);
         Apw = (BitMatrix)ois.readObject();
         Apz = (BitMatrix)ois.readObject();
         ois.close();

      } catch(FileNotFoundException e) {
         System.err.println("Couldn't find GenF2w32.dat");
         e.printStackTrace();
         throw new RuntimeException("  initialisation of GenF2w32");
      } catch(IOException e) {
         e.printStackTrace();
         throw new RuntimeException("  initialisation of GenF2w32");
      } catch(ClassNotFoundException e) {
         e.printStackTrace();
         throw new RuntimeException("  initialisation of GenF2w32");
      }

      initialised = true;
   }



   /*
     Calcule a*z^k dans GF(2^32)
   */
   private static int multiplyZ (int a, int k, int modPoly) {
      for(int i = 0; i < k; i++)
         if((a & 1) != 0)
            a = (a >>> 1) ^ modPoly;
         else
            a >>>= 1;

      return a;
   }

   /*
     Calcule a * b dans GF(2^32)
   */
   private static int multiply (int a, int b, int modPoly) {
      int res = 0;
      int verif = 1;
      for(int i = 0; i < W; i++) {
         if((b & verif) != 0)
            res ^= multiplyZ(a, W - 1 - i, modPoly);
         verif <<= 1;
      }

      return res;
   }

   /*
     Initialise les tables de pre-calculs
   */
   private static void initTables() {
      BrmT1 = new int[1 << W1];
      Br1 = new int[1 << W1];
      BrmT2 = new int[1 << W2];
      Br2 = new int[1 << W2];
      BrmT3 = new int[1 << W3];
      Br3 = new int[1 << W3];

      for(int i = 0; i < Br1.length; i++) {
         BrmT1[i] = multiply(BrmT, i << Wsplit1, modQ);
         Br1[i] = multiply(Br, i << Wsplit1, modQ);
      }
      for(int i = 0; i < Br2.length; i++) {
         BrmT2[i] = multiply(BrmT, i << Wsplit2, modQ);
         Br2[i] = multiply(Br, i << Wsplit2, modQ);
      }
      for(int i = 0; i < Br3.length; i++) {
         BrmT3[i] = multiply(BrmT, i, modQ);
         Br3[i] = multiply(Br, i, modQ);
      }

   }

   private void advanceSeed(int[] seed, BitMatrix bm) {
      BitVector bv = new BitVector(seed, 800);

      bv = bm.multiply(bv);

      for(int i = 0; i < R; i++)
         seed[i] = bv.getInt(i);
   }

   private GenF2w32 (int i) {
      //unit vector (to build the state transition matrice)
      state = new int[R];
      for(int j = 0; j < R; j++)
         state[j] = 0;
      state[i / W] = 1 << (i % W);
      state_i = R - 1;
   }

   /**
    * Constructs a new stream.
    */
   public GenF2w32() {
      if (!initialised)
         initialisation();

      //stream = polyBase.createElem();
      //substream = polyBase.createElem();
      stream = new int[R];
      substream = new int[R];
      state = new int[R];

      for(int i = 0; i < R; i++)
         stream[i] = curr_stream[i];
      //stream.copyFrom(curr_stream);

      advanceSeed(curr_stream, Apz);
      //      curr_stream = curr_stream.multiply(jumpZ);

      resetStartStream();
   }

   /**
    * Constructs a new stream with the identifier `name` (used in the
    * `toString` method).
    *  @param name         name of the stream
    */
   public GenF2w32 (String name) {
      this();
      this.name = name;
   }

   /**
    * Sets the initial seed of the class `GenF2w2r32` to the 25 integers
    * of the vector `seed[0..24]`. This will be the initial seed of the
    * class for the next created stream. At least one of the integers must
    * be non-zero.
    *  @param seed         array of 25 elements representing the seed
    */
   public static void setPackageSeed (int seed[]) {
      if (!initialised)
         initialisation();
      if (seed.length < R)
         throw new IllegalArgumentException("Seed must contain " + R +
                                            "values.");
      boolean goodSeed = false;
      for(int i = 0; i < R; i++)
         if(seed[i] != 0)
            goodSeed = true;
      if(!goodSeed)
         throw new IllegalArgumentException("At least one part of the seed" +
                                            " must be non-zero.");

      for(int i = 0 ; i < R; i++)
         curr_stream[i] = seed[i];
      //curr_stream = polyBase.createElem(seed);
   }

   /**
    * This method is discouraged for normal use. Initializes the stream at
    * the beginning of a stream with the initial seed `seed[0..24]`. The
    * seed must satisfy the same conditions as in `setPackageSeed`. This
    * method only affects the specified stream; the others are not
    * modified. Hence after calling this method, the beginning of the
    * streams will no longer be spaced @f$Z@f$ values apart. For this
    * reason, this method should only be used in very exceptional cases;
    * proper use of the `reset...` methods and of the stream constructor
    * is preferable.
    *  @param seed         array of 25 elements representing the seed
    */
   public void setSeed (int seed[]) {
      if(seed.length != R)
         throw new IllegalArgumentException("Seed must contain " + R +
                                            "values.");
      boolean goodSeed = false;
      for(int i = 0; i < R; i++)
         if(seed[i] != 0)
            goodSeed = true;
      if(!goodSeed)
         throw new IllegalArgumentException("At least one part of the seed" +
                                            " must be non-zero.");

      for(int i = 0 ; i < R; i++)
         stream[i] = seed[i];
      //stream = polyBase.createElem(seed);
      resetStartStream();
   }

   /**
    * Returns the current state of the stream, represented as an array of
    * 25 integers.
    *  @return the current state of the stream
    */
   public int[] getState() {
      int res[] = new int[R];
      for(int i = 0; i < R; i++)
         res[i] = state[(state_i + i) % R];
      return res;
   }

   /**
    * Clones the current generator and return its copy.
    *  @return A deep copy of the current generator
    */
   public GenF2w32 clone() {
      GenF2w32 retour = null;
      retour = (GenF2w32)super.clone();
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
      //substream.copyFrom(stream);

      resetStartSubstream();
   }

   public void resetStartSubstream() {
      state_i = R - 1;
      for(int i = 0; i < R; i++)
         state[i] = substream[i];
      //substream.copyTo(state);
   }

   public void resetNextSubstream() {
      advanceSeed(substream, Apw);
      //      substream = substream.multiply(jumpW);

      resetStartSubstream();
   }

   public String toString()  {
      StringBuffer sb = new StringBuffer();

      sb.append("The state of the ");
      sb.append(name == null ? "GenF2w32" : name);
      sb.append(" is : {");
      for(int i = 0; i < R - 1; i++)
         sb.append(state[(state_i + i) % R] + ", ");
      sb.append(state[(state_i + R - 1) % R] + "}");

      return sb.toString();
   }

   /*
   //LFSR   (plus rapide d'environ 15% sur le LCG)
   protected double nextValue()
   {
      if(state_i >= R)
         state_i = 0;


      if(state_i < R - T)
         state[state_i] = ((BrmT1[(state[state_i + T] & MASK1) >>> Wsplit1] ^
                            BrmT2[(state[state_i + T] & MASK2) >>> Wsplit2] ^
                            BrmT3[(state[state_i + T] & MASK3)]) ^
                           (Br1[(state[state_i] & MASK1) >>> Wsplit1] ^
                            Br2[(state[state_i] & MASK2) >>> Wsplit2] ^
                            Br3[(state[state_i] & MASK3)]));
      else
         state[state_i] = ((BrmT1[(state[state_i +T-R] & MASK1) >>> Wsplit1] ^
                            BrmT2[(state[state_i +T-R] & MASK2) >>> Wsplit2] ^
                            BrmT3[(state[state_i +T-R] & MASK3)]) ^
                           (Br1[(state[state_i] & MASK1) >>> Wsplit1] ^
                            Br2[(state[state_i] & MASK2) >>> Wsplit2] ^
                            Br3[(state[state_i] & MASK3)]));

      long result = state[state_i++];

      return (result <= 0 ? result + 0x100000000L : result) * NORM;
   }
   */


   //LCG (plus lent)
   protected double nextValue() {
      if(state_i < 0)
         state_i = R - 1;

      if(state_i + T < R)
         state[state_i + T] ^= (BrmT1[(state[state_i] & MASK1) >>> Wsplit1] ^
                                BrmT2[(state[state_i] & MASK2) >>> Wsplit2] ^
                                BrmT3[(state[state_i] & MASK3)]);
      else
         state[state_i +T-R] ^= (BrmT1[(state[state_i] & MASK1) >>> Wsplit1] ^
                                 BrmT2[(state[state_i] & MASK2) >>> Wsplit2] ^
                                 BrmT3[(state[state_i] & MASK3)]);

      state[state_i] = (Br1[(state[state_i] & MASK1) >>> Wsplit1] ^
                        Br2[(state[state_i] & MASK2) >>> Wsplit2] ^
                        Br3[(state[state_i] & MASK3)]);

      long result = state[state_i--];

      return (result <= 0 ? result + 0x100000000L : result) * NORM;
   }

/**
 * This method is only meant to be used during the compilation process. It is
 * used to create the resource file the class need in order to run.
 */

   public static void main(String[] args) {
      if(args.length < 1) {
         System.err.println("Must provide the output file.");
         System.exit(1);
      }

      //computes the state transition matrices

      System.out.println("Creating the GenF2w32 state transition matrices.");

      //the state transition matrices
      BitMatrix STp0, STpw, STpz;

      BitVector[] bv = new BitVector[800];
      GenF2w32 gen;
      int[] vect = new int[R];

      for(int i = 0; i < 800; i++) {
         gen = new GenF2w32(i);

         gen.nextValue();
         for(int j = 0; j < R; j++)
            vect[j] = gen.state[(j + R - 1) % R];

         bv[i] = new BitVector(vect, 800);
      }

      STp0 = (new BitMatrix(bv)).transpose();

      STpw = STp0.power2e(w);
      STpz = STpw.power2e(v);


      try {
         FileOutputStream fos = new FileOutputStream(args[0]);
         ObjectOutputStream oos = new ObjectOutputStream(fos);
         oos.writeObject(STpw);
         oos.writeObject(STpz);
         oos.close();
      } catch(FileNotFoundException e) {
         System.err.println("Couldn't create " + args[0]);
         e.printStackTrace();
      } catch(IOException e) {
         e.printStackTrace();
      }

   }
}