/*
 * Class:        MRG32k3aL
 * Description:  Long version of MRG32k3a
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

import umontreal.ssj.util.ArithmeticMod;
import umontreal.ssj.util.PrintfFormat;
import java.io.Serializable;

/**
 * The same generator as  @ref MRG32k3a, except here it is implemented with
 * type `long` instead of `double`. (See  @ref MRG32k3a for more
 * information.)
 *
 * <div class="SSJ-bigskip"></div>
 */
public class MRG32k3aL extends RandomStreamBase {

   private static final long serialVersionUID = 70510L;
   //La date de modification a l'envers, lire 10/05/2007

   // Private constants   %%%%%%%%%%%%%%%%%%%%%%%%%%

   private static final long   m1     = 4294967087L;
   private static final long   m2     = 4294944443L;
   private static final long   a12    =  1403580L;
   private static final long   a13n   =   810728L;
   private static final long   a21    =   527612L;
   private static final long   a23n   =  1370589L;
   private static final long   two17  =   131072L;
   private static final long   two53  =  9007199254740992L;
   private static final double invtwo24 = 5.9604644775390625e-8;
   private static final double norm   = 2.328306549295727688e-10;
   //    private static final double norm   = 1.0 / (m1 + 1.0);


   /*  Unused
   private static final double InvA1[][] = {   // Inverse of A1p0
     { 184888585.0, 0.0, 1945170933.0 },
     {         1.0, 0.0,          0.0 },
     {         0.0, 1.0,          0.0 }
     };
   private static final double InvA2[][] = {   // Inverse of A2p0
     { 0.0, 360363334.0, 4225571728.0 },
     { 1.0,         0.0,          0.0 },
     { 0.0,         1.0,          0.0 }
     };
   */

   private static final long A1p0[][]  =  {
            {       0L,       1L,      0L },
            {       0L,       0L,      1L },
            { -810728L, 1403580L,      0L }
         };
   private static final long A2p0[][]  =  {
            {        0L,   1L,         0L },
            {        0L,   0L,         1L },
            { -1370589L,   0L,    527612L }
         };
   private static final long A1p76[][] = {
       {   82758667L, 1871391091L, 4127413238L },
       { 3672831523L,   69195019L, 1871391091L },
       { 3672091415L, 3528743235L,   69195019L }
                                           };
   private static final long A2p76[][] = {
       { 1511326704L, 3759209742L, 1610795712L },
       { 4292754251L, 1511326704L, 3889917532L },
       { 3859662829L, 4292754251L, 3708466080L }
                                           };
   private static final long A1p127[][] = {
            {    2427906178L, 3580155704L,  949770784L },
            {     226153695L, 1230515664L, 3580155704L },
            {    1988835001L,  986791581L, 1230515664L }
         };
   private static final long A2p127[][] = {
            {    1464411153L,  277697599L, 1610723613L },
            {      32183930L, 1464411153L, 1022607788L },
            {    2824425944L,   32183930L, 2093834863L }
         };


   // Private variables for each stream   %%%%%%%%%%%%%%%%%%%%%%%%

   // Default seed of the package for the first stream
   private static long nextSeed[] = {12345, 12345, 12345,
                                       12345, 12345, 12345};
   private long Cg0, Cg1, Cg2, Cg3, Cg4, Cg5;
   private long Bg[] = new long[6];
   private long Ig[] = new long[6];
   // The arrays Cg, Bg and Ig contain the current state,
   // the starting point of the current substream,
   // and the starting point of the stream, respectively.


   //multiply the first half of v by A with a modulo of m1
   //and the second half by B with a modulo of m2
   private static void multMatVect(long[] v, long[][] A, long m1,
                                   long[][] B, long m2) {
      long[] vv = new long[3];
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
   public MRG32k3aL() {
      name = null;
      anti = false;
      prec53 = false;
      for(int i = 0; i < 6; i++)
         Ig[i] = nextSeed[i];
      resetStartStream();
      multMatVect(nextSeed, A1p127, m1, A2p127, m2);
   }
   public MRG32k3aL (String name) {
      this();
      this.name = name;
   }

   /**
    * @param name         name of the stream
    *
    *  See the description of the same methods in class  @ref MRG32k3a.
    */
   public static void setPackageSeed (long seed[]) {
      // Must use long because there is no unsigned int type.
      validateSeed (seed);
      for (int i = 0; i < 6;  ++i)
         nextSeed[i] = seed[i];
   }

   /**
    * @param seed         array of 6 elements representing the seed
    */
   public void setSeed (long seed[]) {
      // Must use long because there is no unsigned int type.
      validateSeed (seed);
      for (int i = 0; i < 6;  ++i)
         Ig[i] = seed[i];
      resetStartStream();
   }

/**
 * @param seed         array of 6 integers representing the new seed
 */

   public void resetStartStream()  {
      for (int i = 0; i < 6;  ++i)
         Bg[i] = Ig[i];
      resetStartSubstream();
   }

   public void resetStartSubstream()  {
      Cg0 = Bg[0];
      Cg1 = Bg[1];
      Cg2 = Bg[2];
      Cg3 = Bg[3];
      Cg4 = Bg[4];
      Cg5 = Bg[5];
   }

   public void resetNextSubstream()  {
      multMatVect(Bg, A1p76, m1, A2p76, m2);
      resetStartSubstream();
   }
   public long[] getState() {
      return new long[]{Cg0, Cg1, Cg2, Cg3, Cg4, Cg5};
   }

   /**
    * @return the current state of the generator
    */
   public String toString() {
      PrintfFormat str = new PrintfFormat();

      str.append ("The current state of the MRG32k3aL");
      if (name != null && name.length() > 0)
         str.append (" " + name);
      str.append (":" + PrintfFormat.NEWLINE + "   Cg = { ");
      str.append ( Cg0 + ", ");
      str.append ( Cg1 + ", ");
      str.append ( Cg2 + ", ");
      str.append ( Cg3 + ", ");
      str.append ( Cg4 + ", ");
      str.append ( Cg5 + " }" + PrintfFormat.NEWLINE +
           PrintfFormat.NEWLINE);

      return str.toString();
   }

   /**
    * @return the state of the generator, formated as a string
    */
   public String toStringFull() {
      PrintfFormat str = new PrintfFormat();
      str.append ("The MRG32k3aL stream");
      if (name != null && name.length() > 0)
         str.append (" " + name);
      str.append (":" + PrintfFormat.NEWLINE + "   anti = " +
         (anti ? "true" : "false")).append(PrintfFormat.NEWLINE);
      str.append ("   prec53 = " + (prec53 ? "true" : "false")).append(PrintfFormat.NEWLINE);

      str.append ("   Ig = { ");
      for (int i = 0; i < 5; i++)
         str.append ( Ig[i] + ", ");
      str.append ( Ig[5] + " }" + PrintfFormat.NEWLINE);

      str.append ("   Bg = { ");
      for (int i = 0; i < 5; i++)
         str.append ( Bg[i] + ", ");
      str.append ( Bg[5] + " }" + PrintfFormat.NEWLINE);

      str.append ("   Cg = { ");
      str.append ( Cg0 + ", ");
      str.append ( Cg1 + ", ");
      str.append ( Cg2 + ", ");
      str.append ( Cg3 + ", ");
      str.append ( Cg4 + ", ");
      str.append ( Cg5 + " }" + PrintfFormat.NEWLINE +
          PrintfFormat.NEWLINE);

      return str.toString();
   }

   /**
    * @return the detailed state of the generator, formatted as a string
    */
   public MRG32k3aL clone() {
      MRG32k3aL retour = null;

      retour = (MRG32k3aL)super.clone();
      retour.Bg = new long[6];
      retour.Ig = new long[6];
      for (int i = 0; i<6; i++) {
         retour.Bg[i] = Bg[i];
         retour.Ig[i] = Ig[i];
      }
      return retour;
   }

/**
 * @return A deep copy of the current generator
 */

   protected double nextValue() {
      int k;
      long p1, p2;

      /* Component 1 */
      p1 = (a12 * Cg1 - a13n * Cg0) % m1;
      if (p1 < 0)
         p1 += m1;
      Cg0 = Cg1;
      Cg1 = Cg2;
      Cg2 = p1;

      /* Component 2 */
      p2 = (a21 * Cg5 - a23n * Cg3) % m2;
      if (p2 < 0)
         p2 += m2;
      Cg3 = Cg4;
      Cg4 = Cg5;
      Cg5 = p2;

      /* Combination */
      return (double)((p1 > p2) ? (p1 - p2) * norm : (p1 - p2 + m1) * norm);
   }


   private static void validateSeed (long seed[]) {
      if (seed.length < 6)
         throw new IllegalArgumentException ("Seed must contain 6 values");
      if (seed[0] == 0 && seed[1] == 0 && seed[2] == 0)
         throw new IllegalArgumentException
            ("The first 3 values must not be 0");
      if (seed[3] == 0 && seed[4] == 0 && seed[5] == 0)
         throw new IllegalArgumentException
            ("The last 3 values must not be 0");
      final long m1 = 4294967087L;
      if (seed[0] >= m1 || seed[1] >= m1 || seed[2] >= m1)
         throw new IllegalArgumentException
            ("The first 3 values must be less than " + m1);
      final long m2 = 4294944443L;
      if (seed[3] >= m2 || seed[4] >= m2 || seed[5] >= m2)
         throw new IllegalArgumentException
            ("The last 3 values must be less than " + m2);
   }
}