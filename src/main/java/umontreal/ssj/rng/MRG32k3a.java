/*
 * Class:        MRG32k3a
 * Description:  combined multiple recursive generator proposed by L'Ecuyer
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
 * Extends the abstract class  @ref RandomStreamBase by using as a backbone
 * (or main) generator the combined multiple recursive generator (CMRG)
 * `MRG32k3a` proposed by L’Ecuyer @cite rLEC99b&thinsp;, implemented in
 * 64-bit floating-point arithmetic. This backbone generator has a period
 * length of  @f$\rho\approx2^{191}@f$. The values of @f$V@f$, @f$W@f$, and
 * @f$Z@f$ are @f$2^{51}@f$, @f$2^{76}@f$, and @f$2^{127}@f$, respectively.
 * (See  @ref RandomStream for their definition.) The seed of the RNG, and
 * the state of a stream at any given step, are six-dimensional vectors of
 * 32-bit integers, stored in `double`. The default initial seed of the RNG
 * is @f$(12345, 12345, 12345, 12345, 12345, 12345)@f$.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class MRG32k3a extends RandomStreamBase {

   private static final long serialVersionUID = 70510L;
   //La date de modification a l'envers, lire 10/05/2007

   // Private constants   %%%%%%%%%%%%%%%%%%%%%%%%%%

   private static final double m1     = 4294967087.0;
   private static final double m2     = 4294944443.0;
   private static final double a12    =  1403580.0;
   private static final double a13n   =   810728.0;
   private static final double a21    =   527612.0;
   private static final double a23n   =  1370589.0;
   private static final double two17    =  131072.0;
   private static final double two53    =  9007199254740992.0;
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

   private static final double A1p0[][]  =  {
            {       0.0,       1.0,      0.0 },
            {       0.0,       0.0,      1.0 },
            { -810728.0, 1403580.0,      0.0 }
         };
   private static final double A2p0[][]  =  {
            {        0.0,   1.0,         0.0 },
            {        0.0,   0.0,         1.0 },
            { -1370589.0,   0.0,    527612.0 }
         };
   private static final double A1p76[][] = {
       { 82758667.0, 1871391091.0, 4127413238.0 },
       { 3672831523.0,   69195019.0, 1871391091.0 },
       { 3672091415.0, 3528743235.0,   69195019.0 }
                                           };
   private static final double A2p76[][] = {
       { 1511326704.0, 3759209742.0, 1610795712.0 },
       { 4292754251.0, 1511326704.0, 3889917532.0 },
       { 3859662829.0, 4292754251.0, 3708466080.0 }
                                           };
   private static final double A1p127[][] = {
            {    2427906178.0, 3580155704.0,  949770784.0 },
            {     226153695.0, 1230515664.0, 3580155704.0 },
            {    1988835001.0,  986791581.0, 1230515664.0 }
         };
   private static final double A2p127[][] = {
            {    1464411153.0,  277697599.0, 1610723613.0 },
            {      32183930.0, 1464411153.0, 1022607788.0 },
            {    2824425944.0,   32183930.0, 2093834863.0 }
         };


   // Private variables for each stream   %%%%%%%%%%%%%%%%%%%%%%%%

   // Default seed of the package for the first stream
   private static double nextSeed[] = {12345, 12345, 12345,
                                       12345, 12345, 12345};
   private double Cg0, Cg1, Cg2, Cg3, Cg4, Cg5;
   private double Bg[] = new double[6];
   private double Ig[] = new double[6];
   // The arrays Cg, Bg and Ig contain the current state,
   // the starting point of the current substream,
   // and the starting point of the stream, respectively.


   //multiply the first half of v by A with a modulo of m1
   //and the second half by B with a modulo of m2
   private static void multMatVect(double[] v, double[][] A, double m1,
                                   double[][] B, double m2) {
      double[] vv = new double[3];
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
    * Constructs a new stream, initializes its seed @f$I_g@f$, sets
    * @f$B_g@f$ and @f$C_g@f$ equal to @f$I_g@f$, and sets its antithetic
    * switch to `false`. The seed @f$I_g@f$ is equal to the initial seed
    * of the package given by  #setPackageSeed(long[]) if this is the
    * first stream created, otherwise it is @f$Z@f$ steps ahead of that of
    * the stream most recently created in this class.
    */
   public MRG32k3a() {
      name = null;
      anti = false;
      prec53 = false;
      for(int i = 0; i < 6; i++)
         Ig[i] = nextSeed[i];
      resetStartStream();
      multMatVect(nextSeed, A1p127, m1, A2p127, m2);
   }

   /**
    * Constructs a new stream with an identifier `name` (used when
    * printing the stream state).
    *  @param name         name of the stream
    */
   public MRG32k3a (String name) {
      this();
      this.name = name;
   }

   /**
    * Sets the initial seed for the class `MRG32k3a` to the six integers
    * in the vector `seed[0..5]`. This will be the seed (initial state) of
    * the first stream. If this method is not called, the default initial
    * seed is @f$(12345, 12345, 12345, 12345, 12345, 12345)@f$. If it is
    * called, the first 3 values of the seed must all be less than @f$m_1
    * = 4294967087@f$, and not all 0; and the last 3 values must all be
    * less than @f$m_2 = 4294944443@f$, and not all 0.
    *  @param seed         array of 6 elements representing the seed
    */
   public static void setPackageSeed (long seed[]) {
      // Must use long because there is no unsigned int type.
      validateSeed (seed);
      for (int i = 0; i < 6;  ++i)
         nextSeed[i] = seed[i];
   }

   /**
    * Sets the initial seed @f$I_g@f$ of this stream to the vector
    * `seed[0..5]`. This vector must satisfy the same conditions as in
    * `setPackageSeed`. The stream is then reset to this initial seed. The
    * states and seeds of the other streams are not modified. As a result,
    * after calling this method, the initial seeds of the streams are no
    * longer spaced @f$Z@f$ values apart. For this reason, *this method
    * should be used only in very exceptional situations* (I have never
    * used it myself!); proper use of `reset...` and of the stream
    * constructor is preferable.
    *  @param seed         array of 6 integers representing the new seed
    */
   public void setSeed (long seed[]) {
      // Must use long because there is no unsigned int type.
      validateSeed (seed);
      for (int i = 0; i < 6;  ++i)
         Ig[i] = seed[i];
      resetStartStream();
   }

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

/**
 * Returns the current state @f$C_g@f$ of this stream. This is a vector of 6
 * integers. This method is convenient if we want to save the state for
 * subsequent use.
 *  @return the current state of the generator
 */
public long[] getState() {
      return new long[] {(long)Cg0, (long)Cg1, (long)Cg2,
                         (long)Cg3, (long)Cg4, (long)Cg5};
   }

   /**
    * Returns a string containing the name and the current state @f$C_g@f$
    * of this stream.
    *  @return the state of the generator, formated as a string
    */
   public String toString() {
      PrintfFormat str = new PrintfFormat();

      str.append ("The current state of the MRG32k3a");
      if (name != null && name.length() > 0)
         str.append (" " + name);
      str.append (":" + PrintfFormat.NEWLINE + "   Cg = { ");
      str.append ((long) Cg0 + ", ");
      str.append ((long) Cg1 + ", ");
      str.append ((long) Cg2 + ", ");
      str.append ((long) Cg3 + ", ");
      str.append ((long) Cg4 + ", ");
      str.append ((long) Cg5 + " }" + PrintfFormat.NEWLINE +
          PrintfFormat.NEWLINE);

      return str.toString();
   }

   /**
    * Returns a string containing the name of this stream and the values
    * of all its internal variables.
    *  @return the detailed state of the generator, formatted as a string
    */
   public String toStringFull() {
      PrintfFormat str = new PrintfFormat();
      str.append ("The MRG32k3a stream");
      if (name != null && name.length() > 0)
         str.append (" " + name);
      str.append (":" + PrintfFormat.NEWLINE + "   anti = " +
         (anti ? "true" : "false")).append(PrintfFormat.NEWLINE);
      str.append ("   prec53 = " + (prec53 ? "true" : "false")).append(PrintfFormat.NEWLINE);

      str.append ("   Ig = { ");
      for (int i = 0; i < 5; i++)
         str.append ((long) Ig[i] + ", ");
      str.append ((long) Ig[5] + " }" + PrintfFormat.NEWLINE);

      str.append ("   Bg = { ");
      for (int i = 0; i < 5; i++)
         str.append ((long) Bg[i] + ", ");
      str.append ((long) Bg[5] + " }" + PrintfFormat.NEWLINE);

      str.append ("   Cg = { ");
      str.append ((long) Cg0 + ", ");
      str.append ((long) Cg1 + ", ");
      str.append ((long) Cg2 + ", ");
      str.append ((long) Cg3 + ", ");
      str.append ((long) Cg4 + ", ");
      str.append ((long) Cg5 + " }" + PrintfFormat.NEWLINE +
          PrintfFormat.NEWLINE);

      return str.toString();
   }

   /**
    * Clones the current generator and return its copy.
    *  @return A deep copy of the current generator
    */
   public MRG32k3a clone() {
      MRG32k3a retour = null;

      retour = (MRG32k3a)super.clone();
      retour.Bg = new double[6];
      retour.Ig = new double[6];
      for (int i = 0; i<6; i++) {
         retour.Bg[i] = Bg[i];
         retour.Ig[i] = Ig[i];
      }
      return retour;
   }


   protected double nextValue() {
      int k;
      double p1, p2;
      /* Component 1 */
      p1 = a12 * Cg1 - a13n * Cg0;
      k = (int)(p1 / m1);
      p1 -= k * m1;
      if (p1 < 0.0)
         p1 += m1;
      Cg0 = Cg1;
      Cg1 = Cg2;
      Cg2 = p1;
      /* Component 2 */
      p2 = a21 * Cg5 - a23n * Cg3;
      k  = (int)(p2 / m2);
      p2 -= k * m2;
      if (p2 < 0.0)
         p2 += m2;
      Cg3 = Cg4;
      Cg4 = Cg5;
      Cg5 = p2;
      /* Combination */
      return ((p1 > p2) ? (p1 - p2) * norm : (p1 - p2 + m1) * norm);
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