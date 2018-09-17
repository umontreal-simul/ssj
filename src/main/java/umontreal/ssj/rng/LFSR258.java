/*
 * Class:        LFSR258
 * Description:  64-bit composite linear feedback shift register proposed by L'Ecuyer
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
/*
import umontreal.ssj.util.BitVector;
import umontreal.ssj.util.BitMatrix;
*/

/**
 * Extends  @ref RandomStreamBase using a 64-bit composite linear feedback
 * shift register (LFSR) (or Tausworthe) RNG as defined in @cite rLEC96a,
 * @cite rTEZ91b&thinsp;. This generator is the `LFSR258` proposed in
 * @cite rLEC99a&thinsp;. It has five components combined by a bitwise xor.
 * Its period length is  @f$\rho\approx2^{258}@f$. The values of @f$V@f$,
 * @f$W@f$ and @f$Z@f$ are @f$2^{100}@f$, @f$2^{100}@f$ and @f$2^{200}@f$
 * respectively (see  @ref RandomStream for their definition). The seed of
 * the RNG, and the state of a stream at any given step, are five-dimensional
 * vectors of 64-bit integers. The default initial seed <sup title="In
 * previous versions, it was (1234567890, 1234567890, 1234567890, 1234567890,
 * 1234567890).">[1]</sup> of the RNG is (123456789123456789,
 * 123456789123456789, 123456789123456789, 123456789123456789,
 * 123456789123456789). The `nextValue` method returns numbers with 53 bits
 * of precision. This generator is fast for 64-bit machines.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class LFSR258 extends RandomStreamBase {

   private static final long serialVersionUID = 140406L;
   //La date de modification a l'envers, lire 06/04/2014
   //La date de modification a l'envers, lire 10/05/2007

    //private static final double NORM = 5.4210108624275221e-20;

    //equivalent a NORM = 1.0 / 0xFFFFFFFFFFFFF800L
    private static final double NORM = 0.5 / 0x7FFFFFFFFFFFFC00L;
    private static final double MAX = 0xFFFFFFFFFFFFF800L * NORM + 1.0;
    //MAX = plus grand double < 1.0

   private static final long GERME = 123456789123456789L;

    private long z0, z1, z2, z3, z4;       // l'etat

    //stream and substream variables :
    private long[] stream;
    private long[] substream;
    private static long[] curr_stream = {GERME, GERME, GERME, GERME, GERME};

    /**
     * Constructs a new stream.
     */
    public LFSR258() {
        name = null;

        stream = new long[5];
        substream = new long[5];

        for(int i = 0; i < 5; i++)
            stream[i] = curr_stream[i];

        resetStartStream();

        // Les operations qui suivent permettent de faire sauter en avant
        // de 2^200 iterations chacunes des composantes du generateur.
        // L'etat interne apres le saut est cependant legerement different
        // de celui apres 2^200 iterations puisqu'il ignore l'etat dans
        // lequel se retrouvent les premiers bits de chaque composantes,
        // puisqu'ils sont ignores dans la recurrence. L'etat redevient
        // identique a ce que l'on aurait avec des iterations normales
        // apres un appel a nextValue().

        long z, b;

        z = curr_stream[0] & 0xfffffffffffffffeL;
        b = z ^ (z << 1);
        z = (b >>> 58) ^ (b >>> 55) ^ (b >>> 46) ^ (b >>> 43) ^ (z << 5) ^
            (z << 8) ^ (z << 17) ^ (z << 20);
        curr_stream[0] = z;


        z = curr_stream[1] & 0xfffffffffffffe00L;
        b = z ^ (z << 24);
        z = (b >>> 54) ^ (b >>> 53) ^ (b >>> 52) ^ (b >>> 50) ^ (b >>> 49) ^
            (b >>> 48) ^ (b >>> 43) ^ (b >>> 41) ^ (b >>> 38) ^ (b >>> 37) ^
            (b >>> 30) ^ (b >>> 25) ^ (b >>> 24) ^ (b >>> 23) ^ (b >>> 19) ^
            (b >>> 16) ^ (b >>> 15) ^ (b >>> 14) ^ (b >>> 13) ^ (b >>> 11) ^
            (b >>> 8) ^ (b >>> 7) ^ (b >>> 5) ^ (b >>> 3) ^ (z << 0) ^
            (z << 2) ^ (z << 3) ^ (z << 6) ^ (z << 7) ^ (z << 8) ^ (z << 9) ^
            (z << 10) ^ (z << 11) ^ (z << 12) ^ (z << 13) ^ (z << 14) ^
            (z << 16) ^ (z << 18) ^ (z << 19) ^ (z << 21) ^ (z << 25) ^
            (z << 30) ^ (z << 31) ^ (z << 32) ^ (z << 36) ^ (z << 39) ^
            (z << 40) ^ (z << 41) ^ (z << 42) ^ (z << 44) ^ (z << 47) ^
            (z << 48) ^ (z << 50) ^ (z << 52);
        curr_stream[1] = z;


        z = curr_stream[2] & 0xfffffffffffff000L;
        b = z ^ (z << 3);
        z = (b >>> 50) ^ (b >>> 49) ^ (b >>> 46) ^ (b >>> 42) ^ (b >>> 40) ^
            (b >>> 39) ^ (b >>> 38) ^ (b >>> 37) ^ (b >>> 36) ^ (b >>> 32) ^
            (b >>> 29) ^ (b >>> 28) ^ (b >>> 27) ^ (b >>> 25) ^ (b >>> 23) ^
            (b >>> 20) ^ (b >>> 19) ^ (b >>> 15) ^ (b >>> 12) ^ (b >>> 11) ^
            (b >>> 2) ^ (z << 1) ^ (z << 2) ^ (z << 3) ^ (z << 6) ^ (z << 10) ^
            (z << 12) ^ (z << 13) ^ (z << 14) ^ (z << 15) ^ (z << 16) ^
            (z << 20) ^ (z << 23) ^ (z << 24) ^ (z << 25) ^ (z << 27) ^
            (z << 29) ^ (z << 32) ^ (z << 33) ^ (z << 37) ^ (z << 40) ^
            (z << 41) ^ (z << 50);
        curr_stream[2] = z;


        z = curr_stream[3] & 0xfffffffffffe0000L;
        b = z ^ (z << 5);
        z = (b >>> 46) ^ (b >>> 44) ^ (b >>> 42) ^ (b >>> 41) ^ (b >>> 40) ^
            (b >>> 38) ^ (b >>> 36) ^ (b >>> 32) ^ (b >>> 30) ^ (b >>> 25) ^
            (b >>> 18) ^ (b >>> 16) ^ (b >>> 15) ^ (b >>> 14) ^ (b >>> 12) ^
            (b >>> 11) ^ (b >>> 10) ^ (b >>> 9) ^ (b >>> 8) ^ (b >>> 6) ^
            (b >>> 5) ^ (b >>> 4) ^ (b >>> 3) ^ (b >>> 2) ^ (z << 2) ^
            (z << 5) ^ (z << 6) ^ (z << 7) ^ (z << 9) ^ (z << 11) ^ (z << 15) ^
            (z << 17) ^ (z << 22) ^ (z << 29) ^ (z << 31) ^ (z << 32) ^
            (z << 33) ^ (z << 35) ^ (z << 36) ^ (z << 37) ^ (z << 38) ^
            (z << 39) ^ (z << 41) ^ (z << 42) ^ (z << 43) ^ (z << 44) ^
            (z << 45);
        curr_stream[3] = z;


        z = curr_stream[4] & 0xffffffffff800000L;
        b = z ^ (z << 3);
        z = (b >>> 40) ^ (b >>> 29) ^ (b >>> 10) ^ (z << 1) ^ (z << 12) ^
            (z << 31);
        curr_stream[4] = z;

    }

    /**
     * Constructs a new stream with the identifier `name`.
     *  @param name         name of the stream
     */
    public LFSR258 (String name) {
        this();
        this.name = name;
    }

    /**
     * Sets the initial seed for the class `LFSR258` to the five integers
     * of array `seed[0..4]`. This will be the initial state of the next
     * created stream. The default seed <sup title="In previous versions,
     * it was (1234567890, 1234567890, 1234567890, 1234567890,
     * 1234567890).">[2]</sup> for the first stream is
     * (123456789123456789, 123456789123456789, 123456789123456789,
     * 123456789123456789, 123456789123456789). The first, second, third,
     * fourth and fifth integers of `seed` must be either negative, or
     * greater than or equal to 2, 512, 4096, 131072 and 8388608
     * respectively.
     *  @param seed         array of 5 elements representing the seed
     */
    public static void setPackageSeed (long seed[]) {
       checkSeed (seed);
       for(int i = 0; i < 5; i++)
          curr_stream[i] = seed[i];
    }

   private static void checkSeed (long seed[]) {
      if (seed.length < 5)
         throw new IllegalArgumentException("Seed must contain 5 values");
      if ((seed[0] >= 0 && seed[0] < 2)  ||
          (seed[1] >= 0 && seed[1] < 512)  ||
          (seed[2] >= 0 && seed[2] < 4096) ||
          (seed[3] >= 0 && seed[3] < 131072) ||
          (seed[4] >= 0 && seed[4] < 8388608))
         throw new IllegalArgumentException
         ("The seed elements must be either negative or greater than 1, 511, 4095, 131071 and 8388607 respectively");
   }

    /**
     * This method is discouraged for normal use. Initializes the stream
     * at the beginning of a stream with the initial seed `seed[0..4]`.
     * The seed must satisfy the same conditions as in `setPackageSeed`.
     * This method only affects the specified stream; the others are not
     * modified, so the beginning of the streams will not be spaced
     * @f$Z@f$ values apart. For this reason, this method should only be
     * used in very exceptional cases; proper use of the `reset...`
     * methods and of the stream constructor is preferable.
     *  @param seed         array of 5 elements representing the seed
     */
    public void setSeed (long seed[]) {
       checkSeed (seed);
       for(int i = 0; i < 5; i++)
           stream[i] = seed[i];
       resetStartStream();
    }

    /**
     * Returns the current state of the stream, represented as an array
     * of five integers.
     *  @return the current state of the stream
     */
    public long[] getState() {
        return new long[]{z0, z1, z2, z3, z4};
    }

   /**
    * Clones the current generator and return its copy.
    *  @return A deep copy of the current generator
    */
   public LFSR258 clone() {
      LFSR258 retour = null;
      retour = (LFSR258)super.clone();
      retour.stream = new long[5];
      retour.substream = new long[5];
      for (int i = 0; i<5; i++) {
         retour.substream[i] = substream[i];
         retour.stream[i] = stream[i];
      }
      return retour;
   }




    public void resetStartStream() {
        for(int i = 0; i < 5; i++)
            substream[i] = stream[i];
        resetStartSubstream();
    }

    public void resetStartSubstream() {
        z0 = substream[0];
        z1 = substream[1];
        z2 = substream[2];
        z3 = substream[3];
        z4 = substream[4];
    }

    public void resetNextSubstream() {
        // Les operations qui suivent permettent de faire sauter en avant
        // de 2^100 iterations chacunes des composantes du generateur.
        // L'etat interne apres le saut est cependant legerement different
        // de celui apres 2^100 iterations puisqu'il ignore l'etat dans
        // lequel se retrouvent les premiers bits de chaque composantes,
        // puisqu'ils sont ignores dans la recurrence. L'etat redevient
        // identique a ce que l'on aurait avec des iterations normales
        // apres un appel a nextValue().

       long z, b;

        z = substream[0] & 0xfffffffffffffffeL;
        b = z ^ (z << 1);
        z = (b >>> 61) ^ (b >>> 59) ^ (b >>> 58) ^ (b >>> 57) ^ (b >>> 51) ^
            (b >>> 47) ^ (b >>> 46) ^ (b >>> 45) ^ (b >>> 43) ^ (b >>> 39) ^
            (b >>> 30) ^ (b >>> 29) ^ (b >>> 23) ^ (b >>> 15) ^ (z << 2) ^
            (z << 4) ^ (z << 5) ^ (z << 6) ^ (z << 12) ^ (z << 16) ^
            (z << 17) ^ (z << 18) ^ (z << 20) ^ (z << 24) ^ (z << 33) ^
            (z << 34) ^ (z << 40) ^ (z << 48);
        substream[0] = z;


        z = substream[1] & 0xfffffffffffffe00L;
        b = z ^ (z << 24);
        z = (b >>> 52) ^ (b >>> 50) ^ (b >>> 49) ^ (b >>> 46) ^ (b >>> 43) ^
            (b >>> 40) ^ (b >>> 37) ^ (b >>> 34) ^ (b >>> 30) ^ (b >>> 28) ^
            (b >>> 26) ^ (b >>> 25) ^ (b >>> 23) ^ (b >>> 21) ^ (b >>> 20) ^
            (b >>> 19) ^ (b >>> 17) ^ (b >>> 15) ^ (b >>> 13) ^ (b >>> 12) ^
            (b >>> 10) ^ (b >>> 8) ^ (b >>> 7) ^ (b >>> 6) ^ (b >>> 2) ^
            (z << 1) ^ (z << 4) ^ (z << 6) ^ (z << 7) ^ (z << 11) ^ (z << 14) ^
            (z << 15) ^ (z << 16) ^ (z << 17) ^ (z << 21) ^ (z << 22) ^
            (z << 25) ^ (z << 27) ^ (z << 29) ^ (z << 30) ^ (z << 32) ^
            (z << 34) ^ (z << 35) ^ (z << 36) ^ (z << 38) ^ (z << 40) ^
            (z << 42) ^ (z << 43) ^ (z << 45) ^ (z << 47) ^ (z << 48) ^
            (z << 49) ^ (z << 53);
        substream[1] = z;


        z = substream[2] & 0xfffffffffffff000L;
        b = z ^ (z << 3);
        z = (b >>> 49) ^ (b >>> 45) ^ (b >>> 41) ^ (b >>> 40) ^ (b >>> 32) ^
            (b >>> 27) ^ (b >>> 23) ^ (b >>> 14) ^ (b >>> 1) ^ (z << 2) ^
            (z << 3) ^ (z << 7) ^ (z << 11) ^ (z << 12) ^ (z << 20) ^
            (z << 25) ^ (z << 29) ^ (z << 38) ^ (z << 51);
        substream[2] = z;



        z = substream[3] & 0xfffffffffffe0000L;
        b = z ^ (z << 5);
        z = (b >>> 45) ^ (b >>> 32) ^ (b >>> 27) ^ (b >>> 22) ^ (b >>> 17) ^
            (b >>> 13) ^ (b >>> 12) ^ (b >>> 7) ^ (b >>> 3) ^ (b >>> 2) ^
            (z << 3) ^ (z << 15) ^ (z << 20) ^ (z << 25) ^ (z << 30) ^
            (z << 34) ^ (z << 35) ^ (z << 40) ^ (z << 44) ^ (z << 45);
        substream[3] = z;


        z = substream[4] & 0xffffffffff800000L;
        b = z ^ (z << 3);
        z = (b >>> 40) ^ (b >>> 39) ^ (b >>> 38) ^ (b >>> 37) ^ (b >>> 35) ^
            (b >>> 34) ^ (b >>> 31) ^ (b >>> 30) ^ (b >>> 29) ^ (b >>> 28) ^
            (b >>> 27) ^ (b >>> 26) ^ (b >>> 24) ^ (b >>> 23) ^ (b >>> 21) ^
            (b >>> 20) ^ (b >>> 18) ^ (b >>> 15) ^ (b >>> 12) ^ (b >>> 10) ^
            (b >>> 9) ^ (b >>> 7) ^ (b >>> 6) ^ (b >>> 5) ^ (b >>> 4) ^
            (b >>> 3) ^ (z << 1) ^ (z << 2) ^ (z << 3) ^ (z << 4) ^ (z << 6) ^
            (z << 7) ^ (z << 10) ^ (z << 11) ^ (z << 12) ^ (z << 13) ^
            (z << 14) ^ (z << 15) ^ (z << 17) ^ (z << 18) ^ (z << 20) ^
            (z << 21) ^ (z << 23) ^ (z << 26) ^ (z << 29) ^ (z << 31) ^
            (z << 32) ^ (z << 34) ^ (z << 35) ^ (z << 36) ^ (z << 37) ^
            (z << 38);
        substream[4] = z;

        resetStartSubstream();
    }


    public String toString() {
        if (name == null)
            return "The state of the LFSR258 is: " +
                   z0 + "L, " + z1 + "L, " + z2 + "L, " + z3 + "L, " + z4 + "L";
        else
            return "The state of " + name + " is: " +
                   z0 + "L, " + z1 + "L, " + z2 + "L, " + z3 + "L, " + z4 + "L";
    }


    private long nextNumber() {
       long b;
       b  = (((z0 <<  1) ^ z0) >>> 53);
       z0 = (((z0 & 0xFFFFFFFFFFFFFFFEL) << 10) ^ b);
       b  = (((z1 << 24) ^ z1) >>> 50);
       z1 = (((z1 & 0xFFFFFFFFFFFFFE00L) <<  5) ^ b);
       b  = (((z2 <<  3) ^ z2) >>> 23);
       z2 = (((z2 & 0xFFFFFFFFFFFFF000L) << 29) ^ b);
       b  = (((z3 <<  5) ^ z3) >>> 24);
       z3 = (((z3 & 0xFFFFFFFFFFFE0000L) << 23) ^ b);
       b  = (((z4 <<  3) ^ z4) >>> 33);
       z4 = (((z4 & 0xFFFFFFFFFF800000L) <<  8) ^ b);

       return (z0 ^ z1 ^ z2 ^ z3 ^ z4);
    }


    protected double nextValue() {

        long res = nextNumber();
        if (res <= 0)
            return (res * NORM + MAX);
        else
            return res * NORM;
    }


   public int nextInt (int i, int j) {
      if (i > j)
          throw new IllegalArgumentException(i + " is larger than " + j + ".");
      long d = j-i+1;
      long q = 0x4000000000000000L / d;
      long r = 0x4000000000000000L % d;
      long res;

      do {
         res = nextNumber() >>> 2;
      } while (res >= 0x4000000000000000L - r);

      return i + (int) (res / q);
   }


    /*
     Methodes qui permettent de generer les series de shifts des
     methodes de sauts en avant.


     Preuve que la fonction resultante est bien equivalente a la matrice
     fournie :

     Soit M la matrice de transition pour une iteration simple du Tausworthe
     d'une seule composante. (Il est a noter que chaque appel de
     nextValue fait s iterations simples pour chaque composantes.)

     On prend I la matrice identite. (I = M^0)
     Chaque colonne de I represente un vecteur d'etat pour le generateur.
     Puisque les (64 - k) bits les moins significatifs du vecteur d'etat
     ne sont pas repris par la recurrence a long terme, c'est-a-dire que
     la correlation entre leurs valeurs initiales et la valeur de l'etat
     (64 - k) iterations simples plus loin est nulle. Donc, a partir de
     M^(64 - k), les (64 - k) premieres colonnes de la matrice sont nulles.
     Puisque l'on traite qu'avec de tres larges exposants, on peut
     supposer que ces colonnes sont nulles, ce qui est equivalent a
     mettre a 0 les (64 - k) premiers bits.

     Les colonnes qui suivent dans I, jusqu'aux q dernieres, sont toutes
     a une iteration simple d'ecart. Lorsque la matrice I va etre multipliee
     un certain nombre de fois par M pour donner M^x, puisque cette operation
     est equivalente a faire x iterations simples, ces colonnes resteront
     a une iteration simple d'ecart. Puisque l'effet d'une iteration simple
     consiste simplement de faire un shift deplacant tous les bits d'une
     position ainsi que de remplir le bit vide ainsi cree par une nouvelle
     valeur, alors les differentes colonnes ne different que ce shift et
     leurs extremites. La consequence est que, pour ces colonnes, tous les
     bits sur la meme diagonale sont egaux. Il suffit alors de prendre les
     valeurs aux debuts des diagonales (dernieres rangee et colonne) et
     d'en ressortir les shifts puisqu'un shift est equivalent a une
     translation des bits de la matrice identite, donc equivalent a une
     diagonale dans la matrice.

     Le meme argument peut etre applique aux q dernieres colonnes de la
     matrice. Il y a cependant un fracture entre la q-ieme derniere colonne
     et la (q+1)-ieme derniere. Cet ecart vient du fait que, dans I, ces
     deux colonnes ne sont pas a une iteration simple d'ecart. Puisque
     la recurrence qui definit l'iteration simple est
     x_n = x_(n-(k-q)) ^ x_(n-k), le nouveau bit qui s'ajoute au debut
     de la partie significative (k premiers bits) du vecteur d'etat est
     la somme binaire du dernier bit et de (q+1)-ieme dernier bit. Donc,
     la (q+1)-ieme derniere colonne est equivalente a la combinaison
     lineaire de la derniere et du resultat d'une iteration simple sur la
     q-ieme derniere colonne. La consequence est premierement que
     toutes diagonales (donc shift) partant de la derniere colonne continue
     apres la (q+1)-ieme colonne. Secondement, pour chacune de ces diagonales,
     une autre diagonale commence a la meme rangee, mais a la (q+1)-ieme
     colonne, ce qui est equivalent a un autre shift, mais avec un masque
     pour couper les q derniers bits du shift. La combinaison de ces deux
     diagonales est b = z ^ (z << q).



     private static void analyseMat(BitMatrix bm, int decal, int signif) {
         //note : decal = q, signif = k (selon la notation de l'article)

         System.out.println("z = z & 0x" +
                            Long.toHexString(-1 << (64 - signif)) + "L;");
         System.out.println("b = z ^ (z << " + decal + ");");

         StringBuffer sb = new StringBuffer("z =");
         for(int i = (64 - signif); i < 63; i++)
             if(bm.getBool(i, 63))
                 sb.append(" (b >>> " + (63 - i) + ") ^");

         for(int i = 0; i < 64; i++)
             if(bm.getBool(63, 63 - i))
                 sb.append(" (z << " + i + ") ^");
         sb.setCharAt(sb.length() - 1, ';');


         System.out.println(sb);
     }


     public static void main(String[] args) {
         BitMatrix bm;
         BitVector[] bv = new BitVector[64];
         LFSR258 gen = new LFSR258();

         for(int i = 0; i < 64; i++) {
             gen.z0 = 1L << i;
             gen.nextValue();
             bv[i] = new BitVector(new int[]{(int)gen.z0,
                                             (int)(gen.z0 >>> 32)});
         }

         bm = (new BitMatrix(bv)).transpose();

         int W = 100;
         int Z = 200;

         BitMatrix bmPw = bm.power2e(W);
         BitMatrix bmPz = bm.power2e(Z);

         analyseMat(bmPw, 1, 63);
         analyseMat(bmPz, 1, 63);
     }
    */

}