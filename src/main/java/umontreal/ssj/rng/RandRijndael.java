/*
 * Class:        RandRijndael
 * Description:  RNG using the Rijndael block cipher algorithm (AES) with
                 key and block lengths of 128 bits
 * Environment:  Java
 * Software:     SSJ
 * Organization: DIRO, Universite de Montreal
 * @author
 * @since

 * SSJ is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, either version 3 of the License, or
 * any later version.

 * SSJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * A copy of the GNU General Public License is available at
   <a href="http://www.gnu.org/licenses">GPL licence site</a>.
 */
package umontreal.ssj.rng;

import java.io.Serializable;

/**
 * Implements a RNG using the Rijndael block cipher algorithm (AES) with key
 * and block lengths of 128 bits. A block of 128 bits is encrypted by the
 * Rijndael algorithm to generate 128 pseudo-random bits. Those bits are
 * split into four words of 32 bits which are returned successively by the
 * method `nextValue`. The unencrypted block is the state of the generator.
 * It is incremented by 1 at every four calls to `nextValue`. Thus, the
 * period is @f$2^{130}@f$ and jumping ahead is easy. The values of @f$V@f$,
 * @f$W@f$ and @f$Z@f$ are @f$2^{40}@f$, @f$2^{42}@f$ and @f$2^{82}@f$,
 * respectively (see  @ref RandomStream for their definition). Seeds/states
 * must be given as 16-dimensional vectors of bytes (8-bit integers). The
 * default initial seed is a vector filled with zeros.
 *
 * The Rijndael implementation used here is that of the *Cryptix Development
 * Team*, which can be found on the [Rijndael creators’
 * page](http://www.esat.kuleuven.ac.be/~rijmen/rijndael/)
 * [http://www.esat.kuleuven.ac.be/~rijmen/rijndael/](http://www.esat.kuleuven.ac.be/~rijmen/rijndael/).
 *
 * <div class="SSJ-bigskip"></div>
 */
public class RandRijndael extends RandomStreamBase {

   private static final long serialVersionUID = 70510L;
   //La date de modification a l'envers, lire 10/05/2007


   private static final int BLOCK_SIZE = 16;
   private static final int JUMP_STREAM = 10;
   private static final int JUMP_SUBSTREAM = 5;

   //actually a Object[] containing 2 int[][]
   private static Object key;

   private static byte[] curr_stream;
   private byte[] stream;
   private byte[] substream;

   private byte[] state;
   private byte[] output;
   private int outputPos;

   static
   {
      try {
         key = Rijndael_Algorithm.makeKey(new byte[]{1,2,3,4,5,6,7,8,
                                          9,10,11,12,13,14,15,16},
                                          BLOCK_SIZE);
      } catch(Exception e) {
         //pour que Java soit certain que la clef est initialisee
         key = new Object[0];
         e.printStackTrace();
         throw new RuntimeException("  cannot create RandRijndael key");
      }

      curr_stream = new byte[BLOCK_SIZE];
      for(int i = 0; i < BLOCK_SIZE; i++)
         curr_stream[i] = 0;
   }

   private static void iterate (byte[] b, int pos) {
      while((pos < b.length) && (++b[pos++] == 0));
   }

   /**
    * Constructs a new stream.
    */
   public RandRijndael() {
      stream = new byte[BLOCK_SIZE];
      substream = new byte[BLOCK_SIZE];

      state = new byte[BLOCK_SIZE];

      for(int i = 0; i < BLOCK_SIZE; i++)
         stream[i] = curr_stream[i];

      iterate(curr_stream, JUMP_STREAM);

      resetStartStream();
   }

   /**
    * Constructs a new stream with the identifier `name` (used in the
    * `toString` method).
    *  @param name         name of the stream
    */
   public RandRijndael (String name) {
      this();
      this.name = name;
   }

   /**
    * Sets the initial seed for the class `RandRijndael` to the 16 bytes
    * of the vector `seed[0..15]`. This will be the initial state (or
    * seed) of the next created stream. The default seed for the first
    * stream is @f$(0, 0, …, 0, 0)@f$.
    *  @param seed         array of 16 elements representing the seed
    */
   public static void setPackageSeed (byte seed[]) {
      if(seed.length != BLOCK_SIZE)
         throw new IllegalArgumentException("Seed must contain " +
                                            BLOCK_SIZE + " values");
      for(int i = 0; i < BLOCK_SIZE; i++)
         curr_stream[i] = seed[i];
   }

   /**
    * This method is discouraged for normal use. Initializes the stream at
    * the beginning of a stream with the initial seed `seed[0..15]`. This
    * method only affects the specified stream; the others are not
    * modified, so the beginning of the streams will not be spaced @f$Z@f$
    * values apart. For this reason, this method should only be used in
    * very exceptional cases; proper use of the `reset...` methods and of
    * the stream constructor is preferable.
    *  @param seed         array of 16 elements representing the seed
    */
   public void setSeed (byte seed[]) {
      if(seed.length != BLOCK_SIZE)
         throw new IllegalArgumentException("Seed must contain " +
                                            BLOCK_SIZE + " values");
      for(int i = 0; i < BLOCK_SIZE; i++)
         stream[i] = seed[i];
   }

   /**
    * Returns the current state of the stream, represented as an array of
    * four integers. It should be noted that each state of this generator
    * returns 4 successive values. The particular value of these 4 which
    * will be returned next is not given by this method.
    *  @return the current state of the stream
    */
   public byte[] getState() {
      byte[] stateCopy = new byte[BLOCK_SIZE];
      for(int i = 0; i < BLOCK_SIZE; i++)
         stateCopy[i] = state[i];
      return stateCopy;
   }

   /**
    * Clones the current generator and return its copy.
    *  @return A deep copy of the current generator
    */
   public RandRijndael clone() {
      RandRijndael retour = null;

      retour = (RandRijndael)super.clone();
      retour.stream = new byte[BLOCK_SIZE];
      retour.substream = new byte[BLOCK_SIZE];
      retour.state = new byte[BLOCK_SIZE];
      retour.output = new byte[output.length];
      for (int i = 0; i<BLOCK_SIZE; i++) {
         retour.stream[i] = stream[i];
         retour.substream[i] = substream[i];
         retour.state[i] = state[i];
      }
      for (int i=0; i<output.length; i++) {
         retour.output[i] = output[i];
      }

      return retour;
   }


   public void resetStartStream() {
      for(int i = 0; i < BLOCK_SIZE; i++)
         substream[i] = stream[i];

      resetStartSubstream();
   }

   public void resetStartSubstream() {
      for(int i = 0; i < BLOCK_SIZE; i++)
         state[i] = substream[i];
      nextOutput();
   }

   public void resetNextSubstream() {
      iterate(substream, JUMP_SUBSTREAM);
      resetStartSubstream();
   }

   public String toString() {
      StringBuffer sb = new StringBuffer();
      if(name == null)
         sb.append("The state of the RandRijn is : [");
      else
         sb.append("The state of the " + name + " is : [");

      for(int i = 0; i < BLOCK_SIZE - 1; i++)
         sb.append(state[i] + ", ");
      sb.append(state[BLOCK_SIZE - 1] + "]  ");

      sb.append("position : " + outputPos);

      return sb.toString();
   }

   private void nextOutput() {
      output = Rijndael_Algorithm.blockEncrypt(state, 0, key, BLOCK_SIZE);
      outputPos = 0;
      iterate(state,0);
   }

   protected double nextValue() {
      if(outputPos > BLOCK_SIZE - 4)
         nextOutput();


      long val = output[outputPos++] & 0xFF;
      val <<= 8;
      val |= output[outputPos++] & 0xFF;
      val <<= 8;
      val |= output[outputPos++] & 0xFF;
      val <<= 8;
      val |= output[outputPos++] & 0xFF;


      /*
      long val = ((output[outputPos] & 0xFF) << 24) |
         ((output[outputPos + 1] & 0xFF) << 16) |
         ((output[outputPos + 2] & 0xFF) << 8) |
         ((output[outputPos + 3] & 0xFF));
      outputPos += 4;
      */

      return ((double)val + 1) / 0x100000001L;
   }

   /*
   public static void main(String args[]) {
      int num = Integer.parseInt(args[0]);

      RandomStream rng = new RandRijn();

      for(int i = 0; i < num; i++) {
         rng.nextDouble();
         //System.out.println(rng.nextDouble());
      }

      System.out.println(rng.toString());
   }
   */

}