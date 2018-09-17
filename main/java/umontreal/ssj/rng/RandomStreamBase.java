/*
 * Class:        RandomStreamBase
 * Description:  Base class of all random number generators
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
 * This class provides a convenient foundation on which RNGs can be built. It
 * implements all the methods which do not depend directly on the generator
 * itself, but only on its output, which is to be defined by implementing the
 * `abstract` method `nextValue`. In the present class, all methods returning
 * random numbers directly or indirectly (<tt>nextDouble</tt>,
 * `nextArrayOfDouble`, `nextInt` and <tt>nextArrayOfInt</tt>) call
 * `nextValue`. Thus, to define a subclass that implements a RNG, it suffices
 * to implement `nextValue`, in addition to the `reset...` and `toString`
 * methods. Of course, the other methods may also be overridden for improved
 * efficiency.
 *
 * If the `nextValue` already generates numbers with a precision of 53-bits
 * or higher, then  #nextDouble can be overridden to improve the performance.
 * The mechanism for increasing the precision assumes that `nextValue`
 * returns at least 29 bits of precision, in which case the higher precision
 * numbers will have roughly 52 bits of precision. This mechanism was
 * designed primarily for RNGs that return numbers with around 30 to 32 bits
 * of precision.
 *
 *  @ref RandomStreamBase and its subclasses are implementing the
 * Serializable interface. Each class has a serial number wich represent the
 * class version. For instance `70510` means that the last change was the
 * `10th May 2007`.
 *
 * <div class="SSJ-bigskip"></div>
 */
public abstract class RandomStreamBase implements CloneableRandomStream,
                                                  Serializable {

   private static final long serialVersionUID = 70510L;
   //La date de modification a l'envers, lire 10/05/2007
   
   //constants
   protected static double invtwo24 = 5.9604644775390625e-8;  //2^(-24)
   private static double EPSILON = 5.5511151231257827e-17;    //2^(-54)

   protected String name = null;

   // prec53 keeps track if the precision has been increased or not.
   protected boolean prec53 = false;
   protected boolean anti = false;  // Deprecated.

   public abstract void resetStartStream();
   public abstract void resetStartSubstream();
   public abstract void resetNextSubstream();
   public abstract String toString();

   /**
    * After calling this method with `incp = true`, each call to the RNG
    * (direct or indirect) for this stream will return a uniform random
    * number with more bits of precision than what is returned by
    * `nextValue`, and will advance the state of the stream by 2 steps
    * instead of 1 (i.e., `nextValue` will be called twice for each random
    * number).
    *
    * More precisely, if `s` is a stream of a subclass of
    * `RandomStreamBase`, when the precision has been increased, the
    * instruction "<tt>u = s.nextDouble()</tt>", is equivalent to "<tt>u =
    * (s.nextValue() + s.nextValue()*fact) % 1.0</tt>" where the constant
    * `fact` is equal to @f$2^{-24}@f$. This also applies when calling
    * `nextDouble` indirectly (e.g., via `nextInt`, etc.). By default, or
    * if this method is called again with `incp = false`, each call to
    * `nextDouble` for this stream advances the state by 1 step and
    * returns the same number as `nextValue`.
    *  @param incp         if the generator will be set to high precision
    *                      mode
    */
   public void increasedPrecision (boolean incp) {
      prec53 = incp;
   }

   /**
    * This method should return the next random number (between 0 and 1)
    * from the current stream. If the stream is set to the high precision
    * mode (<tt>increasedPrecision(true)</tt> was called), then each call
    * to `nextDouble` will call `nextValue` twice, otherwise it will call
    * it only once.
    *  @return a number in the interval (0,1)
    */
   protected abstract double nextValue();

   /**
    * Returns a uniform random number between 0 and 1 from the stream. Its
    * behavior depends on the last call to  #increasedPrecision. The
    * generators programmed in SSJ never return the values 0 or 1.
    *  @return a number in the interval (0,1)
    */
   public double nextDouble() {
      double u = nextValue();
      if (prec53)
         u = (u + nextValue() * invtwo24) % 1.0 + EPSILON;
      if(anti)
         return 1.0 - u;
      else
         return u;
   }

   /**
    * Calls `nextDouble` `n` times to fill the array `u`.
    *  @param u            the array in which the numbers will be stored
    *  @param start        the first index of `u` to be used
    *  @param n            the number of random numbers to put in `u`
    */
   public void nextArrayOfDouble (double[] u, int start, int n) {
      if(u.length == 0)
         throw new NullPointerException("The array must be initialized.");
      if (u.length < n + start)
         throw new IndexOutOfBoundsException("The array is too small.");
      if(start < 0)
         throw new IndexOutOfBoundsException("Must start at a " +
                                             "non-negative index.");
      if(n < 0)
         throw new IllegalArgumentException("Must have a non-negative " +
                                            "number of elements.");

      for(int ii = start; ii < start + n; ii++)
         u[ii] = nextDouble();
   }

   /**
    * Calls `nextDouble` once to create one integer between `i` and `j`.
    * This method always uses the highest order bits of the random number.
    * It should be overridden if a faster implementation exists for the
    * specific generator.
    *  @param i            the smallest possible returned integer
    *  @param j            the largest possible returned integer
    *  @return a random integer between i and j
    */
   public int nextInt (int i, int j) {
      if(i > j)
         throw new IllegalArgumentException(i + " is larger than " +
                                            j + ".");
      // This works even for an interval [0, 2^31 - 1]. It would not with 
      // return i + (int)(nextDouble() * (j - i + 1));
      return i + (int)(nextDouble() * (j - i + 1.0));
   }

   /**
    * Calls `nextInt` `n` times to fill the array `u`. This method should
    * be overridden if a faster implementation exists for the specific
    * generator.
    *  @param i            the smallest possible integer to put in `u`
    *  @param j            the largest possible integer to put in `u`
    *  @param u            the array in which the numbers will be stored
    *  @param start        the first index of `u` to be used
    *  @param n            the number of random numbers to put in `u`
    */
   public void nextArrayOfInt (int i, int j, int[] u, int start, int n) {
      if(u == null)
         throw new NullPointerException("The array must be " +
                                        "initialized.");
      if (u.length < n + start)
         throw new IndexOutOfBoundsException("The array is too small.");
      if(start < 0)
         throw new IndexOutOfBoundsException("Must start at a " +
                                             "non-negative index.");
      if(n < 0)
         throw new IllegalArgumentException("Must have a non-negative " +
                                            "number of elements.");

      for(int ii = start; ii < start + n; ii++)
         u[ii] = nextInt(i,j);
   }

   /**
    * Use the `toString` method.
    */
   @Deprecated
   public String formatState() {
      return toString();
   }

   /**
    * Use the `toStringFull` method.
    */
   @Deprecated
   public String formatStateFull() {
      throw new UnsupportedOperationException (
            "   call the toStringFull() method instead.");
   }

   /**
    * Clones the current generator and return its copy.
    *  @return A deep copy of the current generator
    */
   public RandomStreamBase clone() {
      RandomStreamBase retour = null;
      try {
         retour = (RandomStreamBase)super.clone();
      }
      catch(CloneNotSupportedException cnse) {
         cnse.printStackTrace(System.err);
      }
      return retour;
   }
   

}