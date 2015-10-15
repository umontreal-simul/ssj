/*
 * Class:        BakerTransformedStream
 * Description:  container class permits one to apply the baker's 
                 transformation to the output of any RandomStream
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
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

/**
 * This container class permits one to apply the baker’s transformation to
 * the output of any  @ref RandomStream. It transforms each @f$u \in[0,1]@f$
 * into @f$2u@f$ if @f$u \le1/2@f$ and @f$2(1-u)@f$ if @f$u > 1/2@f$. The
 * #nextDouble method will return the result of this transformation and the
 * other `next...` methods are affected accordingly. Any instance of this
 * class contains a  @ref RandomStream called its *base stream*, used to
 * generate its numbers and to which the transformation is applied. Any call
 * to one of the `next...` methods of this class will modify the state of the
 * base stream.
 *
 * The baker transformation is often applied when the  @ref RandomStream is
 * actually an iterator over a point set used for quasi-Monte Carlo
 * integration (see the `hups` package).
 *
 * <div class="SSJ-bigskip"></div>
 */
public class BakerTransformedStream implements RandomStream {

   // The base stream.
   private RandomStream st;

   /**
    * Constructs a new baker transformed stream, using the random numbers
    * from the base stream `stream`.
    */
   public BakerTransformedStream (RandomStream stream) {
      st = stream;
   }

   public void resetStartStream() {
      st.resetStartStream();
   }

   public void resetStartSubstream() {
      st.resetStartSubstream();
   }

   public void resetNextSubstream() {
      st.resetNextSubstream();
   }

/**
 * Returns a string starting with `"Baker transformation of "` and finishing
 * with the result of the call to the `toString` method of the generator.
 */
public String toString() {
      return "Baker transformation of " + st.toString();
   }

   /**
    * Returns the baker transformation of `s.nextDouble()` where `s` is
    * the base stream.
    */
   public double nextDouble() {
      double u = st.nextDouble();
      if (u > 0.5) return 2.0 * (1.0 - u);
      else return u + u;
   }

   /**
    * Generates a random integer in @f$\{i,...,j\}@f$ via  #nextDouble (in
    * which the baker transformation is applied).
    */
   public int nextInt (int i, int j) {
      return i + (int)(nextDouble() * (j - i + 1.0));
   }

   /**
    * Calls `nextArrayOfDouble (u, start, n)` for the base stream, then
    * applies the baker transformation.
    *  @param u            the array in which the numbers will be stored
    *  @param start        the first index of `u` to be used
    *  @param n            the number of random numbers to put in `u`
    */
   public void nextArrayOfDouble (double[] u, int start, int n) {
       st.nextArrayOfDouble (u, start, n);
       for (int i = start; i < start + n; i++)
          if (u[i] > 0.5) u[i] = 2.0 * (1.0 - u[i]);
          else u[i] += u[i];
       }

   /**
    * Fills up the array by calling `nextInt (i, j)`.
    *  @param i            the smallest possible integer to put in `u`
    *  @param j            the largest possible integer to put in `u`
    *  @param u            the array in which the numbers will be stored
    *  @param start        the first index of `u` to be used
    *  @param n            the number of random numbers to put in `u`
    */
   public void nextArrayOfInt (int i, int j, int[] u, int start, int n) {
      for(int ii = start; ii < start + n; ii++)
         u[ii] = nextInt(i,j);
   }

}