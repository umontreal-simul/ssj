/*
 * Class:        AntitheticStream
 * Description:  container class allows the user to force any RandomStream
                 to return antithetic variates
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

/**
 * This container class allows the user to force any  @ref RandomStream to
 * return antithetic variates. That is,  #nextDouble returns @f$1-u@f$
 * instead of @f$u@f$ and the corresponding change is made in  #nextInt. Any
 * instance of this class behaves exactly like a  @ref RandomStream, except
 * that it depends on another random number generator stream, called the
 * *base stream*, to generate its numbers. Any call to one of the `next...`
 * methods of this class will modify the state of the base stream.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class AntitheticStream implements RandomStream {

   // The base stream.
   private RandomStream st;

   /**
    * Constructs a new antithetic stream, using the random numbers from
    * the base stream `stream`.
    */
   public AntitheticStream (RandomStream stream) {
      this.st = stream;
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
	 * Returns a string starting with `"Antithetic of "` and finishing with the result of the call
	 * to the `toString` method of the generator.
	 */
	public String toString() {
      return "Antithetic of " + st.toString();
   }

   /**
    * Returns `1.0 - s.nextDouble()` where `s` is the base stream.
    */
   public double nextDouble() {
      return 1.0 - st.nextDouble();
   }

   /**
    * Returns `j - i - s.nextInt(i, j)` where `s` is the base stream.
    */
   public int nextInt (int i, int j) {
      // pas (j - st.nextInt(0,j-i)), au cas ou le resultat varie.
      return j - i - st.nextInt(i, j);
   }

   /**
    * Calls `nextArrayOfDouble (u, start, n)` for the base stream, then
    * replaces each `u[i]` by `1.0 - u[i]`.
    *  @param u            the array in which the numbers will be stored
    *  @param start        the first index of `u` to be used
    *  @param n            the number of random numbers to put in `u`
    */
   public void nextArrayOfDouble (double[] u, int start, int n) {
       st.nextArrayOfDouble (u, start, n);
       for (int ii = start; ii < start + n; ii++)
          u[ii] = 1.0 - u[ii];
   }

   /**
    * Calls `nextArrayOfInt (i, j, u, start, n)` for the base stream, then
    * replaces each `u[i]` by `j - i - u[i]`.
    *  @param i            the smallest possible integer to put in `u`
    *  @param j            the largest possible integer to put in `u`
    *  @param u            the array in which the numbers will be stored
    *  @param start        the first index of `u` to be used
    *  @param n            the number of random numbers to put in `u`
    */
   public void nextArrayOfInt (int i, int j, int[] u, int start, int n) {
       st.nextArrayOfInt (i, j, u, start, n);
       for (int ii = start; ii < start + n; ii++)
          u[ii] = j - i - u[ii];
   }

}