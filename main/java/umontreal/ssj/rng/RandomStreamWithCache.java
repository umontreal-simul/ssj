/*
 * Class:        RandomStreamWithCache
 * Description:  random stream whose uniforms are cached for more efficiency
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

import cern.colt.list.DoubleArrayList;

/**
 * This class represents a random stream whose uniforms are cached for more
 * efficiency when using common random numbers. An object from this class is
 * constructed with a reference to a  @ref RandomStream instance used to get
 * the random numbers. These numbers are stored in an internal array to be
 * retrieved later. The dimension of the array increases as the values are
 * generated. If the  #nextDouble method is called after the object is reset,
 * it gives back the cached values instead of computing new ones. If the
 * cache is exhausted before the stream is reset, new values are computed,
 * and added to the cache.
 *
 * Such caching allows for a better performance with common random numbers,
 * when generating uniforms is time-consuming. It can also help with
 * restoring the simulation to a certain state without setting
 * stream-specific seeds. However, using such caching may lead to memory
 * problems if a large quantity of random numbers are needed.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class RandomStreamWithCache implements RandomStream {
   private RandomStream stream;
   private DoubleArrayList values;
   private int index = 0;
   private boolean caching = true;

   /**
    * Constructs a new cached random stream with internal stream `stream`.
    *  @param stream       the random stream whose uniforms are cached.
    *  @exception NullPointerException if `stream` is `null`.
    */
   public RandomStreamWithCache (RandomStream stream) {
      if (stream == null)
         throw new NullPointerException
            ("The given random stream cannot be null");
      this.stream = stream;
      values = new DoubleArrayList();
   }

   /**
    * Constructs a new cached random stream with internal stream `stream`.
    * The `initialCapacity` parameter is used to set the initial capacity
    * of the internal array which can grow as needed; it does not limit
    * the total size of the cache.
    *  @param stream       the random stream whose values are cached.
    *  @param initialCapacity the initial capacity of the cache.
    *  @exception NullPointerException if `stream` is `null`.
    */
   public RandomStreamWithCache (RandomStream stream, int initialCapacity) {
      if (stream == null)
         throw new NullPointerException
            ("The given random stream cannot be null");
      this.stream = stream;
      values = new DoubleArrayList (initialCapacity);
   }

   /**
    * Determines if the random stream is caching values, default being
    * `true`. When caching is turned OFF, the  #nextDouble method simply
    * calls the corresponding method on the internal random stream,
    * without storing the generated uniforms.
    *  @return the caching indicator.
    */
   public boolean isCaching() {
      return caching;
   }

   /**
    * Sets the caching indicator to `caching`. If caching is turned OFF,
    * this method calls  #clearCache to clear the cached values.
    *  @param caching      the new value of the caching indicator.
    */
   public void setCaching (boolean caching) {
      if (this.caching && !caching)
         clearCache();
      this.caching = caching;
   }

   /**
    * Returns a reference to the random stream whose values are cached.
    *  @return a reference to the random stream whose values are cached.
    */
   public RandomStream getCachedStream() {
      return stream;
   }

   /**
    * Sets the random stream whose values are cached to `stream`. If the
    * stream is changed, the  #clearCache method is called to clear the
    * cache.
    *  @param stream       the new random stream whose values will be
    *                      cached.
    *  @exception NullPointerException if `stream` is `null`.
    */
   public void setCachedStream (RandomStream stream) {
      if (stream == null)
         throw new NullPointerException
            ("The given random stream cannot be null");
      if (stream == this.stream)
         return;
      this.stream = stream;
      clearCache();
   }

   /**
    * Clears the cached values for this random stream. Any subsequent call
    * will then obtain new values from the internal stream.
    */
   public void clearCache() {
      //values.clear();
      // Keep the array previously returned by getCachedValues
      // intact to allow caching values for several
      // replications.
      values = new DoubleArrayList();
      index = 0;
   }

   /**
    * Resets this random stream to recover values from the cache.
    * Subsequent calls to  #nextDouble will return the cached uniforms
    * until all the values are returned. When the array of cached values
    * is exhausted, the internal random stream is used to generate new
    * values which are added to the internal array as well. This method is
    * equivalent to calling  #setCacheIndex(int).
    */
   public void initCache() {
      index = 0;
   }

   /**
    * Returns the total number of values cached by this random stream.
    *  @return the total number of cached values.
    */
   public int getNumCachedValues() {
      return values.size();
   }

   /**
    * Return the index of the next cached value that will be returned by
    * the stream. If the cache is exhausted, the returned value
    * corresponds to the value returned by  #getNumCachedValues, and a
    * subsequent call to  #nextDouble will generate a new variate rather
    * than reading a previous one from the cache. If caching is disabled,
    * this always returns 0.
    *  @return the index of the next cached value.
    */
   public int getCacheIndex() {
      return index;
   }

   /**
    * Sets the index, in the cache, of the next value returned by
    * #nextDouble. If `newIndex` is 0, this is equivalent to calling
    * #initCache. If `newIndex` is  #getNumCachedValues, subsequent calls
    * to  #nextDouble will add new values to the cache.
    *  @param newIndex     the new index.
    *  @exception IllegalArgumentException if `newIndex` is negative or
    * greater than or equal to the cache size.
    */
   public void setCacheIndex (int newIndex) {
      if (newIndex < 0 || newIndex > values.size())
         throw new IllegalArgumentException
         ("newIndex must not be negative or greater than the cache size");
      index = newIndex;
   }

   /**
    * Returns an array list containing the values cached by this random
    * stream.
    *  @return the array of cached values.
    */
   public DoubleArrayList getCachedValues() {
      return values;
   }

   /**
    * Sets the array list containing the cached values to `values`. This
    * resets the cache index to the size of the given array.
    *  @param values       the array list of cached values.
    *  @exception NullPointerException if `values` is `null`.
    */
   public void setCachedValues (DoubleArrayList values) {
      if (values == null)
         throw new NullPointerException();
      this.values = values;
      index = values.size();
   }


   public void resetStartStream () {
      stream.resetStartStream();
   }

   public void resetStartSubstream () {
      stream.resetStartSubstream();
   }

   public void resetNextSubstream () {
      stream.resetNextSubstream();
   }

   public double nextDouble () {
      if (!caching)
         return stream.nextDouble();
      else if (index >= values.size()) {
         double v = stream.nextDouble();
         values.add (v);
         ++index;
         return v;
      }
      else
         return values.getQuick (index++);
   }

   public void nextArrayOfDouble (double[] u, int start, int n) {
      if (!caching) {
         stream.nextArrayOfDouble (u, start, n);
         return;
      }
      int remainingValues = values.size() - index;
      if (remainingValues < 0)
         remainingValues = 0;
      int ncpy = Math.min (n, remainingValues);
      if (ncpy > 0) {
         System.arraycopy (values.elements(), index, u, start, ncpy);
         index += ncpy;
      }
      int ngen = n - ncpy;
      if (ngen > 0) {
         stream.nextArrayOfDouble (u, start + ncpy, ngen);
         for (int i = ncpy; i < n; i++) {
            values.add (u[start + i]);
            ++index;
         }
      }
   }

   public int nextInt (int i, int j) {
      return i + (int) (nextDouble () * (j - i + 1));
   }

   public void nextArrayOfInt (int i, int j, int[] u, int start, int n) {
      for (int x = start; x < start + n; x++)
         u[x] = nextInt (i, j);
   }
}