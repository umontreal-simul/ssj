/*
 * Class:        RandomVariateGenWithCache
 * Description:  random variate generator whose values are cached for efficiency
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
package umontreal.ssj.randvar;

import umontreal.ssj.randvar.RandomVariateGen;
import cern.colt.list.DoubleArrayList;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.probdist.Distribution;

/**
 * This class represents a random variate generator whose values are cached
 * for more efficiency when using common random numbers. An object from this
 * class is constructed with a reference to a  @ref RandomVariateGen instance
 * used to get the random numbers. These numbers are stored in an internal
 * array to be retrieved later. The dimension of the array increases as the
 * values are generated. If the  #nextDouble method is called after the
 * object is reset (by calling  #setCachedValues(DoubleArrayList) ), it gives
 * back the cached values instead of computing new ones. If the cache is
 * exhausted before the generator is reset, new values are computed and added
 * to the cache.
 *
 * Such caching allows for a better performance with common random numbers,
 * when generating random variates is time-consuming. However, using such
 * caching may lead to memory problems if a large quantity of random numbers
 * are needed.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_general
 */
public class RandomVariateGenWithCache extends RandomVariateGen {
   private RandomVariateGen rvg;
   private DoubleArrayList values;
   private int index = 0;
   private boolean caching = true;

   /**
    * Constructs a new cached random variate generator with internal
    * generator `rvg`.
    *  @param rvg          the random variate generator whose values are
    *                      cached.
    *  @exception NullPointerException if `rvg` is `null`.
    */
   public RandomVariateGenWithCache (RandomVariateGen rvg) {
      if (rvg == null)
         throw new NullPointerException
            ("The given random variate generator cannot be null");
      this.rvg = rvg;
      values = new DoubleArrayList();
   }

   /**
    * Constructs a new cached random variate generator with internal
    * generator `rvg`. The `initialCapacity` parameter is used to set the
    * initial capacity of the internal array which can grow as needed; it
    * does not limit the maximal number of cached values.
    *  @param rvg          the random variate generator whose values are
    *                      cached.
    *  @param initialCapacity the number of cached values.
    *  @exception NullPointerException if `rvg` is `null`.
    */
   public RandomVariateGenWithCache (RandomVariateGen rvg,
                                     int initialCapacity) {
      if (rvg == null)
         throw new NullPointerException
            ("The given random variate generator cannot be null");
      this.rvg = rvg;
      values = new DoubleArrayList (initialCapacity);
   }

   /**
    * Determines if the random variate generator is caching values,
    * default being `true`. When caching is turned OFF, the  #nextDouble
    * method simply calls the corresponding method on the internal random
    * variate generator, without storing the generated values.
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
    * Returns a reference to the random variate generator whose values are
    * cached.
    *  @return a reference to the random variate generator whose values
    * are cached.
    */
   public RandomVariateGen getCachedGen() {
      return rvg;
   }

   /**
    * Sets the random variate generator whose values are cached to `rvg`.
    * If the generator is changed, the  #clearCache method is called.
    *  @param rvg          the new random variate generator whose values
    *                      are cached.
    *  @exception NullPointerException if `rvg` is `null`.
    */
   public void setCachedGen (RandomVariateGen rvg) {
      if (rvg == null)
         throw new NullPointerException
            ("The given random variate generator cannot be null");
      if (rvg == this.rvg)
         return;
      this.rvg = rvg;
      clearCache();
   }

   /**
    * Clears the cached values for this cached generator. Any subsequent
    * call will then obtain new values from the internal generator.
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
    * Resets this generator to recover values from the cache. Subsequent
    * calls to  #nextDouble will return the cached random values until all
    * the values are returned. When the array of cached values is
    * exhausted, the internal random variate generator is used to generate
    * new values which are added to the internal array as well. This
    * method is equivalent to calling  #setCacheIndex(int).
    */
   public void initCache() {
      index = 0;
   }

   /**
    * Returns the total number of values cached by this generator.
    *  @return the total number of cached values.
    */
   public int getNumCachedValues() {
      return values.size();
   }

   /**
    * Return the index of the next cached value that will be returned by
    * the generator. If the cache is exhausted, the returned value
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
    * variate generator.
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


   public double nextDouble() {
      if (!caching)
         return rvg.nextDouble();
      else if (index >= values.size()) {
         double v = rvg.nextDouble();
         values.add (v);
         ++index;
         return v;
      }
      else
         return values.getQuick (index++);
   }

   public void nextArrayOfDouble (double[] v, int start, int n) {
      if (!caching) {
         rvg.nextArrayOfDouble (v, start, n);
         return;
      }
      int remainingValues = values.size() - index;
      if (remainingValues < 0)
         remainingValues = 0;
      int ncpy = Math.min (n, remainingValues);
      if (ncpy > 0) {
         System.arraycopy (values.elements(), index, v, start, ncpy);
         index += ncpy;
      }
      int ngen = n - ncpy;
      if (ngen > 0) {
         rvg.nextArrayOfDouble (v, start + ncpy, ngen);
         for (int i = ncpy; i < n; i++) {
            values.add (v[start + i]);
            ++index;
         }
      }
   }

   public RandomStream getStream() {
      return rvg.getStream();
   }

   public Distribution getDistribution() {
      return rvg.getDistribution();
   }
}