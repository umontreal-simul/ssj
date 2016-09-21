/*
 * Class:        RandomPermutation
 * Description:  Provides methods to randomly shuffle arrays or lists
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

import java.util.Collections;
import java.util.List;
import java.util.ListIterator;
import java.util.RandomAccess;

/**
 * Provides methods to randomly shuffle arrays or lists using a random
 * stream.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class RandomPermutation {
   private static final int SHUFFLE_THRESHOLD = 5;

/**
 * Initializes `array` with the first @f$n@f$ positive integers in natural
 * order as <tt>array</tt>@f$[i-1] = i@f$, for @f$i=1,...,n@f$. The size of
 * `array` must be at least @f$n@f$.
 *  @param array        the array to initialize.
 *  @param n            number of elements initialized.
 */
@Deprecated
   public static void init (byte[] array, int n) {
      for (byte k = 1; k <= n; k++)
         array[k-1] = k;
   }

   /**
    * Similar to  {@link #init() init(byte[], int)}.
    *  @param array        the array to initialize.
    *  @param n            number of elements initialized.
    */
   @Deprecated
   public static void init (short[] array, int n) {
      for (short k = 1; k <= n; k++)
         array[k-1] = k;
   }

   /**
    * Similar to  {@link #init() init(byte[], int)}.
    *  @param array        the array to initialize.
    *  @param n            number of elements initialized.
    */
   @Deprecated
   public static void init (int[] array, int n) {
      for (int k = 1; k <= n; k++)
         array[k-1] = k;
   }

   /**
    * Similar to  {@link #init() init(byte[], int)}.
    *  @param array        the array to initialize.
    *  @param n            number of elements initialized.
    */
   @Deprecated
   public static void init (long[] array, int n) {
      for (int k = 1; k <= n; k++)
         array[k-1] = k;
   }

   /**
    * Similar to  {@link #init() init(byte[], int)}.
    *  @param array        the array to initialize.
    *  @param n            number of elements initialized.
    */
   @Deprecated
   public static void init (float[] array, int n) {
      for (int k = 1; k <= n; k++)
         array[k-1] = k;
   }

   /**
    * Similar to  {@link #init() init(byte[], int)}. <div
    * class="SSJ-bigskip"></div>
    *  @param array        the array to initialize.
    *  @param n            number of elements initialized.
    */
   @Deprecated
   public static void init (double[] array, int n) {
      for (int k = 1; k <= n; k++)
         array[k-1] = k;
   }
/**
 * Same as `java.util.Collections.shuffle(List<?>, Random)`, but uses a
 * @ref RandomStream instead of `java.util.Random`.
 *  @param list         the list being shuffled.
 *  @param stream       the random stream used to generate integers.
 */
@SuppressWarnings("unchecked")
public static void shuffle (List<?> list, RandomStream stream) {
      // The implementation is inspired from Sun's Collections.shuffle
      final int size = list.size ();
      if (size < SHUFFLE_THRESHOLD || list instanceof RandomAccess) {
         for (int i = size; i > 1; i--)
            Collections.swap (list, i - 1, stream.nextInt (0, i - 1));

      } else {
         final Object arr[] = list.toArray ();

         // Shuffle array<
         shuffle (arr, stream);

         // Dump array back into list
         final ListIterator it = list.listIterator ();
         for (Object element : arr) {
            it.next ();
            it.set (element);
         }
      }
   }

   /**
    * Randomly permutes `array` using `stream`. This method permutes the
    * whole array.
    *  @param array        the array being shuffled.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (Object[] array, RandomStream stream) {
      final int size = array.length;
      for (int i = size - 1; i > 0; i--) {
         final int j = stream.nextInt (0, i);
         final Object tmp = array[i];
         array[i] = array[j];
         array[j] = tmp;
      }
   }

   /**
    * Randomly permutes `array` using `stream`. This method permutes the
    * whole array.
    *  @param array        the array being shuffled.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (byte[] array, RandomStream stream) {
      final int size = array.length;
      for (int i = size - 1; i > 0; i--) {
         final int j = stream.nextInt (0, i);
         final byte tmp = array[i];
         array[i] = array[j];
         array[j] = tmp;
      }
   }

   /**
    * Similar to  {@link #shuffle() shuffle(byte[], RandomStream)}.
    *  @param array        the array being shuffled.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (short[] array, RandomStream stream) {
      final int size = array.length;
      for (int i = size - 1; i > 0; i--) {
         final int j = stream.nextInt (0, i);
         final short tmp = array[i];
         array[i] = array[j];
         array[j] = tmp;
      }
   }

   /**
    * Similar to  {@link #shuffle() shuffle(byte[], RandomStream)}.
    *  @param array        the array being shuffled.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (int[] array, RandomStream stream) {
      final int size = array.length;
      for (int i = size - 1; i > 0; i--) {
         final int j = stream.nextInt (0, i);
         final int tmp = array[i];
         array[i] = array[j];
         array[j] = tmp;
      }
   }

   /**
    * Similar to  {@link #shuffle() shuffle(byte[], RandomStream)}.
    *  @param array        the array being shuffled.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (long[] array, RandomStream stream) {
      final int size = array.length;
      for (int i = size - 1; i > 0; i--) {
         final int j = stream.nextInt (0, i);
         final long tmp = array[i];
         array[i] = array[j];
         array[j] = tmp;
      }
   }

   /**
    * Similar to  {@link #shuffle() shuffle(byte[], RandomStream)}.
    *  @param array        the array being shuffled.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (char[] array, RandomStream stream) {
      final int size = array.length;
      for (int i = size - 1; i > 0; i--) {
         final int j = stream.nextInt (0, i);
         final char tmp = array[i];
         array[i] = array[j];
         array[j] = tmp;
      }
   }

   /**
    * Similar to  {@link #shuffle() shuffle(byte[], RandomStream)}.
    *  @param array        the array being shuffled.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (boolean[] array, RandomStream stream) {
      final int size = array.length;
      for (int i = size - 1; i > 0; i--) {
         final int j = stream.nextInt (0, i);
         final boolean tmp = array[i];
         array[i] = array[j];
         array[j] = tmp;
      }
   }

   /**
    * Similar to  {@link #shuffle() shuffle(byte[], RandomStream)}.
    *  @param array        the array being shuffled.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (float[] array, RandomStream stream) {
      final int size = array.length;
      for (int i = size - 1; i > 0; i--) {
         final int j = stream.nextInt (0, i);
         final float tmp = array[i];
         array[i] = array[j];
         array[j] = tmp;
      }
   }

   /**
    * Similar to  {@link #shuffle() shuffle(byte[], RandomStream)}. <div
    * class="SSJ-bigskip"></div>
    *  @param array        the array being shuffled.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (double[] array, RandomStream stream) {
      final int size = array.length;
      for (int i = size - 1; i > 0; i--) {
         final int j = stream.nextInt (0, i);
         final double tmp = array[i];
         array[i] = array[j];
         array[j] = tmp;
      }
   }
/**
 * Partially permutes `list` as follows using `stream`: draws the first
 * @f$k@f$ new elements of `list` randomly among the @f$n@f$ old elements of
 * `list`, assuming that @f$k \le n = @f$ `list.size()`. In other words,
 * @f$k@f$ elements are selected at random without replacement from the
 * @f$n@f$ `list` entries and are placed in the first @f$k@f$ positions, in
 * random order.
 *  @param list         the list being shuffled.
 *  @param k            number of elements selected.
 *  @param stream       the random stream used to generate integers.
 */
@SuppressWarnings("unchecked")
public static void shuffle (List<?> list, int k, RandomStream stream) {
      // @precondition 0 <= k <= n <= size.

      // The implementation is inspired from Sun's Collections.shuffle
      final int size = list.size ();
      if (k < 0 || k > size)
         throw new IllegalArgumentException("k must be   0 <= k <= list.size()");
      if (0 == k) return;
      if (size < SHUFFLE_THRESHOLD || list instanceof RandomAccess) {
         for (int i = 0; i < k; i++) {
            // Get random j in {i,...,n-1} and interchange a[i] with a[j].
            int j = stream.nextInt (i, size-1);
            Collections.swap (list, i, j);
         }

      } else {
         final Object arr[] = list.toArray ();

         // Shuffle array<
         shuffle (arr, size, k, stream);

         // Dump array back into list
         final ListIterator it = list.listIterator ();
         for (Object element : arr) {
            it.next ();
            it.set (element);
         }
      }
   }

   /**
    * Partially permutes `array` as follows using `stream`: draws the new
    * @f$k@f$ elements, `array[0]` to `array[k-1]`, randomly among the old
    * @f$n@f$ elements, `array[0]` to `array[n-1]`, assuming that @f$k
    * \le n \le@f$ `array.length`. In other words, @f$k@f$ elements are
    * selected at random without replacement from the first @f$n@f$ array
    * elements and are placed in the first @f$k@f$ positions, in random
    * order.
    *  @param array        the array being shuffled.
    *  @param n            selection amongst the first n elements.
    *  @param k            number of elements selected.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (Object[] array, int n, int k,
                               RandomStream stream) {
      // @precondition 0 <= k <= n <= a.length.
      // Replace by 
      // if (k < 0 || k > n) throw new IllegalArgumentException();
      // or at least assert 0 <= k && k <= n;
      if (k < 0 || k > n)
         throw new IllegalArgumentException("k must be   0 <= k <= n");
      for (int i = 0; i < k; i++) {
         // Get random j in {i,...,n-1} and interchange a[i] with a[j].
         int j = stream.nextInt (i, n-1);
         Object temp = array[j];
         array[j] = array[i];
         array[i] = temp;
      }
   }

   /**
    * Similar to  {@link #shuffle() shuffle(Object[], n, k,
    * RandomStream)}.
    *  @param array        the array being shuffled.
    *  @param n            selection amongst the first n elements.
    *  @param k            number of elements selected.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (byte[] array, int n, int k,
                               RandomStream stream) {
      // @precondition 0 <= k <= n <= a.length.
      if (k < 0 || k > n)
         throw new IllegalArgumentException("k must be   0 <= k <= n");
      for (int i = 0; i < k; i++) {
         // Get random j in {i,...,n-1} and interchange a[i] with a[j].
         int j = stream.nextInt (i, n-1);
         byte temp = array[j];
         array[j] = array[i];
         array[i] = temp;
      }
   }

   /**
    * Similar to  {@link #shuffle() shuffle(Object[], n, k,
    * RandomStream)}.
    *  @param array        the array being shuffled.
    *  @param n            selection amongst the first n elements.
    *  @param k            number of elements selected.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (short[] array, int n, int k,
                               RandomStream stream) {
      // @precondition 0 <= k <= n <= a.length.
      if (k < 0 || k > n)
         throw new IllegalArgumentException("k must be   0 <= k <= n");
      for (int i = 0; i < k; i++) {
         // Get random j in {i,...,n-1} and interchange a[i] with a[j].
         int j = stream.nextInt (i, n-1);
         short temp = array[j];
         array[j] = array[i];
         array[i] = temp;
      }
   }

   /**
    * Similar to  {@link #shuffle() shuffle(Object[], n, k,
    * RandomStream)}.
    *  @param array        the array being shuffled.
    *  @param n            selection amongst the first n elements.
    *  @param k            number of elements selected.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (int[] array, int n, int k,
                               RandomStream stream) {
      // @precondition 0 <= k <= n <= a.length.
      if (k < 0 || k > n)
         throw new IllegalArgumentException("k must be   0 <= k <= n");
      for (int i = 0; i < k; i++) {
         // Get random j in {i,...,n-1} and interchange a[i] with a[j].
         int j = stream.nextInt (i, n-1);
         int temp = array[j];
         array[j] = array[i];
         array[i] = temp;
      }
   }

   /**
    * Similar to  {@link #shuffle() shuffle(Object[], n, k,
    * RandomStream)}.
    *  @param array        the array being shuffled.
    *  @param n            selection amongst the first n elements.
    *  @param k            number of elements selected.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (long[] array, int n, int k,
                               RandomStream stream) {
      // @precondition 0 <= k <= n <= a.length.
      if (k < 0 || k > n)
         throw new IllegalArgumentException("k must be   0 <= k <= n");
      for (int i = 0; i < k; i++) {
         // Get random j in {i,...,n-1} and interchange a[i] with a[j].
         int j = stream.nextInt (i, n-1);
         long temp = array[j];
         array[j] = array[i];
         array[i] = temp;
      }
   }

   /**
    * Similar to  {@link #shuffle() shuffle(Object[], n, k,
    * RandomStream)}.
    *  @param array        the array being shuffled.
    *  @param n            selection amongst the first n elements.
    *  @param k            number of elements selected.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (char[] array, int n, int k,
                               RandomStream stream) {
      // @precondition 0 <= k <= n <= a.length.
      if (k < 0 || k > n)
         throw new IllegalArgumentException("k must be   0 <= k <= n");
      for (int i = 0; i < k; i++) {
         // Get random j in {i,...,n-1} and interchange a[i] with a[j].
         int j = stream.nextInt (i, n-1);
         char temp = array[j];
         array[j] = array[i];
         array[i] = temp;
      }
   }

   /**
    * Similar to  {@link #shuffle() shuffle(Object[], n, k,
    * RandomStream)}.
    *  @param array        the array being shuffled.
    *  @param n            selection amongst the first n elements.
    *  @param k            number of elements selected.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (boolean[] array, int n, int k,
                               RandomStream stream) {
      // @precondition 0 <= k <= n <= a.length.
      if (k < 0 || k > n)
         throw new IllegalArgumentException("k must be   0 <= k <= n");
      for (int i = 0; i < k; i++) {
         // Get random j in {i,...,n-1} and interchange a[i] with a[j].
         int j = stream.nextInt (i, n-1);
         boolean temp = array[j];
         array[j] = array[i];
         array[i] = temp;
      }
   }

   /**
    * Similar to  {@link #shuffle() shuffle(Object[], n, k,
    * RandomStream)}.
    *  @param array        the array being shuffled.
    *  @param n            selection amongst the first n elements.
    *  @param k            number of elements selected.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (float[] array, int n, int k,
                               RandomStream stream) {
      // @precondition 0 <= k <= n <= a.length.
      if (k < 0 || k > n)
         throw new IllegalArgumentException("k must be   0 <= k <= n");
      for (int i = 0; i < k; i++) {
         // Get random j in {i,...,n-1} and interchange a[i] with a[j].
         int j = stream.nextInt (i, n-1);
         float temp = array[j];
         array[j] = array[i];
         array[i] = temp;
      }
   }

   /**
    * Similar to  {@link #shuffle() shuffle(Object[], n, k,
    * RandomStream)}.
    *  @param array        the array being shuffled.
    *  @param n            selection amongst the first n elements.
    *  @param k            number of elements selected.
    *  @param stream       the random stream used to generate random
    *                      numbers.
    */
   public static void shuffle (double[] array, int n, int k,
                               RandomStream stream) {
      // @precondition 0 <= k <= n <= a.length.
      if (k < 0 || k > n)
         throw new IllegalArgumentException("k must be   0 <= k <= n");
      for (int i = 0; i < k; i++) {
         // Get random j in {i,...,n-1} and interchange a[i] with a[j].
         int j = stream.nextInt (i, n-1);
         double temp = array[j];
         array[j] = array[i];
         array[i] = temp;
      }
   }

}