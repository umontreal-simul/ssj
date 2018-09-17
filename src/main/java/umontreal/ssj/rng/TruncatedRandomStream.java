/*
 * Class:        TruncatedRandomStream
 * Description:  container random stream generating numbers in an interval
                 (a,b) instead of in (0,1)
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

import umontreal.ssj.rng.RandomStream;

/**
 * Represents a container random stream generating numbers in an interval
 * @f$(a,b)@f$ instead of in @f$(0,1)@f$, where @f$0\le a < b \le1@f$, by
 * using the contained stream. If `nextDouble` returns @f$u@f$ for the
 * contained stream, it will return @f$v = a + (b-a)u@f$, which is uniform
 * over @f$(a,b)@f$, for the truncated stream. The method `nextInt` returns
 * the integer that corresponds to @f$v@f$ (by inversion); this integer is no
 * longer uniformly distributed in general.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class TruncatedRandomStream implements RandomStream {
   private RandomStream stream;
   private double a;
   private double bminusa;
   public TruncatedRandomStream (RandomStream stream, double a, double b) {
      if (stream == null)
         throw new NullPointerException ("The given stream must not be null");
      if (a >= b)
         throw new IllegalArgumentException ("a must be smaller than b");
      if (a < 0 || b < 0 || a > 1 || b > 1)
         throw new IllegalArgumentException ("a and b must be in [0, 1]");
      this.stream = stream;
      this.a = a;
      bminusa = b - a;
   }

   public void resetStartStream () {
      stream.resetStartStream ();
   }

   public void resetStartSubstream () {
      stream.resetStartSubstream ();
   }

   public void resetNextSubstream () {
      stream.resetNextSubstream ();
   }

   public double nextDouble () {
      double v = stream.nextDouble ();
      return a + v * bminusa;
   }

   public void nextArrayOfDouble (double[] u, int start, int n) {
      stream.nextArrayOfDouble (u, start, n);
      for (int i = start; i < start + n; i++)
         u[i] = a + u[i] * bminusa;
   }

   public int nextInt (int i, int j) {
      return i + (int) (nextDouble () * (j - i + 1));
   }

   public void nextArrayOfInt (int i, int j, int[] u, int start, int n) {
      for (int x = start; x < start + n; x++)
         u[x] = nextInt (i, j);
   }
}