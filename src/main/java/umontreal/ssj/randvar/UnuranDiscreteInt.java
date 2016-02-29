/*
 * Class:        UnuranDiscreteInt
 * Description:  create a discrete generator using UNURAN
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
import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.rng.RandomStream;

/**
 * This class permits one to create a discrete univariate generator using
 * UNURAN via its string API.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_unuran
 */
public class UnuranDiscreteInt extends RandomVariateGenInt {

   private RandUnuran unuran = new RandUnuran();

   /**
    * Same as  {@link #UnuranDiscreteInt() UnuranDiscreteInt(s, s,
    * genStr)}.
    */
   public UnuranDiscreteInt (RandomStream s, String genStr) {
      if (s == null)
         throw new IllegalArgumentException ("mainStream must not be null.");

      unuran.mainStream = unuran.auxStream = s;
      unuran.init (genStr);
      if (!unuran.isDiscrete()) {
         unuran.close();
         throw new IllegalArgumentException ("not a discrete distribution");
      }
   }

   /**
    * Constructs a new discrete random number generator using the UNURAN
    * generator specification string `genStr`, main stream `s`, and
    * auxiliary stream `aux`.
    */
   public UnuranDiscreteInt (RandomStream s, RandomStream aux,
                             String genStr) {
      if (s == null)
         throw new IllegalArgumentException ("mainStream must not be null.");
      if (aux == null)
         throw new IllegalArgumentException ("auxStream must not be null.");

      unuran.mainStream = s;
      unuran.auxStream = aux;
      unuran.init (genStr);
      if (!unuran.isDiscrete()) {
         unuran.close();
         throw new IllegalArgumentException ("not a discrete distribution");
      }
   }


   public int nextInt() {
      if (unuran.nativeParams == 0)
         throw new IllegalStateException();
      return unuran.getRandDisc (unuran.mainStream.nextDouble(), unuran.nativeParams);
   }

   public void nextArrayOfInt (int[] v, int start, int n) {
      if (v == null || start < 0 || n < 0 || (start+n) > v.length)
         throw new IllegalArgumentException();
      if (unuran.unifArray == null || unuran.unifArray.length < n)
         unuran.unifArray = new double[n];
      if (unuran.mainStream != unuran.auxStream &&
         (unuran.unifAuxArray == null || unuran.unifAuxArray.length < n))
         unuran.unifAuxArray = new double[n];
      unuran.getRandDiscArray (unuran.nativeParams, unuran.unifArray,
                        unuran.unifAuxArray, v, start, n);
   }

   protected void finalize() {
      unuran.close();
   }

   public DiscreteDistributionInt getDistribution() { return null; }

   public RandomStream getStream() { return unuran.mainStream; }

/**
 * Returns the auxiliary random number stream.
 */
public RandomStream getAuxStream() {
      return unuran.auxStream;
   }
}