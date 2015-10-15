/*
 * Class:        UnuranDiscreteInt
 * Description:  create a discrete generator using UNURAN
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