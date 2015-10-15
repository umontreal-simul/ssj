/*
 * Class:        UnuranContinuous
 * Description:  create continuous univariate generators using UNURAN
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
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.rng.RandomStream;

/**
 * This class permits one to create continuous univariate generators using
 * UNURAN via its string API.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_unuran
 */
public class UnuranContinuous extends RandomVariateGen {

   private RandUnuran unuran = new RandUnuran();

   /**
    * Same as  {@link #UnuranContinuous() UnuranContinuous(s, s, genStr)}.
    */
   public UnuranContinuous (RandomStream s, String genStr) {
      if (s == null)
         throw new IllegalArgumentException ("s must not be null.");

      unuran.mainStream = unuran.auxStream = s;
      unuran.init (genStr);
      if (!unuran.isContinuous()) {
         unuran.close();
         throw new IllegalArgumentException ("not a continuous distribution");
      }
   }

   /**
    * Constructs a new continuous random number generator using the UNURAN
    * generator specification string `genStr`, main stream `s`, and
    * auxiliary stream `aux`.
    */
   public UnuranContinuous (RandomStream s, RandomStream aux,
                            String genStr) {
      if (s == null)
         throw new IllegalArgumentException ("s must not be null.");
      if (aux == null)
         throw new IllegalArgumentException ("aux must not be null.");

      unuran.mainStream = s;
      unuran.auxStream = aux;
      unuran.init (genStr);
      if (!unuran.isContinuous()) {
         unuran.close();
         throw new IllegalArgumentException ("not a continuous distribution");
      }
   }


   public double nextDouble() {
      if (unuran.nativeParams == 0)
         throw new IllegalStateException();
      return unuran.getRandCont (unuran.mainStream.nextDouble(), unuran.nativeParams);
   }

   public void nextArrayOfDouble (double[] v, int start, int n) {
      if (v == null || start < 0 || n < 0 || (start+n) > v.length)
         throw new IllegalArgumentException();
      if (unuran.unifArray == null || unuran.unifArray.length < n)
         unuran.unifArray = new double[n];
      if (unuran.mainStream != unuran.auxStream &&
         (unuran.unifAuxArray == null || unuran.unifAuxArray.length < n))
         unuran.unifAuxArray = new double[n];
      unuran.getRandContArray (unuran.nativeParams, unuran.unifArray,
                        unuran.unifAuxArray, v, start, n);
   }

   protected void finalize() {
      unuran.close();
   }

   public Distribution getDistribution() { return null; }

   public RandomStream getStream() { return unuran.mainStream; }

/**
 * Returns the auxiliary random number stream.
 */
public RandomStream getAuxStream() {
      return unuran.auxStream;
   }
}