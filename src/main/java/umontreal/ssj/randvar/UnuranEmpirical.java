/*
 * Class:        UnuranEmpirical
 * Description:  create generators for empirical distributions using UNURAN
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
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PiecewiseLinearEmpiricalDist;
import umontreal.ssj.rng.RandomStream;

/**
 * This class permits one to create generators for empirical and
 * quasi-empirical univariate distributions using UNURAN via its string
 * interface. The empirical data can be read from a file, from an array, or
 * simply encoded into the generator specification string. When reading from
 * a file or an array, the generator specification string must *not* contain
 * a distribution specification string.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_unuran
 */
public class UnuranEmpirical extends RandomVariateGen {
   private RandUnuran unuran = new RandUnuran();

   /**
    * Constructs a new empirical univariate generator using the
    * specification string `genStr` and stream `s`.
    */
   public UnuranEmpirical (RandomStream s, String genStr) {
      if (s == null)
         throw new IllegalArgumentException ("s must not be null.");

      unuran.mainStream = unuran.auxStream = s;
      unuran.init (genStr);
      if (!unuran.isEmpirical()) {
         unuran.close();
         throw new IllegalArgumentException ("not an empirical distribution");
      }
   }

   /**
    * Constructs a new empirical univariate generator using the
    * specification string `genStr`, with main stream `s` and auxiliary
    * stream `aux`.
    */
   public UnuranEmpirical (RandomStream s, RandomStream aux, String genStr) {
      if (s == null)
         throw new IllegalArgumentException ("s must not be null.");
      if (aux == null)
         throw new IllegalArgumentException ("aux must not be null.");

      unuran.mainStream = s;
      unuran.auxStream = aux;
      unuran.init (genStr);
      if (!unuran.isEmpirical()) {
         unuran.close();
         throw new IllegalArgumentException ("not an empirical distribution");
      }
   }

   /**
    * Same as  {@link #UnuranEmpirical() UnuranEmpirical(s, s, dist,
    * genStr)}.
    */
   public UnuranEmpirical (RandomStream s,
                           PiecewiseLinearEmpiricalDist dist, String genStr) {
      if (s == null)
         throw new IllegalArgumentException ("s must not be null.");

      unuran.mainStream = unuran.auxStream = s;
      String gstr = readDistr (dist) + 
        (genStr == null || genStr.equals ("") ? "" : "&" + genStr);
      unuran.init (gstr);
      if (!unuran.isEmpirical()) {
         unuran.close();
         throw new IllegalArgumentException ("not an empirical distribution");
      }
   }

   /**
    * Same as  {@link #UnuranEmpirical() UnuranEmpirical(s, aux, genStr)},
    * but reading the observations from the empirical distribution `dist`.
    * The `genStr` argument must not contain a distribution part because
    * the distribution will be generated from the input stream reader.
    */
   public UnuranEmpirical (RandomStream s, RandomStream aux,
                           PiecewiseLinearEmpiricalDist dist, String genStr) {
      if (s == null)
         throw new IllegalArgumentException ("s must not be null.");
      if (aux == null)
         throw new IllegalArgumentException ("aux must not be null.");

      unuran.mainStream = s;
      unuran.auxStream = aux;
      String gstr = readDistr (dist) + 
        (genStr == null || genStr.equals ("") ? "" : "&" + genStr);
      unuran.init (gstr);
      if (!unuran.isEmpirical()) {
         unuran.close();
         throw new IllegalArgumentException ("not an empirical distribution");
      }
   }

   // Constructs and returns a distribution string for empirical distr.
   private String readDistr (PiecewiseLinearEmpiricalDist dist) {
      StringBuffer sb = new StringBuffer ("distr=cemp; data=(");
      boolean first = true;

      for (int i = 0; i < dist.getN(); i++) {
         if (first)
            first = false;
         else
            sb.append (",");
         sb.append (dist.getObs (i));
      }
      sb.append (")");
      return sb.toString();
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