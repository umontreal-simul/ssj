/*
 * Class:        F2wStructure
 * Description:  Tools for point sets and sequences based on field F_{2^w}
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
package umontreal.ssj.hups;
import java.io.*;
import java.util.*;

/**
 * This class implements methods and fields needed by the classes
 * @ref umontreal.ssj.hups.F2wNetLFSR,
 * @ref umontreal.ssj.hups.F2wNetPolyLCG,
 * @ref umontreal.ssj.hups.F2wCycleBasedLFSR and
 * @ref umontreal.ssj.hups.F2wCycleBasedPolyLCG. It also stores the
 * parameters of these point sets which will contain @f$2^{rw}@f$ points (see
 * the meaning of @f$r@f$ and @f$w@f$ below). The parameters can be stored as
 * a polynomial @f$P(z)@f$ over @f$\mathbb F_{2^w}[z]@f$
 * @f[
 *   P(z) = z^r + \sum_{i=1}^r b_i z^{r-i}
 * @f]
 * where @f$b_i\in\mathbb F_{2^w}@f$ for @f$i=1,…,r@f$. Let @f$\zeta@f$ be
 * the root of an irreducible polynomial @f$Q(z)\in\mathbb F_2[z]@f$. It is
 * well known that @f$\zeta@f$ is a generator of the finite field @f$\mathbb
 * F_{2^w}@f$. The elements of @f$\mathbb F_{2^w}@f$ are represented using
 * the polynomial ordered basis @f$(1,\zeta,…,\zeta^{w-1})@f$.
 *
 * In this class, only the non-zero coefficients of @f$P(z)@f$ are stored. It
 * is stored as
 * @f[
 *   P(z) = z^{\mathtt{r}} + \sum_{i=0}^{\mathtt{nbcoeff}} {\mathtt{coeff[}}i{\mathtt{]}} z^{{\mathtt{nocoeff[}}i{\mathtt{]}}}
 * @f]
 * where the coefficients in `coeff[]` represent the non-zero coefficients
 * @f$b_i@f$ of @f$P(z)@f$ using the polynomial basis. The finite field
 * @f$\mathbb F_{2^w}@f$ used is defined by the polynomial
 * @f[
 *   Q(z) = z^w + \sum_{i=1}^w a_i z^{w-i}
 * @f]
 * where @f$a_i\in\mathbb F_2@f$, for @f$i=1,…,w@f$. Polynomial @f$Q@f$ is
 * stored as the bit vector `modQ` = @f$(a_w,…,a_1)@f$.
 *
 * The class also stores the parameter `step` that is used by the classes
 * @ref umontreal.ssj.hups.F2wNetLFSR,
 * @ref umontreal.ssj.hups.F2wNetPolyLCG,
 * @ref umontreal.ssj.hups.F2wCycleBasedLFSR and
 * @ref umontreal.ssj.hups.F2wCycleBasedPolyLCG. This parameter is such that
 * the implementation of the recurrence will output a value at every `step`
 * iterations.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 *
 * ## Additional Information
 *
 * <tt>
 * <div class="SSJ-minipage">
 * <table class="SSJ-table SSJ-has-hlines">
 * <tr>
 *   <td colspan="2" class="c">Directory LFSR_equid_max</td>
 * </tr><tr class="bt">
 *   <td class="c bl br">Filename</td>
 *   <td class="c bl br">Num of poly.</td>
 * </tr><tr class="bt">
 *   <td class="c bl br">f2wR2_W5.dat</td>
 *   <td class="c bl br">2358</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W6.dat</td>
 *   <td class="c bl br">1618</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W7.dat</td>
 *   <td class="c bl br">507</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W8.dat</td>
 *   <td class="c bl br">26</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W9.dat</td>
 *   <td class="c bl br">3</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W4.dat</td>
 *   <td class="c bl br">369</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W5.dat</td>
 *   <td class="c bl br">26</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W6.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR4_W3.dat</td>
 *   <td class="c bl br">117</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR4_W4.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR5_W2.dat</td>
 *   <td class="c bl br">165</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR5_W3.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR6_W2.dat</td>
 *   <td class="c bl br">36</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR6_W3.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR7_W2.dat</td>
 *   <td class="c bl br">10</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR8_W2.dat</td>
 *   <td class="c bl br">11</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR9_W2.dat</td>
 *   <td class="c bl br">1</td>
 * </tr>
 * </table>
 *  </div> <div class="SSJ-minipage">
 * <table class="SSJ-table SSJ-has-hlines">
 * <tr>
 *   <td colspan="2" class="c">Directory LFSR_equid_sum</td>
 * </tr><tr class="bt">
 *   <td class="c bl br">Filename</td>
 *   <td class="c bl br">Num of poly.</td>
 * </tr><tr class="bt">
 *   <td class="c bl br">f2wR2_W5.dat</td>
 *   <td class="c bl br">2276</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W6.dat</td>
 *   <td class="c bl br">1121</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W7.dat</td>
 *   <td class="c bl br">474</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W8.dat</td>
 *   <td class="c bl br">37</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W9.dat</td>
 *   <td class="c bl br">6</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W4.dat</td>
 *   <td class="c bl br">381</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W5.dat</td>
 *   <td class="c bl br">65</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W6.dat</td>
 *   <td class="c bl br">7</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR4_W3.dat</td>
 *   <td class="c bl br">154</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR4_W4.dat</td>
 *   <td class="c bl br">2</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR5_W2.dat</td>
 *   <td class="c bl br">688</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR5_W3.dat</td>
 *   <td class="c bl br">5</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR6_W2.dat</td>
 *   <td class="c bl br">70</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR6_W3.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR7_W2.dat</td>
 *   <td class="c bl br">9</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR8_W2.dat</td>
 *   <td class="c bl br">3</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR9_W2.dat</td>
 *   <td class="c bl br">3</td>
 * </tr>
 * </table>
 *  </div> <div class="SSJ-bigskip"></div>
 *
 * <div class="SSJ-minipage">
 * <table class="SSJ-table SSJ-has-hlines">
 * <tr>
 *   <td colspan="2" class="c">Directory LFSR_mindist_max</td>
 * </tr><tr class="bt">
 *   <td class="c bl br">Filename</td>
 *   <td class="c bl br">Num of poly.</td>
 * </tr><tr class="bt">
 *   <td class="c bl br">f2wR2_W5.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W6.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W7.dat</td>
 *   <td class="c bl br">2</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W8.dat</td>
 *   <td class="c bl br">2</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W9.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W4.dat</td>
 *   <td class="c bl br">2</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W5.dat</td>
 *   <td class="c bl br">2</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W6.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR4_W3.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR4_W4.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR5_W2.dat</td>
 *   <td class="c bl br">2</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR5_W3.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR6_W2.dat</td>
 *   <td class="c bl br">4</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR6_W3.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR7_W2.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR8_W2.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR9_W2.dat</td>
 *   <td class="c bl br">1</td>
 * </tr>
 * </table>
 *  </div> <div class="SSJ-minipage">
 * <table class="SSJ-table SSJ-has-hlines">
 * <tr>
 *   <td colspan="2" class="c">Directory LFSR_mindist_sum</td>
 * </tr><tr class="bt">
 *   <td class="c bl br">Filename</td>
 *   <td class="c bl br">Num of poly.</td>
 * </tr><tr class="bt">
 *   <td class="c bl br">f2wR2_W5.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W6.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W7.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W8.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W9.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W4.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W5.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W6.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR4_W3.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR4_W4.dat</td>
 *   <td class="c bl br">2</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR5_W2.dat</td>
 *   <td class="c bl br">2</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR5_W3.dat</td>
 *   <td class="c bl br">2</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR6_W2.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR6_W3.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR7_W2.dat</td>
 *   <td class="c bl br">2</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR8_W2.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR9_W2.dat</td>
 *   <td class="c bl br">2</td>
 * </tr>
 * </table>
 *  </div>
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-minipage">
 * <table class="SSJ-table SSJ-has-hlines">
 * <tr>
 *   <td colspan="2" class="c">Directory LFSR_tvalue_max</td>
 * </tr><tr class="bt">
 *   <td class="c bl br">Filename</td>
 *   <td class="c bl br">Num of poly.</td>
 * </tr><tr class="bt">
 *   <td class="c bl br">f2wR2_W5.dat</td>
 *   <td class="c bl br">7</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W6.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W7.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W8.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W9.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W4.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W5.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W6.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR4_W3.dat</td>
 *   <td class="c bl br">2</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR4_W4.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR5_W2.dat</td>
 *   <td class="c bl br">14</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR5_W3.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR6_W2.dat</td>
 *   <td class="c bl br">2</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR6_W3.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR7_W2.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR8_W2.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR9_W2.dat</td>
 *   <td class="c bl br">1</td>
 * </tr>
 * </table>
 *  </div> <div class="SSJ-minipage">
 * <table class="SSJ-table SSJ-has-hlines">
 * <tr>
 *   <td colspan="2" class="c">Directory LFSR_tvalue_sum</td>
 * </tr><tr class="bt">
 *   <td class="c bl br">Filename</td>
 *   <td class="c bl br">Num of poly.</td>
 * </tr><tr class="bt">
 *   <td class="c bl br">f2wR2_W5.dat</td>
 *   <td class="c bl br">15</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W6.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W7.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W8.dat</td>
 *   <td class="c bl br">2</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR2_W9.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W4.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W5.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR3_W6.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR4_W3.dat</td>
 *   <td class="c bl br">2</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR4_W4.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR5_W2.dat</td>
 *   <td class="c bl br">13</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR5_W3.dat</td>
 *   <td class="c bl br">2</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR6_W2.dat</td>
 *   <td class="c bl br">12</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR6_W3.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR7_W2.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR8_W2.dat</td>
 *   <td class="c bl br">1</td>
 * </tr><tr>
 *   <td class="c bl br">f2wR9_W2.dat</td>
 *   <td class="c bl br">1</td>
 * </tr>
 * </table>
 *  </div>
 * </tt>
 */
public class F2wStructure {

   private final int ALLONES = 2147483647; // 2^31-1 --> 01111...1
   int w;
   int r;
   int numBits;
   private int modQ;
   private int step;
   private int[] coeff;
   private int[] nocoeff;
   private int nbcoeff;
   int S;
   private int maskw;
   private int maskrw;
   private int maskZrm1;
   private int mask31;
   private int t;
   private int masktrw;
   private int[] maskv;
   int state;
   int output;            // augmented state
   double normFactor;
   double EpsilonHalf;
   final static int MBL = 140; //maximum of bytes in 1 line
   //92 bytes for a number of coeff = 15


   private void init (int w, int r, int modQ, int step,
      int nbcoeff, int coeff[], int nocoeff[])
   {
      normFactor = 1.0 / (1L << 31); // 4.65661287307739258e-10;
      EpsilonHalf = 0.5*normFactor;
      numBits = 31;
      this.step = step;
      this.w = w;
      this.r = r;
      S = 31 - r * w;
      mask31 = ~(1 << 31);
      maskw = (1 << w) - 1;
      maskrw = ((1 << (r * w)) - 1) << S;
      maskZrm1 = (ALLONES >> (r * w)) ^ (ALLONES >> ((r - 1) * w));
      this.modQ = modQ;
      this.nbcoeff = nbcoeff;
      this.nocoeff = new int[nbcoeff];
      this.coeff = new int[nbcoeff];
      for (int j = 0; j < nbcoeff; j++) {
         this.nocoeff[j] = nocoeff[j];
         this.coeff[j] = coeff[j];
      }
   }

   void initParamLFSR ()
   {
      t = (31 - r * w) / w;
      masktrw = (~0) << (31 - (t + r) * w) & mask31;
      maskv = new int[r];
      for (int j = 0; j < r; j++) {
         maskv[j] = maskw << (S + ((r - 1 - j) * w));
      }
   }

   /**
    * Constructs a `F2wStructure` object that contains the parameters of a
    * polynomial in @f$\mathbb F_{2^w}[z]@f$, as well as a stepping
    * parameter.
    */
   F2wStructure (int w, int r, int modQ, int step, int nbcoeff,
                 int coeff[], int nocoeff[])
   {
      init (w, r, modQ, step, nbcoeff, coeff, nocoeff);
   }

   /**
    * Constructs a polynomial in @f$\mathbb F_{2^w}[z]@f$ after reading
    * its parameters from file `filename`; the parameters of this
    * polynomial are stored at line number `no` of `filename`. The files
    * are kept in different directories depending on the criteria used in
    * the searches for the parameters defining the polynomials. The
    * different criteria for the searches and the theory behind it are
    * described in @cite rPAN04d, @cite rPAN04t&thinsp;. The existing
    * files and the number of polynomials they contain are given in the
    * following tables. The first table below contains files in
    * subdirectory <tt>LFSR_equid_max</tt>. The name of each file
    * indicates the value of @f$r@f$ and @f$w@f$ for the polynomials. For
    * example, file <tt>f2wR2_W5.dat</tt> in directory
    * <tt>LFSR_equid_max</tt> contains the parameters of 2358 polynomials
    * with @f$r=2@f$ and @f$w=5@f$. For example, to use the 5<em>-th</em>
    * polynomial of file <tt>f2wR2_W5.dat</tt>, one may call
    * <tt>F2wStructure("f2wR2_W5.dat", 5)</tt>. The files of parameters
    * have been stored at the address
    * [http://simul.iro.umontreal.ca/ssj/dataF2w/Panneton/](http://simul.iro.umontreal.ca/ssj/dataF2w/Panneton/).
    * The files should be copied in the user directory, and passed as
    * parameter to the constructor.
    */
   F2wStructure (String filename, int no)
   {
     // If filename can be found starting from the program's directory,
     // it will be used; otherwise, the filename in the Jar archive will
     // be used.
     BufferedReader input;
     try {
       if ((new File (filename)).exists()) {
          input = new BufferedReader (new FileReader (filename));
       } else {
          // does not work anymore since the files and directories have been removed
          // from package hups and put instead on the WWW page.
          // Should be read with protocol http as in class DigitalNetFromFile
          DataInputStream dataInput;
          dataInput = new DataInputStream (
             F2wStructure.class.getClassLoader().getResourceAsStream (
                 "umontreal/ssj/hups/dataF2w/Panneton/" + filename));
          input = new BufferedReader (new InputStreamReader (dataInput));
       }
       initFromReader (filename, input, no);
       input.close ();

     } catch (Exception e) {
       System.out.println ("IO Error: problems finding file " + filename);
       System.exit (1);
     }
   }


   private int multiplyz (int a, int k)
   {
      int i;
      if (k == 0)
         return a & maskw;
      else {
         for (i = 0; i < k; i++) {
            if ((1 & a) == 1) {
               a = (a >> 1) ^ modQ;
            } else
               a = a >> 1;
         }
         return a & maskw;
      }
   }

   /**
    * This method returns the product @f$rw@f$.
    */
   int getLog2N ()
   {
      return r * w;
   }

   /**
    * Method that multiplies two elements in @f$\mathbb F_{2^w}@f$.
    */
   int multiply (int a, int b)

   {
      int i;
      int res = 0, verif = 1;
      for (i = 0; i < w; i++) {
         if ((b & verif) == verif)
            res ^= multiplyz (a, w - 1 - i);
         verif <<= 1;
      }
      return res & maskw;
   }


   void initF2wLFSR ()     // Initialisation de l'etat d'un LFSR
   {
      int v, d = 0;
      int tempState;

      tempState = state << S;
      output = tempState;
      for (int i = 1; i <= t; i++) {
         d = 0;
         for (int j = 0; j < nbcoeff; j++) {
            v = (tempState & maskv[nocoeff[j]]) >>
                 (S + (r - 1 - nocoeff[j]) * w);
            d ^= multiply (coeff[j], v);
         }
         output |= d << (S - i * w);
         tempState = (output << (i * w)) & maskrw;
      }
   }


   void F2wLFSR ()       // Une iteration d'un LFSR
   {
      int v, d = 0;
      int tempState;
      for (int i = 0; i < step; i++) {
         tempState = (output << (t * w)) & maskrw;
         d = 0;
         for (int j = 0; j < nbcoeff; j++) {
            v = (tempState & maskv[nocoeff[j]]) >>
                (S + (r - 1 - nocoeff[j]) * w);
            d ^= multiply (coeff[j], v);
         }
         output = ((output << w) & masktrw) |
                  (d << (31 - (r + t) * w));
      }
      state = (output & maskrw) >> S;
   }


   int F2wPolyLCG ()    // Une iteration d'un PolyLCG
   {
      int Zrm1, d;
      for (int i = 0; i < step; i++) {
         Zrm1 = (state & maskZrm1) >> S;
         state = (state >> w) & maskrw;
         for (int j = 0; j < nbcoeff; j++)
            state ^=
               multiply (coeff[j], Zrm1) << (S + (r - 1 - nocoeff[j]) * w);
      }
      return state;
   }

/**
 * Prints the content of file `filename`. See the constructor above for the
 * conditions on `filename`.
 */
public static void print (String filename)
   {
     BufferedReader input;
     try {
       if ((new File (filename)).exists()) {
          input = new BufferedReader (new FileReader (filename));
       } else {
          DataInputStream dataInput;
          dataInput = new DataInputStream (
             F2wStructure.class.getClassLoader().getResourceAsStream (
                 "umontreal/ssj/hups/dataF2w/" + filename));
          input = new BufferedReader (new InputStreamReader (dataInput));
       }

     String s;
     for (int i = 0; i < 4; i++)
        input.readLine ();
     while ((s = input.readLine ()) != null)
        System.out.println (s);
     input.close ();

     } catch (Exception e) {
       System.out.println ("IO Error: problems reading file " + filename);
       System.exit (1);
     }
   }

   /**
    * This method returns a string containing the polynomial @f$P(z)@f$
    * and the stepping parameter.
    */
   public String toString ()
   {
      StringBuffer sb = new StringBuffer ("z^");
      sb.append (r);
      for (int j = nbcoeff - 1; j >= 0; j--)
         sb.append (" + (" + coeff[j] + ") z^" + nocoeff[j]);
      sb.append ("   modQ = " + modQ + "    w = " + w + "   step = " + step);
      return sb.toString ();
   }

    private void initFromReader (String filename, BufferedReader input, int no)
    {
      int w, r, modQ, step, nbcoeff;
      int coeff[], nocoeff[];
      StringTokenizer line;
      int nl = no + 4;

      try {
        for (int j = 1; j < nl ; j++)
          input.readLine ();

        line = new StringTokenizer (input.readLine ());
        w = Integer.parseInt (line.nextToken ());
        r = Integer.parseInt (line.nextToken ());
        modQ = Integer.parseInt (line.nextToken ());
        step = Integer.parseInt (line.nextToken ());
        nbcoeff = Integer.parseInt (line.nextToken ());
        nocoeff = new int[nbcoeff];
        coeff = new int[nbcoeff];
        for (int i = 0; i < nbcoeff; i++) {
          coeff[i] = Integer.parseInt (line.nextToken ());
          nocoeff[i] = Integer.parseInt (line.nextToken ());
        }
        init (w, r, modQ, step, nbcoeff, coeff, nocoeff);
        input.close ();

      } catch (Exception e) {
        System.out.println ("Input Error: problems reading file " + filename);
        System.exit (1);
      }
    }
  }