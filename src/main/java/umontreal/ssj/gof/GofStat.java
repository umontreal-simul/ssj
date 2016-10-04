/*
 * Class:        GofStat
 * Description:  Goodness-of-fit test statistics
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
package umontreal.ssj.gof;
import cern.colt.list.*;
import umontreal.ssj.util.*;
import umontreal.ssj.probdist.*;
import java.util.Arrays;

/**
 * This class provides methods to compute several types of EDF
 * goodness-of-fit test statistics and to apply certain transformations to a
 * set of observations. This includes the probability integral transformation
 * @f$U_i = F(X_i)@f$, as well as the power ratio and iterated spacings
 * transformations @cite tSTE86a&thinsp;. Here, @f$U_{(0)}, …, U_{(n-1)}@f$
 * stand for @f$n@f$ observations @f$U_0,…,U_{n-1}@f$ sorted by increasing
 * order, where @f$0\le U_i\le1@f$.
 *
 * Note: This class uses the Colt library.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class GofStat {
   private GofStat() {}

/**
 * @name Transforming the observations
 * @{
 */

   // Used in discontinuous distributions
   private static double EPSILOND = 1.0E-15;

/**
 * Applies the probability integral transformation @f$U_i = F (V_i)@f$ for
 * @f$i = 0, 1, …, n-1@f$, where @f$F@f$ is a *continuous* distribution
 * function, and returns the result as an array of length @f$n@f$. @f$V@f$
 * represents the @f$n@f$ observations contained in `data`, and @f$U@f$, the
 * returned transformed observations. If `data` contains random variables
 * from the distribution function `dist`, then the result will contain
 * uniform random variables over @f$[0,1]@f$.
 *  @param data         array of observations to be transformed
 *  @param dist         assumed distribution of the observations
 *  @return the array of transformed observations
 */
public static DoubleArrayList unifTransform (DoubleArrayList data,
                                                ContinuousDistribution dist) {
      double[] v = data.elements();
      int n = data.size();

      double[] u = new double[n];
      for (int i = 0; i < n; i++)
         u[i] = dist.cdf (v[i]);
      return new DoubleArrayList(u);
   }

   /**
    * Applies the transformation @f$U_i = F (V_i)@f$ for @f$i = 0, 1, …,
    * n-1@f$, where @f$F@f$ is a *discrete* distribution function, and
    * returns the result as an array of length @f$n@f$. @f$V@f$ represents
    * the @f$n@f$ observations contained in `data`, and @f$U@f$, the
    * returned transformed observations.
    *
    * Note: If @f$V@f$ are the values of random variables with
    * distribution function `dist`, then the result will contain the
    * values of *discrete* random variables distributed over the set of
    * values taken by `dist`, not uniform random variables over
    * @f$[0,1]@f$.
    *  @param data         array of observations to be transformed
    *  @param dist         assumed distribution of the observations
    *  @return the array of transformed observations
    */
   public static DoubleArrayList unifTransform (DoubleArrayList data,
                                                DiscreteDistribution dist) {
       double[] v = data.elements();
       int n = data.size();

       double[] u = new double[n];
       for (int i = 0; i < n; i++)
          u[i] = dist.cdf ((int)v[i]);
       return new DoubleArrayList (u);
   }

   /**
    * Assumes that the real-valued observations @f$U_0,…,U_{n-1}@f$
    * contained in `sortedData` are already sorted in increasing order and
    * computes the differences between the successive observations. Let
    * @f$D@f$ be the differences returned in `spacings`. The difference
    * @f$U_i - U_{i-1}@f$ is put in @f$D_i@f$ for `n1 < i <= n2`, whereas
    * @f$U_{n1} - a@f$ is put into @f$D_{n1}@f$ and @f$b - U_{n2}@f$ is
    * put into @f$D_{n2+1}@f$. The number of observations must be greater
    * or equal than `n2`, we must have `n1 < n2`, and `n1` and `n2` are
    * greater than 0. The size of `spacings` will be at least @f$n+1@f$
    * after the call returns.
    *  @param sortedData   array of sorted observations
    *  @param spacings     pointer to an array object that will be filled
    *                      with spacings
    *  @param n1           starting index, in `sortedData`, of the
    *                      processed observations
    *  @param n2           ending index, in `sortedData` of the processed
    *                      observations
    *  @param a            minimum value of the observations
    *  @param b            maximum value of the observations
    */
   public static void diff (IntArrayList sortedData, IntArrayList spacings,
                            int n1, int n2, int a, int b) {
      if (n1 < 0 || n2 < 0 || n1 >= n2 || n2 >= sortedData.size())
         throw new IllegalArgumentException ("n1 and n2 not valid.");
      int[] u = sortedData.elements();
      int n = sortedData.size();
      if (spacings.size() <= (n2 + 2))
         spacings.setSize (n2 + 2);
      int[] d = spacings.elements();

      d[n1] = u[n1] - a;
      for (int i = n1 + 1; i <= n2; i++)
         d[i] = u[i] - u[i - 1];
      d[n2+1] = b - u[n2];
   }

   /**
    * Same as method  {@link #diff()
    * diff(IntArrayList,IntArrayList,int,int,int,int)}, but for the
    * continuous case.
    *  @param sortedData   array of sorted observations
    *  @param spacings     pointer to an array object that will be filled
    *                      with spacings
    *  @param n1           starting index of the processed observations in
    *                      `sortedData`
    *  @param n2           ending index, in `sortedData` of the processed
    *                      observations
    *  @param a            minimum value of the observations
    *  @param b            maximum value of the observations
    */
   public static void diff (DoubleArrayList sortedData,
                            DoubleArrayList spacings,
                            int n1, int n2, double a, double b) {

      if (n1 < 0 || n2 < 0 || n1 >= n2 || n2 >= sortedData.size())
         throw new IllegalArgumentException ("n1 and n2 not valid.");
      double[] u = sortedData.elements();
      int n = sortedData.size();
      if (spacings.size() <= (n2 + 2))
         spacings.setSize (n2 + 2);
      double[] d = spacings.elements();

      d[n1] = u[n1] - a;
      for (int i = n1 + 1; i <= n2; i++)
         d[i] = u[i] - u[i - 1];
      d[n2+1] = b - u[n2];
   }

   /**
    * Applies one iteration of the *iterated spacings* transformation
    * @cite rKNU98a, @cite tSTE86a&thinsp;. Let @f$U@f$ be the @f$n@f$
    * observations contained into `data`, and let @f$S@f$ be the spacings
    * contained into `spacings`, Assumes that @f$S[0..n]@f$ contains the
    * *spacings* between @f$n@f$ real numbers @f$U_0,…,U_{n-1}@f$ in the
    * interval @f$[0,1]@f$. These spacings are defined by
    * @f[
    *   S_i = U_{(i)} - U_{(i-1)}, \qquad1\le i < n,
    * @f]
    * where @f$U_{(0)}=0@f$, @f$U_{(n-1)}=1@f$, and
    * @f$U_{(0)},…,U_{(n-1)}@f$, are the @f$U_i@f$ sorted in increasing
    * order. These spacings may have been obtained by calling
    * #diff(DoubleArrayList,DoubleArrayList,int,int,double,double). This
    * method transforms the spacings into new spacings, by a variant of
    * the method described in section 11 of @cite rMAR85a&thinsp; and also
    * by Stephens @cite tSTE86a&thinsp;: it sorts @f$S_0,…,S_n@f$ to
    * obtain @f$S_{(0)} \le S_{(1)} \le S_{(2)} \le\cdots\le S_{(n)}@f$,
    * computes the weighted differences
    * @f{align*}{
    *    S_0 
    *    & 
    *   =
    *    (n+1) S_{(0)}, 
    *    \\ 
    *   S_1 
    *    & 
    *   =
    *    n (S_{(1)}-S_{(0)}), 
    *    \\ 
    *   S_2 
    *    & 
    *   =
    *    (n-1) (S_{(2)}-S_{(1)}),
    *    \\  & 
    *    \vdots
    *      \\ 
    *   S_n 
    *    & 
    *   =
    *    S_{(n)}-S_{(n-1)},
    * @f}
    * and computes @f$V_i = S_0 + S_1 + \cdots+ S_i@f$ for @f$0\le i <
    * n@f$. It then returns @f$S_0,…,S_n@f$ in `S[0..n]` and
    * @f$V_1,…,V_n@f$ in `V[1..n]`.
    *
    * Under the assumption that the @f$U_i@f$ are i.i.d. @f$U (0,1)@f$,
    * the new @f$S_i@f$ can be considered as a new set of spacings having
    * the same distribution as the original spacings, and the @f$V_i@f$
    * are a new sample of i.i.d. @f$U (0,1)@f$ random variables, sorted by
    * increasing order.
    *
    * This transformation is useful to detect *clustering* in a data set:
    * A pair of observations that are close to each other is transformed
    * into an observation close to zero. A data set with unusually
    * clustered observations is thus transformed to a data set with an
    * accumulation of observations near zero, which is easily detected by
    * the Anderson-Darling GOF test.
    *  @param data         array of observations
    *  @param spacings     spacings between the observations, will be
    *                      filled with the new spacings
    */
   public static void iterateSpacings (DoubleArrayList data,
                                       DoubleArrayList spacings) {
      if (spacings.size() < (data.size()+1))
         throw new IllegalArgumentException ("Invalid array sizes.");
      double[] v = data.elements();
      spacings.quickSortFromTo (0, data.size());
      double[] s = spacings.elements();
      int n = data.size();

      for (int i = 0; i < n; i++)
         s[n - i] = (i + 1) *  (s[n - i] - s[n - i - 1]);
      s[0] = (n + 1) * s[0];
      v[0] = s[0];
      for (int i = 1; i < n; i++)
         v[i] = v[i - 1] + s[i];
   }

   /**
    * Applies the *power ratios* transformation @f$W@f$ described in
    * section 8.4 of Stephens @cite tSTE86a&thinsp;. Let @f$U@f$ be the
    * @f$n@f$ observations contained into `sortedData`. Assumes that
    * @f$U@f$ contains @f$n@f$ real numbers @f$U_{(0)},…,U_{(n-1)}@f$ from
    * the interval @f$[0,1]@f$, already sorted in increasing order, and
    * computes the transformations:
    * @f[
    *   U’_i = (U_{(i)} / U_{(i+1)})^{i+1}, \qquad i=0,…,n-1,
    * @f]
    * with @f$U_{(n)} = 1@f$. These @f$U’_i@f$ are sorted in increasing
    * order and put back in `U[1...n]`. If the @f$U_{(i)}@f$ are i.i.d.
    * @f$U (0,1)@f$ sorted by increasing order, then the @f$U’_i@f$ are
    * also i.i.d. @f$U (0,1)@f$.
    *
    * This transformation is useful to detect clustering, as explained in
    * #iterateSpacings(DoubleArrayList,DoubleArrayList), except that here
    * a pair of observations close to each other is transformed into an
    * observation close to 1. An accumulation of observations near 1 is
    * also easily detected by the Anderson-Darling GOF test.
    *  @param sortedData   sorted array of real-valued observations in the
    *                      interval @f$[0,1]@f$ that will be overwritten
    *                      with the transformed observations
    */
   public static void powerRatios (DoubleArrayList sortedData) {

      double[] u = sortedData.elements();
      int n = sortedData.size();

      for (int i = 0; i < (n-1); i++) {
         if (u[i + 1] == 0.0 || u[i + 1] == -0.0)
            u[i] = 1.0;
         else
            u[i] = Math.pow (u[i] / u[i + 1], (double) i + 1);
      }

      u[n-1] = Math.pow (u[n-1], (double) n);
      sortedData.quickSortFromTo (0, sortedData.size() - 1);
   }

   /**
    * @}
    */

   /**
    * @name Partitions for the chi-square tests
    * @{
    */

   /**
    * This class helps managing the partitions of possible outcomes into
    * categories for applying chi-square tests. It permits one to
    * automatically regroup categories to make sure that the expected
    * number of observations in each category is large enough. To use this
    * facility, one must first construct an `OutcomeCategoriesChi2` object
    * by passing to the constructor the expected number of observations
    * for each original category. Then, calling the method
    * #regroupCategories will regroup categories in a way that the
    * expected number of observations in each category reaches a given
    * threshold `minExp`. Experts in statistics recommend that `minExp` be
    * always larger than or equal to 5 for the chi-square test to be
    * valid. Thus, `minExp` = 10 is a safe value to use. After the call,
    * `nbExp` gives the expected numbers in the new categories and
    * `loc[i]` gives the relocation of category @f$i@f$, for each @f$i@f$.
    * That is, `loc[i] = j` means that category @f$i@f$ has been merged
    * with category @f$j@f$ because its original expected number was too
    * small, and `nbExp[i]` has been added to `nbExp[j]` and then set to
    * zero. In this case, all observations that previously belonged to
    * category @f$i@f$ are redirected to category @f$j@f$. The variable
    * `nbCategories` gives the final number of categories, `smin` contains
    * the new index of the lowest category, and `smax` the new index of
    * the highest category.
    */
   public static class OutcomeCategoriesChi2 {

      /**
       * Total number of categories.
       */
      public int nbCategories;

      /**
       * Minimum index for valid expected numbers in the array `nbExp`.
       */
      public int smin;

      /**
       * Maximum index for valid expected numbers in the array `nbExp`.
       */
      public int smax;

      /**
       * Expected number of observations for each category.
       */
      public double[] nbExp;

      /**
       * `loc[i]` gives the relocation of the category `i` in the
       * `nbExp` array.
       */
      public int[] loc;

      /**
       * Constructs an `OutcomeCategoriesChi2` object using the array
       * `nbExp` for the number of expected observations in each
       * category. The `smin` and `smax` fields are set to 0 and
       * @f$(n-1)@f$ respectively, where @f$n@f$ is the length of array
       * `nbExp`. The `loc` field is set such that `loc[i]=i` for each
       * `i`. The field `nbCategories` is set to @f$n@f$.
       *  @param nbExp        array of expected observations for each
       *                      category
       */
      public OutcomeCategoriesChi2 (double[] nbExp) {
         this.nbExp = nbExp;
         smin = 0;
         smax = nbExp.length - 1;
         nbCategories = nbExp.length;
         loc = new int[nbExp.length];
         for (int i = 0; i < nbExp.length; i++)
            loc[i] = i;
      }

      /**
       * Constructs an `OutcomeCategoriesChi2` object using the given
       * `nbExp` expected observations array. Only the expected numbers
       * from the `smin` to `smax` (inclusive) indices will be
       * considered valid. The `loc` field is set such that `loc[i]=i`
       * for each `i` in the interval `[smin, smax]`. All `loc[i]` for
       * <tt>i </tt>@f$\le@f$<tt> smin</tt> are set to `smin`, and all
       * `loc[i]` for <tt>i </tt>@f$\ge@f$<tt> smax</tt> are set to
       * `smax`. The field `nbCategories` is set to (<tt>smax - smin +
       * 1</tt>).
       *  @param nbExp        array of expected observations for each
       *                      category
       *  @param smin         Minimum index for valid expected number
       *                      of observations
       *  @param smax         Maximum index for valid expected number
       *                      of observations
       */
      public OutcomeCategoriesChi2 (double[] nbExp, int smin, int smax) {
         this.nbExp = nbExp;
         this.smin = smin;
         this.smax = smax;
         nbCategories = smax - smin + 1;
         loc = new int[nbExp.length];
         for (int i = 0; i < smin; i++)
            loc[i] = smin;
         for (int i = smin; i < smax; i++)
            loc[i] = i;
         for (int i = smax; i < nbExp.length; i++)
            loc[i] = smax;
      }

      /**
       * Constructs an `OutcomeCategoriesChi2` object. The field
       * `nbCategories` is set to `nbCat`.
       *  @param nbExp        array of expected observations for each
       *                      category
       *  @param smin         Minimum index for valid expected number
       *                      of observations
       *  @param smax         Maximum index for valid expected number
       *                      of observations
       *  @param loc          array for which `loc[i]` gives the
       *                      relocation of the category `i`
       */
      public OutcomeCategoriesChi2 (double[] nbExp, int[] loc,
                                    int smin, int smax, int nbCat) {
         this.nbExp = nbExp;
         this.smin = smin;
         this.smax = smax;
         this.nbCategories = nbCat;
         this.loc = loc;
      }

      /**
       * Regroup categories as explained earlier, so that the expected
       * number of observations in each category is at least `minExp`.
       * We usually choose `minExp` = 10.
       *  @param minExp       mininum number of expected observations
       *                      in each category
       */
      public void regroupCategories (double minExp) {
         int s0 = 0, j;
         double somme;

         nbCategories = 0;
         int s = smin;
         while (s <= smax) {
            /* Merge categories to ensure that the number expected
               in each category is >= minExp. */
            if (nbExp[s] < minExp) {
               s0 = s;
               somme = nbExp[s];
               while (somme < minExp && s < smax) {
                  nbExp[s] = 0.0;
                  ++s;
                  somme += nbExp[s];
               }
               nbExp[s] = somme;
               for (j = s0; j <= s; j++)
                  loc[j] = s;

            } else
               loc[s] = s;

            ++nbCategories;
            ++s;
         }
         smin = loc[smin];

         // Special case: the last category, if nbExp < minExp
         if (nbExp[smax] < minExp) {
            if (s0 > smin)
               --s0;
            nbExp[s0] += nbExp[smax];
            nbExp[smax] = 0.0;
            --nbCategories;
            for (j = s0 + 1; j <= smax; j++)
               loc[j] = s0;
            smax = s0;
         }
         if (nbCategories <= 1)
           throw new IllegalStateException ("nbCategories < 2");
         }

      /**
       * Provides a report on the categories.
       *  @return the categories represented as a string
       */
      public String toString() {
         int s, s0;
         double somme;
         final double EPSILON = 5.0E-16;
         StringBuffer sb = new StringBuffer();
         sb.append ("-----------------------------------------------" +
                     PrintfFormat.NEWLINE);
         if (nbExp[smin] < EPSILON)
            sb.append ("Only expected numbers larger than " +
                       PrintfFormat.g (6, 1, EPSILON) + "  are printed" +
                                       PrintfFormat.NEWLINE);
         sb.append ("Number of categories: " +
               PrintfFormat.d (4, nbCategories) + PrintfFormat.NEWLINE +
               "Expected numbers per category:" + PrintfFormat.NEWLINE +
                PrintfFormat.NEWLINE + "Category s      nbExp[s]" +
                PrintfFormat.NEWLINE);

         // Do not print values < EPSILON
         s = smin;
         while (nbExp[s] < EPSILON)
            s++;
         int s1 = s;

         s = smax;
         while (nbExp[s] < EPSILON)
            s--;
         int s2 = s;

         somme = 0.0;
         for (s = s1 ; s <= s2; s++)
            if (loc[s] == s) {
               somme += nbExp[s];
               sb.append (PrintfFormat.d (4, s) + " " +
                          PrintfFormat.f (18, 4, nbExp[s]) +
                          PrintfFormat.NEWLINE);
            }
         sb.append (PrintfFormat.NEWLINE + "Total expected number = " +
                    PrintfFormat.f (18, 2, somme) + PrintfFormat.NEWLINE +
                    PrintfFormat.NEWLINE +
                    "The groupings:" + PrintfFormat.NEWLINE +
                    " Category s      loc[s]" + PrintfFormat.NEWLINE);
         for (s = smin; s <= smax; s++) {
            if ((s == smin) && (s > 0))
               sb.append ("<= ");
            else if ((s == smax) && (s < loc.length - 1))
               sb.append (">= ");
            else
               sb.append ("   ");
            sb.append (PrintfFormat.d (4, s) + " " +
                       PrintfFormat.d (12, loc[s]) + PrintfFormat.NEWLINE);
         }

         sb.append (PrintfFormat.NEWLINE + PrintfFormat.NEWLINE);
         return sb.toString();
      }
   }
   

   /**
    * @}
    */

   /**
    * @name Computing EDF test statistics
    * @{
    */

   /**
    * Computes and returns the chi-square statistic for the observations
    * @f$o_i@f$ in `count[smin...smax]`, for which the corresponding
    * expected values @f$e_i@f$ are in `nbExp[smin...smax]`. Assuming that
    * @f$i@f$ goes from 1 to @f$k@f$, where @f$k =@f$ `smax-smin+1` is the
    * number of categories, the chi-square statistic is defined as
    * @anchor REF_gof_GofStat_eq_chi_square
    * @f[
    *   X^2 = \sum_{i=1}^k \frac{(o_i - e_i)^2}{e_i}. \tag{chi-square}
    * @f]
    * Under the hypothesis that the @f$e_i@f$ are the correct expectations
    * and if these @f$e_i@f$ are large enough, @f$X^2@f$ follows
    * approximately the chi-square distribution with @f$k-1@f$ degrees of
    * freedom. If some of the @f$e_i@f$ are too small, one can use
    * `OutcomeCategoriesChi2` to regroup categories.
    *  @param nbExp        numbers expected in each category
    *  @param count        numbers observed in each category
    *  @param smin         index of the first valid data in `count` and
    *                      `nbExp`
    *  @param smax         index of the last valid data in `count` and
    *                      `nbExp`
    *  @return the @f$X^2@f$ statistic
    */
   public static double chi2 (double[] nbExp, int[] count,
                              int smin, int smax) {
      double diff, khi = 0.0;

      for (int s = smin; s <= smax; s++) {
         if (nbExp[s] <= 0.0) {
            if (count[s] != 0)
              throw new IllegalArgumentException (
                             "nbExp[s] = 0 and count[s] > 0");
         }
         else {
            diff = count[s] - nbExp[s];
            khi += diff * diff / nbExp[s];
         }
      }
      return khi;
   }

   /**
    * Computes and returns the chi-square statistic for the observations
    * @f$o_i@f$ in `count`, for which the corresponding expected values
    * @f$e_i@f$ are in `cat`. This assumes that `cat.regroupCategories`
    * has been called before to regroup categories in order to make sure
    * that the expected numbers in each category are large enough for the
    * chi-square test.
    *  @param cat          numbers expected in each category
    *  @param count        numbers observed in each category
    *  @return the @f$X^2@f$ statistic
    */
   public static double chi2 (OutcomeCategoriesChi2 cat, int[] count) {
      int[] newcount = new int[1 + cat.smax];
      for (int s = cat.smin; s <= cat.smax; s++) {
         newcount[cat.loc[s]] += count[s];
      }

      double diff, khi = 0.0;

      for (int s = cat.smin; s <= cat.smax; s++) {
         if (cat.nbExp[s] > 0.0) {
            diff = newcount[s] - cat.nbExp[s];
            khi += diff * diff / cat.nbExp[s];
         }
      }
      newcount = null;
      return khi;
   }

   /**
    * Computes and returns the chi-square statistic for the observations
    * stored in `data`, assuming that these observations follow the
    * discrete distribution `dist`. For `dist`, we assume that there is
    * one set @f$S=\{a, a+1,…, b-1, b\}@f$, where @f$a<b@f$ and
    * @f$a\ge0@f$, for which @f$p(s)>0@f$ if @f$s\in S@f$ and
    * @f$p(s)=0@f$ otherwise.
    *
    * Generally, it is not possible to divide the integers in intervals
    * satisfying @f$nP(a_0\le s< a_1)=nP(a_1\le s<
    * a_2)=\cdots=nP(a_{j-1}\le s< a_j)@f$ for a discrete distribution,
    * where @f$n@f$ is the sample size, i.e., the number of observations
    * stored into `data`. To perform a general chi-square test, the method
    * starts from `smin` and finds the first non-negligible probability
    * @f$p(s)\ge\epsilon@f$, where @f$\epsilon=@f$
    * DiscreteDistributionInt.EPSILON. It uses `smax` to allocate an array
    * storing the number of expected observations (@f$np(s)@f$) for each
    * @f$s\ge@f$ `smin`. Starting from @f$s=@f$ `smin`, the @f$np(s)@f$
    * terms are computed and the allocated array grows if required until a
    * negligible probability term is found. This gives the number of
    * expected elements for each category, where an outcome category
    * corresponds here to an interval in which sample observations could
    * lie. The categories are regrouped to have at least `minExp`
    * observations per category. The method then counts the number of
    * samples in each categories and calls  #chi2(double[],int[],int,int)
    * to get the chi-square test statistic. If `numCat` is not `null`, the
    * number of categories after regrouping is returned in `numCat[0]`.
    * The number of degrees of freedom is equal to `numCat[0]-1`. We
    * usually choose `minExp` = 10.
    *  @param data         observations, not necessarily sorted
    *  @param dist         assumed probability distribution
    *  @param smin         estimated minimum value of @f$s@f$ for which
    *                      @f$p(s)>0@f$
    *  @param smax         estimated maximum value of @f$s@f$ for which
    *                      @f$p(s)>0@f$
    *  @param minExp       minimum number of expected observations in each
    *                      interval
    *  @param numCat       one-element array that will be filled with the
    *                      number of categories after regrouping
    *  @return the chi-square statistic for a discrete distribution
    */
   public static double chi2 (IntArrayList data, DiscreteDistributionInt dist,
                              int smin, int smax, double minExp, int[] numCat) {
      int i;
      int n = data.size();

      // Find the first non-negligible probability term and fix
      // the real smin.  The linear search starts from the given smin.
      i = smin;
      while (dist.prob (i)*n <= DiscreteDistributionInt.EPSILON)
         i++;
      smin = i--;

      // smax > smin is required
      while (smax <= smin)
         smax = 2*smax + 1;

      // Allocate and fill the array of expected observations
      // Each category s corresponds to a value s for which p(s)>0.
      double[] nbExp = new double[smax+1];
      do {
         i++;
         if (i > smax) {
            smax *= 2;
            double[] newNbExp = new double[smax + 1];
            System.arraycopy (nbExp, smin, newNbExp, smin, nbExp.length - smin);
            nbExp = newNbExp;
         }
         nbExp[i] = dist.prob (i)*n;
      }
      while (nbExp[i] > DiscreteDistributionInt.EPSILON);
      smax = i - 1;

      // Regroup the expected observations intervals
      // satisfying np(s)>=minExp
      OutcomeCategoriesChi2 cat = new OutcomeCategoriesChi2
         (nbExp, smin, smax);
      cat.regroupCategories (minExp);
      if (numCat != null)
         numCat[0] = cat.nbCategories;

      // Count the number of observations in each categories.
      int[] count = new int[cat.smax+1];
      for (i = 0; i < count.length; i++)
         count[i] = 0;
      for (i = 0; i < n; i++) {
         int s = data.get (i);
         while (cat.loc[s] != s)
            s = cat.loc[s];
         count[s]++;
      }

      // Perform the chi-square test
      return chi2 (cat.nbExp, count, cat.smin, cat.smax);
   }

   /**
    * Similar to  #chi2(double[],int[],int,int), except that the expected
    * number of observations per category is assumed to be the same for
    * all categories, and equal to `nbExp`.
    *  @param nbExp        number of expected observations in each
    *                      category (or interval)
    *  @param count        number of counted observations in each category
    *  @param smin         index of the first valid data in `count` and
    *                      `nbExp`
    *  @param smax         index of the last valid data in `count` and
    *                      `nbExp`
    *  @return the @f$X^2@f$ statistic
    */
   public static double chi2Equal (double nbExp, int[] count,
                                   int smin, int smax) {

      double diff, khi = 0.0;
      for (int s = smin; s <= smax; s++) {
         diff = count[s] - nbExp;
         khi += diff * diff;
      }
      return khi / nbExp;
   }

   /**
    * Computes the chi-square statistic for a continuous distribution.
    * Here, the equiprobable case can be used. Assuming that `data`
    * contains observations coming from the uniform distribution, the
    * @f$[0,1]@f$ interval is divided into @f$1/p@f$ subintervals, where
    * @f$p=@f$ <tt>minExp</tt>@f$/n@f$, @f$n@f$ being the sample size,
    * i.e., the number of observations stored in `data`. For each
    * subinterval, the method counts the number of contained observations
    * and the chi-square statistic is computed using
    * #chi2Equal(double,int[],int,int). We usually choose `minExp` = 10.
    *  @param data         array of observations in @f$[0,1)@f$
    *  @param minExp       minimum number of expected observations in each
    *                      subintervals
    *  @return the chi-square statistic for a continuous distribution
    */
   public static double chi2Equal (DoubleArrayList data, double minExp) {
      int n = data.size();
      if (n < (int)Math.ceil (minExp))
         throw new IllegalArgumentException ("Not enough observations");
      double p = minExp/n;
      int m = (int)Math.ceil (1.0/p);
      // to avoid an exception when data[i] = 1/p, reserve one element more
      int[] count = new int[m + 1];
      for (int i = 0; i < n; i++) {
         int j = (int)Math.floor (data.get (i)/p);
         count[j]++;
      }
      // put the elements in count[m] where they belong: in count[m-1]
      count[m - 1] += count[m];
      return chi2Equal (minExp, count, 0, m - 1);
   }

   /**
    * Equivalent to `chi2Equal (data, 10)`.
    *  @param data         array of observations in @f$[0,1)@f$
    *  @return the chi-square statistic for a continuous distribution
    */
   public static double chi2Equal (DoubleArrayList data) {
   return chi2Equal (data, 10.0);
}

   /**
    * Computes and returns the scan statistic @f$S_n (d)@f$, defined in (
    * {@link REF_gof_FBar_eq_scan scan} ). Let @f$U@f$ be
    * the @f$n@f$ observations contained into `sortedData`. The @f$n@f$
    * observations in @f$U[0..n-1]@f$ must be real numbers in the interval
    * @f$[0,1]@f$, sorted in increasing order. (See  FBar.scan for the
    * distribution function of @f$S_n (d)@f$).
    *  @param sortedData   sorted array of real-valued observations in the
    *                      interval @f$[0,1]@f$
    *  @param d            length of the test interval (@f$\in(0,1)@f$)
    *  @return the scan statistic
    */
   public static int scan (DoubleArrayList sortedData, double d) {

      double[] u = sortedData.elements();
      int n = sortedData.size();

      int m = 1, j = 0, i = -1;
      double High = 0.0;

      while (j < (n-1) && High < 1.0) {
         ++i;

         High = u[i] + d;
         while (j < n && u[j] < High)
            ++j;
         // j is now the index of the first obs. to the right of High.
         if (j - i > m)
            m = j - i;
      }
      return m;
   }

   /**
    * Computes and returns the Cramér-von Mises statistic @f$W_n^2@f$ (see
    * @cite tDUR73a, @cite tSTE70a, @cite tSTE86b&thinsp;), defined by
    * @anchor REF_gof_GofStat_eq_CraMis
    * @f[
    *   W_n^2 = \frac{1}{ 12n} + \sum_{j=0}^{n-1} \left(U_{(j)} - \frac{(j+0.5) }{ n}\right)^2, \tag{CraMis}
    * @f]
    * assuming that `sortedData` contains @f$U_{(0)},…,U_{(n-1)}@f$ sorted
    * in increasing order.
    *  @param sortedData   array of sorted real-valued observations in the
    *                      interval @f$[0,1]@f$
    *  @return the Cramér-von Mises statistic
    */
   public static double cramerVonMises (DoubleArrayList sortedData) {
      double w, w2;
      double[] u = sortedData.elements();
      int n = sortedData.size();

      if (n <= 0) {
         System.err.println ("cramerVonMises:  n <= 0");
         return 0.0;
      }

      w2 = 1.0 / (12 * n);
      for (int i = 0; i < n; i++) {
         w = u[i] - (i + 0.5) / n;
         w2 += w * w;
      }
      return w2;
   }

   /**
    * Computes and returns the Watson statistic @f$G_n@f$ (see
    * @cite tWAT76a, @cite tDAR83a&thinsp;), defined by
    * @anchor REF_gof_GofStat_eq_WatsonG
    * @f{align}{
    *    G_n 
    *    & 
    *   =
    *    \sqrt{n} \max_{\Rule{0.0pt}{7.0pt}{0.0pt} 0\le j \le n-1} \left\{ (j+1)/n - U_{(j)} + \overline{U}_n - 1/2 \right\} \tag{WatsonG} 
    *    \\  & 
    *   =
    *    \sqrt{n}\left(D_n^+ + \overline{U}_n - 1/2\right), \nonumber
    * @f}
    * where @f$\overline{U}_n@f$ is the average of the observations
    * @f$U_{(j)}@f$, assuming that `sortedData` contains the sorted
    * @f$U_{(0)},…,U_{(n-1)}@f$.
    *  @param sortedData   array of sorted real-valued observations in the
    *                      interval @f$[0,1]@f$
    *  @return the Watson statistic @f$G_n@f$
    */
   public static double watsonG (DoubleArrayList sortedData) {
      double[] u = sortedData.elements();
      int n = sortedData.size();
      double sumZ;
      double d2;
      double dp, g;
      double unSurN = 1.0 / n;

      if (n <= 0) {
         System.err.println ("watsonG: n <= 0");
         return 0.0;
      }

      // degenerate case n = 1
      if (n == 1)
         return 0.0;

      // We assume that u is already sorted.
      dp = sumZ = 0.0;
      for (int i = 0; i < n; i++) {
         d2 = (i + 1) * unSurN - u[i];
         if (d2 > dp)
            dp = d2;
         sumZ += u[i];
      }
      sumZ = sumZ * unSurN - 0.5;
      g = Math.sqrt ((double) n) * (dp + sumZ);
      return g;
   }

   /**
    * Computes and returns the Watson statistic @f$U_n^2@f$ (see
    * @cite tDUR73a, @cite tSTE70a, @cite tSTE86b&thinsp;), defined by
    * @anchor REF_gof_GofStat_eq_WatsonU
    * @f{align}{
    *    W_n^2 
    *    & 
    *   =
    *    \frac{1}{ 12n} + \sum_{j=0}^{n-1} \left\{U_{(j)} - \frac{(j + 0.5)}{ n} \right\}^2, 
    *    \\ 
    *   U_n^2 
    *    & 
    *   =
    *    W_n^2 - n\left(\overline{U}_n - 1/2\right)^2. \tag{WatsonU}
    * @f}
    * where @f$\overline{U}_n@f$ is the average of the observations
    * @f$U_{(j)}@f$, assuming that `sortedData` contains the sorted
    * @f$U_{(0)},…,U_{(n-1)}@f$.
    *  @param sortedData   array of sorted real-valued observations in the
    *                      interval @f$[0,1]@f$
    *  @return the Watson statistic @f$U_n^2@f$
    */
   public static double watsonU (DoubleArrayList sortedData) {
      double sumZ, w, w2, u2;
      double[] u = sortedData.elements();
      int n = sortedData.size();

      if (n <= 0) {
         System.err.println ("watsonU: n <= 0");
         return 0.0;
      }

      // degenerate case n = 1
      if (n == 1)
         return 1.0 / 12.0;

      sumZ = 0.0;
      w2 = 1.0 / (12 * n);
      for (int i = 0; i < n; i++) {
         sumZ += u[i];
         w = u[i] - (i + 0.5) / n;
         w2 += w * w;
      }
      sumZ = sumZ / n - 0.5;
      u2 = w2 - sumZ * sumZ * n;
      return u2;
   }

   /**
    * Used by  #andersonDarling(DoubleArrayList). <tt>Num.DBL_EPSILON</tt>
    * is usually @f$2^{-52}@f$.
    */
   public static double EPSILONAD = Num.DBL_EPSILON/2;

   /**
    * Computes and returns the Anderson-Darling statistic @f$A_n^2@f$ (see
    * method  #andersonDarling(double[]) ).
    *  @param sortedData   array of sorted real-valued observations in the
    *                      interval @f$[0,1]@f$
    *  @return the Anderson-Darling statistic
    */
   public static double andersonDarling (DoubleArrayList sortedData) {
      double[] v = sortedData.elements();
      return andersonDarling (v);
   }

   /**
    * Computes and returns the Anderson-Darling statistic @f$A_n^2@f$ (see
    * @cite tLEW61a, @cite tSTE86b, @cite tAND52a&thinsp;), defined by
    * @anchor REF_gof_GofStat_eq_Andar
    * @f{align*}{
    *    A_n^2 
    *    & 
    *   =
    *    -n -\frac{1}{ n} \sum_{j=0}^{n-1} \left\{ (2j+1)\ln(U_{(j)}) + (2n-1-2j) \ln(1-U_{(j)}) \right\}, \tag{Andar}
    * @f}
    * assuming that `sortedData` contains @f$U_{(0)},…,U_{(n-1)}@f$ sorted
    * in increasing order.  When computing @f$A_n^2@f$, all observations
    * @f$U_i@f$ are projected on the interval @f$[\epsilon,
    * 1-\epsilon]@f$ for some @f$\epsilon> 0@f$, in order to avoid
    * numerical overflow when taking the logarithm of @f$U_i@f$ or
    * @f$1-U_i@f$. The variable `EPSILONAD` gives the value of
    * @f$\epsilon@f$.
    *  @param sortedData   array of sorted real-valued observations in the
    *                      interval @f$[0,1]@f$
    *  @return the Anderson-Darling statistic
    */
   public static double andersonDarling (double[] sortedData) {
      double u1;
      double u, a2;
      int n = sortedData.length;

      if (n <= 0) {
         System.err.println ("andersonDarling: n <= 0");
         return 0.0;
      }

      a2 = 0.0;
      for (int i = 0; i < n; i++) {
         u = sortedData[i];
         u1 = 1.0 - u;
         if (u < EPSILONAD)
            u = EPSILONAD;
         else if (u1 < EPSILONAD)
            u1 = EPSILONAD;
         a2 += (2*i + 1)*Math.log (u) + (1 + 2*(n - i - 1))*
                    Math.log (u1);
      }
      a2 = -n - a2 / n;
      return a2;
   }

   /**
    * Computes the Anderson-Darling statistic @f$A_n^2@f$ and the
    * corresponding @f$p@f$-value @f$p@f$. The @f$n@f$ (unsorted)
    * observations in `data` are assumed to be independent and to come
    * from the continuous distribution `dist`. Returns the 2-elements
    * array [@f$A_n^2@f$, @f$p@f$].
    *  @param data         array of observations
    *  @param dist         assumed distribution of the observations
    *  @return the array @f$[A_n^2@f$, @f$p]@f$.
    */
   public static double[] andersonDarling (double[] data,
                                           ContinuousDistribution dist) {
      int n = data.length;
      double[] U = new double[n];
      for (int i = 0; i < n; i++) {
         U[i] = dist.cdf(data[i]);
      }

      Arrays.sort(U);
      double x = GofStat.andersonDarling(U);
      double v = AndersonDarlingDistQuick.barF(n, x);
      double[] res = {x, v};
      return res;
   }

   /**
    * Computes the Kolmogorov-Smirnov (KS) test statistics @f$D_n^+@f$,
    * @f$D_n^-@f$, and @f$D_n@f$ (see method
    * #kolmogorovSmirnov(DoubleArrayList) ). Returns the array
    * [@f$D_n^+@f$, @f$D_n^-@f$, @f$D_n@f$].
    *  @param sortedData   array of sorted real-valued observations in the
    *                      interval @f$[0,1]@f$
    *  @return the array [@f$D_n^+@f$, @f$D_n^-@f$, @f$D_n@f$]
    */
   public static double[] kolmogorovSmirnov (double[] sortedData) {
      DoubleArrayList v = new DoubleArrayList(sortedData);
      return kolmogorovSmirnov (v);
   }

   /**
    * Computes the Kolmogorov-Smirnov (KS) test statistics @f$D_n^+@f$,
    * @f$D_n^-@f$, and @f$D_n@f$ defined by
    * @anchor REF_gof_GofStat_eq_DNp
    * @anchor REF_gof_GofStat_eq_DNm
    * @anchor REF_gof_GofStat_eq_DN
    * @f{align}{
    *    D_n^+ 
    *    & 
    *   =
    *    \max_{0\le j\le n-1} \left((j+1)/n - U_{(j)}\right), \tag{DNp} 
    *    \\ 
    *   D_n^- 
    *    & 
    *   =
    *    \max_{0\le j\le n-1} \left(U_{(j)} - j/n\right), \tag{DNm} 
    *    \\ 
    *   D_n 
    *    & 
    *   =
    *    \max (D_n^+, D_n^-). \tag{DN}
    * @f}
    * and returns an array of length 3 that contains [@f$D_n^+@f$,
    * @f$D_n^-@f$, @f$D_n@f$]. These statistics compare the empirical
    * distribution of @f$U_{(1)},…,U_{(n)}@f$, which are assumed to be in
    * `sortedData`, with the uniform distribution over @f$[0,1]@f$.
    *  @param sortedData   array of sorted real-valued observations in the
    *                      interval @f$[0,1]@f$
    *  @return the array [@f$D_n^+@f$, @f$D_n^-@f$, @f$D_n@f$]
    */
   public static double[] kolmogorovSmirnov (DoubleArrayList sortedData) {
      double[] ret = new double[3];
      int n = sortedData.size();

      if (n <= 0) {
         ret[0] = ret[1] = ret[2] = 0.0;
         System.err.println ("kolmogorovSmirnov:   n <= 0");
         return ret;
      }

      double[] retjo = kolmogorovSmirnovJumpOne (sortedData, 0.0);
      ret[0] = retjo[0];
      ret[1] = retjo[1];
      if (ret[1] > ret[0])
         ret[2] = ret[1];
      else
         ret[2] = ret[0];

      return ret;
   }

   /**
    * Computes the KolmogorovSmirnov (KS) test statistics and their
    * @f$p@f$-values. This is to compare the empirical distribution of the
    * (unsorted) observations in `data` with the theoretical distribution
    * `dist`. The KS statistics @f$D_n^+@f$, @f$D_n^-@f$ and @f$D_n@f$ are
    * returned in `sval[0]`, `sval[1]`, and `sval[2]` respectively, and
    * their corresponding @f$p@f$-values are returned in `pval[0]`,
    * `pval[1]`, and `pval[2]`.
    *  @param data         array of observations to be tested
    *  @param dist         assumed distribution of the observations
    *  @param sval         values of the 3 KS statistics
    *  @param pval         @f$p@f$-values for the 3 KS statistics
    */
   public static void kolmogorovSmirnov (double[] data,
                                         ContinuousDistribution dist,
                                         double[] sval,
                                         double[] pval) {
      int n = data.length;
      double[] T = new double[n];
      for (int i = 0; i < n; i++) {
         T[i] = dist.cdf (data[i]);
      }

      Arrays.sort (T);
      double[] statks = GofStat.kolmogorovSmirnov (T);
      for (int i = 0; i < 3; i++) {
         sval[i] = statks[i];
      }
      pval[2] = KolmogorovSmirnovDistQuick.barF (n, sval[2]);
      pval[1] = KolmogorovSmirnovPlusDist.barF (n, sval[1]);
      pval[0] = KolmogorovSmirnovPlusDist.barF (n, sval[0]);
   }

   /**
    * Compute the KS statistics @f$D_n^+(a)@f$ and @f$D_n^-(a)@f$ defined
    * in the description of the method
    * FDist.kolmogorovSmirnovPlusJumpOne, assuming that @f$F@f$ is the
    * uniform distribution over @f$[0,1]@f$ and that
    * @f$U_{(1)},…,U_{(n)}@f$ are in `sortedData`. Returns the array
    * [@f$D_n^+@f$, @f$D_n^-@f$].
    *  @param sortedData   array of sorted real-valued observations in the
    *                      interval @f$[0,1]@f$
    *  @param a            size of the jump
    *  @return the array [@f$D_n^+@f$, @f$D_n^-@f$]
    */
   public static double[] kolmogorovSmirnovJumpOne (DoubleArrayList sortedData,
                                                    double a) {
      /* Statistics KS+ and KS-. Case with 1 jump at a, near the lower tail of
         the distribution. */

      double[] u = sortedData.elements();
      int n = sortedData.size();
      int j, i;
      double d2, d1, unSurN;
      double[] ret = new double[2];

      if (n <= 0) {
         ret[0] = ret[1] = 0.0;
         System.err.println ("kolmogorovSmirnovJumpOne: n <= 0");
         return ret;
      }

      ret[0] = 0.0;
      ret[1] = 0.0;
      unSurN = 1.0 / n;
      j = 0;

      while (j < n && u[j] <= a + EPSILOND) ++j;

      for (i = j - 1; i < n; i++) {
         if (i >= 0) {
            d1 = (i + 1) * unSurN - u[i];
            if (d1 > ret[0])
               ret[0] = d1;
         }
         if (i >= j) {
            d2 = u[i] - i * unSurN;
            if (d2 > ret[1])
               ret[1] = d2;
         }
      }
      return ret;
   }

   /**
    * Computes a variant of the @f$p@f$-value @f$p@f$ whenever a test
    * statistic has a *discrete* probability distribution. This
    * @f$p@f$-value is defined as follows:
    * @f{align*}{
    *    p_L 
    *    & 
    *    = 
    *    P[Y \le y] 
    *    \\ 
    *   p_R 
    *    & 
    *    = 
    *    P[Y \ge y] 
    *    \\ 
    *    p 
    *    & 
    *    = 
    *    \left\{ 
    *   \begin{array}{l@{qquad}l}
    *    p_R, 
    *    & 
    *    \mbox{if } p_R < p_L 
    *    \\ 
    *    1 - p_L, 
    *    \mbox{if } p_R \ge p_L \mbox{ and } p_L < 0.5 
    *    \\ 
    *    0.5 
    *    & 
    *    \mbox{otherwise.} 
    *   \end{array}
    *    \right.
    * @f}
    * @f[
    *   \begin{array}{rll}
    *    p =
    *    & 
    *    p_R, 
    *    & 
    *   \qquad\mbox{if } p_R < p_L, 
    *    \\ 
    *   p =
    *    & 
    *    1 - p_L, 
    *    & 
    *   \qquad\mbox{if } p_R \ge p_L \mbox{ and } p_L < 0.5, 
    *    \\ 
    *   p =
    *    & 
    *    0.5 
    *    & 
    *   \qquad\mbox{otherwise.} 
    *   \end{array}
    * @f]
    * The function takes @f$p_L@f$ and @f$p_R@f$ as input and returns
    * @f$p@f$.
    *  @param pL           left @f$p@f$-value
    *  @param pR           right @f$p@f$-value
    *  @return the @f$p@f$-value for a test on a discrete distribution
    */
   public static double pDisc (double pL, double pR) {
      double p;

      if (pR < pL)
         p = pR;
      else if (pL > 0.5)
         p = 0.5;
      else
         p = 1.0 - pL;
      // Note: si p est tres proche de 1, on perd toute la precision ici!
      // Note2: je ne pense pas que cela puisse se produire a cause des if (RS)
      return p;
   }
}

/**
 * @}
 */
