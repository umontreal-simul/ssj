/*
 * Class:        ListOfTalliesWithCovariance
 * Description:  List of tallies with covariance
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Éric Buist 
 * @since        2007
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
package umontreal.ssj.stat.list;

import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.TallyStore;
import java.util.logging.Level;
import java.util.logging.Logger;
import cern.colt.matrix.DoubleMatrix1D;

/**
 * Extends  @ref ListOfTallies to add support for the computation of the
 * sample covariance between each pair of elements in a list, without storing
 * all observations. This list of tallies contains internal structures to
 * keep track of @f$\bar{X}_{n, i}@f$ for @f$i=0, …, d-1@f$, and
 * @f$\sum_{k=0}^{n-1} (X_{i, k} - \bar{X}_{k, i})(X_{j, k} - \bar{X}_{k,
 * j})/n@f$, for @f$i=0,…, d-2@f$ and @f$j=1,…,d-1@f$, with @f$j>i@f$. Here,
 * @f$\bar{X}_{n, i}@f$ is the @f$i@f$th component of
 * @f$\bar{\mathbf{X}}_n@f$, the average vector, and @f$\bar{X}_{0, i}=0@f$
 * for @f$i=0,…,d-1@f$. The value @f$X_{i,k}@f$ corresponds to the @f$i@f$th
 * component of the @f$k@f$th observation @f$\mathbf{X}_k@f$. These sums are
 * updated every time a vector is added to this list, and are used to
 * estimate the covariances.
 *
 * Note: the size of the list of tallies must remain fixed because of the
 * data structures used for computing sample covariances. As a result, the
 * first call to `init` makes this list unmodifiable.
 *
 * Note: for the sample covariance to be computed between a pair of tallies,
 * the number of observations in each tally should be the same. It is
 * therefore recommended to always add complete vectors of observations to
 * this list. Moreover, one must use the  #add method in this class to add
 * vectors of observations for the sums used for covariance estimation to be
 * updated correctly. Failure to use this method, e.g., adding observations
 * to each individual tally in the list, will result in an incorrect estimate
 * of the covariances, unless the tallies in the list can store observations.
 * For example, the following code, which adds the vector `v` in the list of
 * tallies `list`, works correctly only if the list contains instances of
 * @ref umontreal.ssj.stat.TallyStore:
 *
 * @code
 *
 *    for (int i = 0; i < v.length; i++)
 *       list.get (i).add (v[i]);
 *
 * @endcode
 *
 *  But the following code is always correct:
 *
 * @code
 *
 *    list.add (v);
 *
 * @endcode
 *
 * <div class="SSJ-bigskip"></div>
 */
public class ListOfTalliesWithCovariance<E extends Tally>
       extends ListOfTallies<E> {
   private double[] tempArray;
   private double[][] sxy;

   // Determines if we use a numerically stable covariance
   // formula.
   private boolean isStable = true;

   // The average of the first observations, for each tally
   private double[] curAverages;

   // The sum (xi - average)(yi - average) of the first observations
   private double[][] curSum2;
   private Logger log = Logger.getLogger ("umontreal.ssj.stat.list");

   /**
    * Creates an empty list of tallies with covariance support. One must
    * fill the list with tallies, and call  #init before adding any
    * observation.
    */
   public ListOfTalliesWithCovariance() {
      super();
   }

   /**
    * Creates an empty list of tallies with covariance support and name
    * `name`. One must fill the list with tallies, and call  #init before
    * adding any observation.
    *  @param name         the name of the new list.
    */
   public ListOfTalliesWithCovariance (String name) {
      super (name);
   }

   /**
    * This factory method constructs and returns a list of tallies with
    * `size` instances of  @ref umontreal.ssj.stat.Tally.
    *  @param size         the size of the list.
    *  @return the created list.
    */
   public static ListOfTalliesWithCovariance<Tally> createWithTally (int size) {
      ListOfTalliesWithCovariance<Tally> list = new ListOfTalliesWithCovariance<Tally>();
      for (int i = 0; i < size; i++)
         list.add (new Tally());
      list.init();
      return list;
   }

   /**
    * This factory method constructs and returns a list of tallies with
    * `size` instances of  @ref umontreal.ssj.stat.TallyStore.
    *  @param size         the size of the list.
    *  @return the created list.
    */
   public static ListOfTalliesWithCovariance<TallyStore> createWithTallyStore
                                             (int size) {
      ListOfTalliesWithCovariance<TallyStore> list = new ListOfTalliesWithCovariance<TallyStore>();
      for (int i = 0; i < size; i++)
         list.add (new TallyStore());
      list.init();
      return list;
   }


   private void createSxy() {
      int l = size();
      if (isStable) {
         curAverages = new double[l];
         curSum2 = new double[l-1][];
         for (int i = 0; i < l - 1; i++)
            curSum2[i] = new double[l - 1 - i];
      }
      else {
         sxy = new double[l - 1][];
         for (int i = 0; i < l - 1; i++)
            sxy[i] = new double[l - 1 - i];
      }
      tempArray = new double[l];
   }

   public void init() {
      super.init();

      if (isModifiable()) {
         setUnmodifiable();
         createSxy();
      }
      if (isStable) {
         for (int i = 0; i < curAverages.length; i++)
            curAverages[i] = 0;
         for (int i = 0; i < curSum2.length; i++)
            for (int j = 0; j < curSum2[i].length; j++)
               curSum2[i][j] = 0;
      }
      else
         for (int i = 0; i < sxy.length; i++)
            for (int j = 0; j < sxy[i].length; j++)
               sxy[i][j] = 0;
   }

/**
 * Adds a new vector of observations `x` to this list of tallies, and updates
 * the internal data structures computing averages, and sums of products. One
 * must use this method instead of adding observations to individual tallies
 * to get a covariance estimate.
 *  @param x            the new vector of observations.
 */
public void add (double[] x) {
      int l = size();

      int structSize = 0;
      structSize = (isStable) ? curSum2.length : sxy.length;
      if (structSize != l - 1)
            throw new IllegalArgumentException ("The structure's size mismatches the list's size");

      super.add (x);
      if (isStable) {
         int numObs = get (0).numberObs();
         // get (i1).average() would return the average over the n
         // observations, but we need the average over the last n-1 observations.
         for (int i1 = 0; i1 < l - 1; i1++)
            for (int i2 = i1 + 1; i2 < l; i2++)
               curSum2[i1][i2 - i1 - 1] += (numObs - 1)*
                  (x[i1] - curAverages[i1])*(x[i2] - curAverages[i2])/numObs;
         for (int i = 0; i < l; i++)
            curAverages[i] += (x[i] - curAverages[i])/numObs;
         // Now, curAverages[i] == get (i).average()
      }
      else
         for (int i1 = 0; i1 < l - 1; i1++)
            for (int i2 = i1 + 1; i2 < l; i2++)
               sxy[i1][i2 - i1 - 1] += x[i1]*x[i2];
   }


   public void add (DoubleMatrix1D x) {
      x.toArray (tempArray);
      add (tempArray);
   }

   public double covariance (int i, int j) {
      if (i == j)
         return get (i).variance();
      if (i > j) {
         // Make sure that i1 < i2, to have a single case
         int tmp = i;
         i = j;
         j = tmp;
      }

      Tally tallyi = get (i);
      Tally tallyj = get (j);
      if (tallyi == null || tallyj == null)
         return Double.NaN;
      int n = tallyi.numberObs();
      if (n != tallyj.numberObs()) {
         log.logp (Level.WARNING, "ListOfTalliesWithCovariance", "covariance",
            "Tally " + i  + ", with name " + tallyi.getName() + ", contains " 
            + n + " observations while " +
              "tally " + j + ", with name " + tallyj.getName() + ", contains " + tallyj.numberObs() + "observations");
         return Double.NaN;
      }

      if (n < 2) {
         log.logp (Level.WARNING, "ListOfTalliesWithCovariance", "covariance",
            "Tally " + i + ", with name " + tallyi.getName() + ", contains " + n + " observation");
         return Double.NaN;
      }
      if (tallyi instanceof TallyStore && tallyj instanceof TallyStore)
         return ((TallyStore) tallyi).covariance ((TallyStore) tallyj);
      else if (isStable)
         return curSum2[i][j - i - 1]/(n-1);
      else {
         double sum1 = tallyi.sum();
         double sum2 = tallyj.sum();
         double sum12 = sxy[i][j - i - 1];
         return (sum12 - sum1*sum2/n)/(n-1);
      }
   }

/**
 * Clones this object. This clones the list of tallies and the data
 * structures holding the sums of products but not the tallies comprising the
 * list. The created clone is modifiable, even though the original list is
 * unmodifiable.
 */
public ListOfTalliesWithCovariance<E> clone() {
      ListOfTalliesWithCovariance<E> ta = (ListOfTalliesWithCovariance<E>)super.clone();
      ta.tempArray = new double[size()];
      if (curAverages != null)
         ta.curAverages = curAverages.clone();
      if (sxy != null) {
         ta.sxy = new double[sxy.length][];
         for (int i = 0; i < sxy.length; i++)
            ta.sxy[i] = sxy[i].clone();
      }
      if (curSum2 != null) {
         ta.curSum2 = new double[curSum2.length][];
         for (int i = 0; i < curSum2.length; i++)
            ta.curSum2[i] = curSum2[i].clone();
      }
      return ta;
   }
}