/*
 * Class:        ListOfTallies
 * Description:  List of statistical collectors.
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

import umontreal.ssj.util.PrintfFormat;
import cern.colt.list.DoubleArrayList;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.TallyStore;

/**
 * Represents a list of tally statistical collectors. Each element of the
 * list is an instance of  @ref umontreal.ssj.stat.Tally, and a vector of
 * observations can be added with the  #add(double[]) method. This class
 * defines factory methods to fill a newly-constructed list with `Tally` or
 * `TallyStore` instances.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class ListOfTallies<E extends Tally> extends ListOfStatProbes<E> {

   /**
    * Constructs a new empty list of tallies.
    */
   public ListOfTallies() {
      super();
   }

   /**
    * Constructs a new empty list of tallies with name `name`.
    *  @param name         the name of the new list.
    */
   public ListOfTallies (String name) {
      super (name);
   }

   /**
    * This factory method constructs and returns a list of tallies with
    * `size` instances of  @ref umontreal.ssj.stat.Tally.
    *  @param size         the size of the list.
    *  @return the created list.
    */
   public static ListOfTallies<Tally> createWithTally (int size) {
      ListOfTallies<Tally> list = new ListOfTallies<Tally>();
      for (int i = 0; i < size; i++)
         list.add (new Tally());
      return list;
   }

   /**
    * This factory method constructs and returns a list of tallies with
    * `size` instances of  @ref umontreal.ssj.stat.TallyStore.
    *  @param size         the size of the list.
    *  @return the created list.
    */
   public static ListOfTallies<TallyStore> createWithTallyStore (int size) {
      ListOfTallies<TallyStore> list = new ListOfTallies<TallyStore>();
      for (int i = 0; i < size; i++)
         list.add (new TallyStore());
      return list;
   }
   
   /**
    * This factory method constructs and returns a list of tallies with
    * `size` instances of  @ref umontreal.ssj.stat.TallyStore, each of size \f$t\f$.
    * @param size the size of the list.
    * @param t size of each TallyStore in this list.
    * @return the created list.
    */
   
   public static ListOfTallies<TallyStore> createWithTallyStore (int size,int t) {
	      ListOfTallies<TallyStore> list = new ListOfTallies<TallyStore>();
	      for (int i = 0; i < t; i++)
	         list.add (new TallyStore(size));
	      return list;
	   }


   /**
    * Adds the observation `x[i]` in tally `i` of this list, for <tt>i =
    * 0,…, size() - 1</tt>. No observation is added if the value is
    * `Double.NaN`, or if collecting is turned OFF. If broadcasting is ON,
    * the given array is notified to all registered observers. The given
    * array `x` not being stored by this object, it can be freely used and
    * modified after the call to this method.
    *  @param x            the array of observations.
    *  @exception NullPointerException if `x` is `null`.
    *  @exception IllegalArgumentException if the length of `x` does not
    * correspond to `size()`.
    */
   public void add (double[] x) {
      int l = size();
      if (x.length != l)
         throw new IllegalArgumentException
            ("Incompatible array length: given " +
            x.length + ", required " + l);
      if (collect)
         for (int i = 0; i < l; i++) {
            double v = x[i];
            Tally ta = get (i);
            if (!Double.isNaN (v) && ta != null)
               ta.add (v);
         }
         notifyListeners (x);
   }

   /**
    * Assuming that each tally in this list contains the same number of
    * observations, returns the number of observations in tally&nbsp;0, or
    * 0 if this list is empty.
    *  @return the number of observations.
    */
   public int numberObs() {
      if (size() == 0)
         return 0;
      Tally t0 = get (0);
      return t0 == null ? 0 : t0.numberObs();
   }

   /**
    * Tests that every tally in this list contains the same number of
    * observations. This returns `true` if and only if all tallies have
    * the same number of observations, or if this list is empty. If
    * observations are always added using the  #add(double[]) method from
    * this class, and not  umontreal.ssj.stat.Tally.add(double) from
    * @ref umontreal.ssj.stat.Tally, this method always returns `true`.
    *  @return the success indicator of the test.
    */
   public boolean areAllNumberObsEqual() {
      final int l = size();
      int n = numberObs();
      for (int i = 1; i < l; i++) {
         Tally t = get (i);
         if (t.numberObs() != n)
            return false;
      }
      return true;
   }

   /**
    * Computes the average for each tally in this list, and stores the
    * averages in the array `r`. If the tally `i` has no observation, the
    * `Double.NaN` value is stored in the array, at index&nbsp;`i`.
    */
   public void average (double[] r) {
      final int l = size();
      for (int i = 0; i < l; i++) {
          // Manual check to avoid repetitive logs when all tallies
          // have 0 observation.
         Tally ta = get (i);
         double v = ta == null || ta.numberObs() == 0 ? Double.NaN : ta.average();
         r[i] = v;
      }
   }

   /**
    * For each tally in this list, computes the sample variance, and
    * stores the variances into the array `v`. If, for some
    * tally&nbsp;`i`, there are not enough observations for estimating the
    * variance, `Double.NaN` is stored in the array.
    *  @param v            the array to be filled with sample variances.
    *  @exception NullPointerException if `v` is `null`.
    *  @exception IllegalArgumentException if `v.length` does not
    * correspond to  umontreal.ssj.stat.list.ListOfStatProbes.size.
    */
   public void variance (double[] v) {
      if (size() != v.length)
         throw new IllegalArgumentException
            ("Invalid length of given array");
      for (int i = 0; i < v.length; i++) {
         Tally tally = get (i);
         if (tally == null || tally.numberObs() < 2)
            v[i] = Double.NaN;
         else
            v[i] = tally.variance();
      }
   }

   /**
    * For each tally in this list, computes the sample standard deviation,
    * and stores the standard deviations into the array `std`. This is
    * equivalent to calling  #variance(double[]) and performing a square
    * root on every element of the filled array.
    *  @param std          the array to be filled with standard
    *                      deviations.
    *  @exception NullPointerException if `std` is `null`.
    *  @exception IllegalArgumentException if `std.length` does not
    * correspond to `size()`.
    */
   public void standardDeviation (double[] std) {
      if (size() != std.length)
         throw new IllegalArgumentException
            ("Invalid length of given array");
      for (int i = 0; i < std.length; i++) {
         Tally tally = get (i);
         if (tally == null || tally.numberObs() < 2)
            std[i] = Double.NaN;
         else
            std[i] = tally.standardDeviation();
      }
   }

   /**
    * Returns the empirical covariance of the observations in tallies with
    * indices `i` and `j`. If @f$x_1,…,x_n@f$ represent the observations
    * in tally `i` whereas @f$y_1,…,y_n@f$ represent the observations in
    * tally `j`, then the covariance is given by
    * @f[
    *   S_{X, Y} = \frac{1}{n-1}\sum_{k=1}^n (x_k - \bar{X}_n)(y_k - \bar{Y}_n) = \frac{1}{n-1}\left(\sum_{k=1}^n x_ky_k - \frac{1}{n} \sum_{k=1}^n x_k\sum_{r=1}^n y_r\right).
    * @f]
    * This returns `Double.NaN` if the tallies do not contain the same
    * number of observations, or if they contain less than two
    * observations. This method throws an exception if the underlying
    * tallies are not capable of storing observations, i.e. if the tallies
    * are not TallyStores. The  @ref ListOfTalliesWithCovariance subclass
    * provides an alternative implementation of this method which does not
    * require the observations to be stored.
    *  @param i            the index of the first tally.
    *  @param j            the index of the second tally.
    *  @return the value of the covariance.
    *
    *  @exception ArrayIndexOutOfBoundsException if one or both indices
    * are out of bounds.
    */
   public double covariance (int i, int j) {
      if (i == j)
         return get (i).variance();

      TallyStore tallyi = (TallyStore)get (i);
      TallyStore tallyj = (TallyStore)get (j);
      return tallyi.covariance (tallyj);
   }

   /**
    * Returns the empirical correlation between the observations in
    * tallies with indices `i` and `j`. If the tally `i` contains a sample
    * of the random variate @f$X@f$ and the tally `j` contains a sample of
    * @f$Y@f$, this corresponds to
    * @f[
    *   \mathrm{Cor}(X, Y)=\mathrm{Cov}(X, Y) /\sqrt{\mathrm{Var}(X)\mathrm{Var}(Y)}.
    * @f]
    * This method uses  #covariance(int,int) to obtain an estimate of the
    * covariance, and  umontreal.ssj.stat.Tally.variance in class
    * @ref umontreal.ssj.stat.Tally to obtain the sample variances.
    *  @param i            the index of the first tally.
    *  @param j            the index of the second tally.
    *  @return the value of the correlation.
    *
    *  @exception ArrayIndexOutOfBoundsException if one or both indices
    * are out of bounds.
    */
   public double correlation (int i, int j) {
      if (i == j)
         return 1.0;
      double cov = covariance (i, j);
      Tally tallyi = get (i);
      Tally tallyj = get (j);
      if (tallyi == null || tallyj == null)
         return Double.NaN;
      return cov/Math.sqrt (tallyi.variance()*tallyj.variance());
   }

   /**
    * Constructs and returns the sample covariance matrix for the tallies
    * in this list. The given @f$d\times d@f$ matrix `c`, where
    * @f$d=@f$&nbsp;`size()`, is filled with the computed sample
    * covariances. Element `c.get (i, j)` corresponds to the result of
    * `covariance (i, j)`.
    *  @param c            the matrix to be filled with the sample
    *                      covariances.
    *  @exception NullPointerException if `c` is `null`.
    *  @exception IllegalArgumentException if the number of rows or
    * columns in `c` does not correspond to `size()`.
    */
   public void covariance (DoubleMatrix2D c) {
      int l = size();
      if (c.rows() != l)
         throw new IllegalArgumentException
            ("Invalid number of rows in covariance matrix");
      if (c.columns() != l)
         throw new IllegalArgumentException
            ("Invalid number of columns in covariance matrix");
      for (int i1 = 0; i1 < l; i1++)
         c.setQuick (i1, i1, get (i1).variance());
      for (int i1 = 0; i1 < l - 1; i1++)
         for (int i2 = i1 + 1; i2 < l; i2++) {
            double cov = covariance (i1, i2);
            c.setQuick (i1, i2, cov);
            c.setQuick (i2, i1, cov);
         }
   }

   /**
    * Similar to  #covariance(DoubleMatrix2D) for computing the sample
    * correlation matrix.
    *  @param c            the matrix to be filled with the correlations.
    *  @exception NullPointerException if `c` is `null`.
    *  @exception IllegalArgumentException if the number of rows or
    * columns in `c` does not correspond to
    * umontreal.ssj.stat.list.ListOfStatProbes.size.
    */
   public void correlation (DoubleMatrix2D c) {
      int l = size();
      if (c.rows() != l)
         throw new IllegalArgumentException
            ("Invalid number of rows in correlation matrix");
      if (c.columns() != l)
         throw new IllegalArgumentException
            ("Invalid number of columns in correlation matrix");
      for (int i1 = 0; i1 < l; i1++)
         c.setQuick (i1, i1, 1.0);
      for (int i1 = 0; i1 < l - 1; i1++)
         for (int i2 = i1 + 1; i2 < l; i2++) {
            double cor = correlation (i1, i2);
            c.setQuick (i1, i2, cor);
            c.setQuick (i2, i1, cor);
         }
   }

   /**
    * Clones this object. This makes a shallow copy of this list, i.e.,
    * this does not clone all the tallies in the list. The created clone
    * is modifiable, even if the original list is unmodifiable.
    */
   public ListOfTallies<E> clone() {
      return (ListOfTallies<E>)super.clone();
   }
}