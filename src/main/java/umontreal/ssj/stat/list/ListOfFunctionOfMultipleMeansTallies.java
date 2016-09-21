/*
 * Class:        ListOfFunctionOfMultipleMeansTallies
 * Description:  List of statistical collectors for a vector of functions
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Éric Buist 
 * @since        2007

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
package umontreal.ssj.stat.list;

import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.stat.FunctionOfMultipleMeansTally;
import umontreal.ssj.util.MultivariateFunction;
import umontreal.ssj.stat.Tally;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

/**
 * Represents a list of tally statistical collectors for a vector of
 * functions of multiple means. Each element of such a list is an instance of
 * @ref umontreal.ssj.stat.FunctionOfMultipleMeansTally, and observations can
 * be added with the  #add(double[][]) method. This class defines a factory
 * method to construct a list of tallies computing the same function and
 * sharing the same dimension.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class ListOfFunctionOfMultipleMeansTallies<E extends
             FunctionOfMultipleMeansTally> extends ListOfStatProbes<E> {
   double[][] temp = null;

   /**
    * Constructs a new empty list of tallies.
    */
   public ListOfFunctionOfMultipleMeansTallies() {
      super();
   }

   /**
    * Constructs a new empty list of tallies with name `name`.
    *  @param name         the name of this list.
    */
   public ListOfFunctionOfMultipleMeansTallies (String name) {
      super (name);
   }

   /**
    * This factory method constructs and returns a list of tallies with
    * `size` instances of<br>
    * @ref umontreal.ssj.stat.FunctionOfMultipleMeansTally. Each tally
    * computes the multivariate function `func`, in `d` dimensions.
    *  @param func         the multivariate function computed by the
    *                      tallies.
    *  @param d            the dimension of the tallies.
    *  @param size         the size of the list.
    *  @return the created list.
    */
   public static ListOfFunctionOfMultipleMeansTallies<FunctionOfMultipleMeansTally>
      create (MultivariateFunction func, int d, int size) {
      ListOfFunctionOfMultipleMeansTallies<FunctionOfMultipleMeansTally> list =
          new ListOfFunctionOfMultipleMeansTallies<FunctionOfMultipleMeansTally>();
      for (int i = 0; i < size; i++)
         list.add (new FunctionOfMultipleMeansTally (func, d));
      return list;
   }

   /**
    * For each tally `i` in this list, adds the vector `x[i]`. Since each
    * element of the 2D array `x` can have a different dimension, it is
    * not necessary for `x` to be rectangular, although this is generally
    * the case. If collecting is turned ON, `x[i]` is added into the
    * <tt>i</tt>th tally of this list.
    *  @param x            the array of vectors of observations being
    *                      added.
    *  @exception IllegalArgumentException if the length of `x` does not
    * correspond to `size()`, or the length of `x[i]` does not correspond
    * to the dimension of the underlying tally&nbsp;`i`, for at least one
    * `i`.
    */
   public void add (double[][] x) {
      if (x.length != size())
         throw new IllegalArgumentException
            ("Invalid number of vectors of observations: the given length is " +
             x.length + ", but the required length is " + size());
      if (collect)
         for (int i = 0; i < size(); i++) {
            FunctionOfMultipleMeansTally ta = get (i);
            if (ta != null)
               ta.add (x[i]);
         }
   }

   /**
    * Equivalent to  {@link #add(double[][]) add(x.toArray())}, without
    * copying the elements of `x` into a temporary 2D array. This can be
    * used only when all the function of multiple means tallies in this
    * list share the same dimension.
    *  @param x            the matrix of observations being added, each
    *                      row corresponding to a vector added to a tally.
    *  @exception IllegalArgumentException if the number of rows in `x`
    * does not correspond to `size()`, or the number of columns does not
    * correspond to the dimension of the underlying tallies.
    */
   public void add (DoubleMatrix2D x) {
      if (x.rows() != size())
         throw new IllegalArgumentException
            ("Invalid number of vectors of observations: the given number is " +
             x.rows() + ", but the required number is " + size());
      if (collect)
         for (int i = 0; i < size(); i++) {
            FunctionOfMultipleMeansTally ta = get (i);
            if (ta != null)
               ta.add (x.viewRow (i).toArray());
         }
   }

  /**
   * For each element `i` of this list of tallies, adds the vector of
   * observations `x[0][i]`, …, `x[d-1][i]`. This method can be used only
   * when all tallies in this list share the same dimension. It creates a
   * transposed 2D array compatible with  #add(double[][]), and calls the
   * latter method with this array. For example, let `lt` be a list of
   * tallies whose `l` elements correspond to ratios. If a program can
   * generate two arrays `num` and `den` of length `l` representing
   * observations for the numerator and denominator of the ratios,
   * respectively, the observations can be added using `lt.add (num, den)`
   * instead of creating an intermediate matrix of observations.
   *  @param x            the 2D array of observations.
   *  @exception IllegalArgumentException if the length of `x` does not
   * correspond to the dimension of the tally, or the length of the arrays
   * `x[j]` are not equal for all <tt>j=0,…,d-1</tt>.
   */
  public void addSameDimension (double[]... x) {
     final int l = size();
     final int d = getDimension();
     if (x.length != d)
        throw new IllegalArgumentException
        ("The length of the given array must be " + d +
              " while its actual length is " + x.length);

     if (x.length == 0 || l == 0)
        return;

     if (l != x[0].length)
        throw new IllegalArgumentException ("The given arrays must have the same length");

     for (int i = 0; i < x.length - 1; i++)
        if (x[i].length != x[i+1].length)
           throw new IllegalArgumentException ("The given arrays must have the same length");

     if ((temp == null) || (temp.length != l) || (temp[0].length != d))
        temp = new double[l][d];

     for (int i = 0; i < l; i++)
        for (int j = 0; j < d; j++)
           temp[i][j] = x[j][i];
     add (temp);
  }

   /**
    * Equivalent to  #addSameDimension(double[][]) `x.toArray()`, without
    * copying the elements of `x` into a temporary 1D array. This can be
    * used only when all the function of multiple means tallies in this
    * list share the same dimension.
    */
   public void addSameDimension (DoubleMatrix1D... x) {
     final int l = size();
     final int d = getDimension();
     if (x.length != d)
        throw new IllegalArgumentException
        ("The length of the given array must be " + d +
              " while its actual length is " + x.length);

     if (x.length == 0 || l == 0)
        return;

     if (l != x[0].size())
        throw new IllegalArgumentException ("The given arrays must have the same length");

     for (int i = 0; i < x.length - 1; i++)
        if (x[i].size() != x[i+1].size())
           throw new IllegalArgumentException ("The given arrays must have the same length");

     if ((temp == null) || (temp.length != l) || (temp[0].length != d))
        temp = new double[l][d];

     for (int i = 0; i < l; i++)
        for (int j = 0; j < d; j++)
           temp[i][j] = x[j].get(i);
     add (temp);
   }

   /**
    * Assuming that each tally in this list has the same dimension,
    * returns the dimension of tally&nbsp;0, or 0 if this list is empty.
    *  @return the dimension.
    */
   public int getDimension() {
      if (size() == 0)
         return 0;
      FunctionOfMultipleMeansTally t0 = get (0);
      return t0 == null ? 0 : t0.getDimension();
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
      FunctionOfMultipleMeansTally t0 = get (0);
      return t0 == null ? 0 : t0.numberObs();
   }

   /**
    * Tests that every tally in this list contains the same number of
    * observations. This returns `true` if and only if all tallies have
    * the same number of observations or if the list is empty.
    *  @return the success indicator of the test.
    */
   public boolean areAllNumberObsEqual() {
      final int l = size();
      int n = numberObs();
      for (int i = 1; i < l; i++) {
         FunctionOfMultipleMeansTally t = get (i);
         if (t.numberObs() != n)
            return false;
      }
      return true;
   }

   /**
    * Computes the function of averages for each tally in this list. If
    * the tally `i` has no vector of observations, the `Double.NaN` value
    * is stored at index `i` of the array `a`.
    */
   public void average (double[] a) {
      super.average (a);
      for (int i = 0; i < a.length; i++)
         if (!Double.isNaN (a[i]) && get(i).numberObs() == 0)
            a[i] = Double.NaN;
   }

   /**
    * For each tally in this list, computes the sample variance and stores
    * it into `v`. If, for some tally&nbsp;`i`, there are not enough
    * observations for estimating the variance, `Double.NaN` is stored at
    * index `i` of the array `v`.
    *  @param v            the array to be filled with sample variances.
    *  @exception NullPointerException if `v` is `null`.
    *  @exception IllegalArgumentException if `v.length` does not
    * correspond to  umontreal.ssj.stat.list.ListOfStatProbes.size.
    */
   public void variance (double[] v) {
      if (v.length != size())
         throw new IllegalArgumentException
            ("Incompatible array length: the given length is " + v.length +
             " while the required length is " + size());
      for (int i = 0; i < v.length; i++) {
         FunctionOfMultipleMeansTally tally = get (i);
         if (tally == null || tally.numberObs() < 2)
            v[i] =  Double.NaN;
         else
            v[i] = tally.variance();
      }
   }

   /**
    * For each tally in this list, computes the standard deviation, and
    * stores it into `v`. This is equivalent to calling
    * #variance(double[]) and performing a square root on every element of
    * the filled array.
    *  @param v            the array to be filled with standard
    *                      deviations.
    *  @exception NullPointerException if `v` is `null`.
    *  @exception IllegalArgumentException if `v.length` does not
    * correspond to `size()`.
    */
   public void standardDeviation (double[] v) {
      if (v.length != size())
         throw new IllegalArgumentException
            ("Incompatible array length: the given length is " + v.length +
             " while the required length is " + size());
      for (int i = 0; i < v.length; i++) {
         FunctionOfMultipleMeansTally tally = get (i);
         if (tally == null || tally.numberObs() < 2)
            v[i] = Double.NaN;
         else
            v[i] = tally.standardDeviation();
      }
   }

   /**
    * Clones this object. This makes a shallow copy of this list, i.e.,
    * this does not clone all the objects in the list. The created clone
    * is modifiable, even if the original list is unmodifiable.
    */
   public ListOfFunctionOfMultipleMeansTallies<E> clone() {
      ListOfFunctionOfMultipleMeansTallies<E> clone = (ListOfFunctionOfMultipleMeansTallies<E>)super.clone();
      if (temp != null) {
         clone.temp = new double[temp.length][temp[0].length];
         for (int i = 0; i < temp.length; i++)
            for (int j = 0; j < temp[i].length; i++)
               clone.temp[i][j] = temp[i][j];
      }

      return clone;
   }
}