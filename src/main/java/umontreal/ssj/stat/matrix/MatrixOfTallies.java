/*
 * Class:        MatrixOfTallies
 * Description:  Matrix of statistical collectors
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Éric Buist
 * @since        2006

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
package umontreal.ssj.stat.matrix;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.TallyStore;
import cern.colt.matrix.DoubleMatrix2D;
import umontreal.ssj.util.PrintfFormat;
import cern.colt.list.DoubleArrayList;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import umontreal.ssj.stat.list.ListOfTallies;

/**
 * Represents a matrix of tally statistical collectors. Each element of such
 * a matrix is an instance of  @ref umontreal.ssj.stat.Tally, and
 * observations can be added with the  #add(DoubleMatrix2D) method. This
 * class defines factory methods to fill a newly-constructed matrix with
 * `Tally` or `TallyStore` instances.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class MatrixOfTallies<E extends Tally> extends MatrixOfStatProbes<E> {

   /**
    * Constructs a new unnamed matrix of tallies with `numRows` rows, and
    * `numColumns` columns, and filled with *null* references.
    *  @param numRows      the number of rows in the matrix.
    *  @param numColumns   the number of columns in the matrix.
    *  @exception NegativeArraySizeException if `numRows` or `numColumns`
    * are negative.
    */
   public MatrixOfTallies (int numRows, int numColumns) {
      super (numRows, numColumns);
   }

   /**
    * Constructs a new matrix of tallies with name `name`, `numRows` rows,
    * and `numColumns` columns, and filled with `null` references.
    *  @param name         the global name of the matrix.
    *  @param numRows      the number of rows in the matrix.
    *  @param numColumns   the number of columns in the matrix.
    *  @exception NegativeArraySizeException if `numRows` or `numColumns`
    * are negative.
    */
   public MatrixOfTallies (String name, int numRows, int numColumns) {
      super (name, numRows, numColumns);
   }

   /**
    * This factory method constructs and returns a matrix of tallies with
    * `numRows` rows, `numColumns` columns, and filled with new instances
    * of  @ref umontreal.ssj.stat.Tally.
    *  @param numRows      the number of rows in the matrix.
    *  @param numColumns   the number of columns in the matrix.
    *  @return the created matrix.
    */
   public static MatrixOfTallies<Tally> createWithTally
                                 (int numRows, int numColumns) {
      MatrixOfTallies<Tally> matrix = new MatrixOfTallies<Tally> (numRows, numColumns);
      for (int r = 0; r < numRows; r++)
         for (int c = 0; c < numColumns; c++)
            matrix.set (r, c, new Tally());
      return matrix;
   }

   /**
    * This factory method constructs and returns a matrix of tallies with
    * `numRows` rows, `numColumns` columns, and filled with new instances
    * of  @ref umontreal.ssj.stat.TallyStore.
    *  @param numRows      the number of rows in the matrix.
    *  @param numColumns   the number of columns in the matrix.
    *  @return the created matrix.
    */
   public static MatrixOfTallies<TallyStore> createWithTallyStore
                                 (int numRows, int numColumns) {
      MatrixOfTallies<TallyStore> matrix = new MatrixOfTallies<TallyStore>
                     (numRows, numColumns);
      for (int r = 0; r < numRows; r++)
         for (int c = 0; c < numColumns; c++)
            matrix.set (r, c, new TallyStore());
      return matrix;
   }

   /**
    * Adds the observation `x.get(r, c)` in the tally whose row is `r` and
    * column is `c`, for <tt>r = 0, …, </tt>
    * umontreal.ssj.stat.matrix.MatrixOfStatProbes.rows ` - 1`, and <tt>c
    * = 0, …, </tt> umontreal.ssj.stat.matrix.MatrixOfStatProbes.columns
    * ` - 1`. No observation is added if the value is `Double.NaN`, or if
    * observations collecting is turned OFF. If broadcasting is ON, the
    * given matrix is notified to all registered observers. The given
    * matrix `x` not being stored by this object, it can be freely used
    * and modified after the call to this method.
    *  @param x            the matrix of observations.
    *  @exception NullPointerException if `x` is `null`.
    *  @exception IllegalArgumentException if the dimensions of `x` do not
    * correspond to the dimensions of this matrix of tallies.
    */
   public void add (DoubleMatrix2D x) {
      int rows = rows ();
      int columns = columns ();
      if (x.rows () != rows || x.columns () != columns)
         throw new IllegalArgumentException (
               "Incompatible matrix dimensions: given " + x.rows () + "x"
                     + x.columns () + ", required " + rows () + "x"
                     + columns ());
      if (collect)
         for (int r = 0; r < rows; r++)
            for (int c = 0; c < columns; c++) {
               double v = x.getQuick (r, c);
               Tally ta = get (r, c);
               if (!Double.isNaN (v) && ta != null)
                  ta.add (v);
            }
      notifyListeners (x);
   }

   /**
    * Same as  #add(DoubleMatrix2D) for a 2D array.
    */
   public void add (double[][] x) {
      int rows = rows ();
      int columns = columns ();
      if (x.length != rows)
         throw new IllegalArgumentException (
               "Incompatible number of rows: given " + x.length + ", required "
                     + rows);
      if (collect)
         for (int r = 0; r < rows; r++) {
            if (x[r].length != columns)
               throw new IllegalArgumentException (
                     "Incompatible number of columns in row " + r + ": given "
                           + x[r].length + ", but required " + columns);
            for (int c = 0; c < columns; c++) {
               double v = x[r][c];
               Tally ta = get (r, c);
               if (!Double.isNaN (v) && ta != null)
                  ta.add (v);
            }
         }
      notifyListeners (new DenseDoubleMatrix2D (x));
   }

   /**
    * Assuming that each tally in this matrix contains the same number of
    * observations, returns the number of observations in tally (0, 0), or
    * 0 if this matrix has 0 row or column.
    *  @return the number of observations.
    */
   public int numberObs () {
      if (rows () == 0 || columns () == 0)
         return 0;
      Tally t0 = get (0, 0);
      return t0 == null ? 0 : t0.numberObs ();
   }

   /**
    * Tests that every tally in this matrix contains the same number of
    * observations. This returns `true` if and only if all tallies have
    * the same number of observations.
    *  @return the success indicator of the test.
    */
   public boolean areAllNumberObsEqual () {
      final int nr = rows ();
      final int nc = columns ();
      int n = numberObs ();
      for (int r = 0; r < nr; r++)
         for (int c = 0; c < nc; c++) {
            Tally t = get (r, c);
            if (t.numberObs () != n)
               return false;
         }
      return true;
   }

   /**
    * Computes the average for each tally in the matrix. If the tally
    * @f$(r, c)@f$ has no observation, the `Double.NaN` value is stored in
    * the given matrix, at row&nbsp;`r`, and column&nbsp;`c`.
    */
   public void average (DoubleMatrix2D m) {
      for (int i = 0; i < m.rows (); i++)
         for (int j = 0; j < m.columns (); j++) {
            Tally ta = get (i, j);
            double v = ta == null || ta.numberObs() == 0 ? Double.NaN : ta.average();
            m.setQuick (i, j, v);
         }
   }

   /**
    * For each tally in the matrix, computes the sample variance, and
    * stores it into the given matrix. If, for some tally @f$(r, c)@f$,
    * there are not enough observations for estimating the variance,
    * `Double.NaN` is stored in the corresponding element of the given
    * matrix `m`.
    *
    *  @param m   the matrix to be filled with sample variances.
    *  @exception NullPointerException if `m` is `null`.
    *  @exception IllegalArgumentException if `m.rows()` does not
    * correspond to  #rows, or `m.columns()` does not correspond to
    *  @param m   the matrix to be filled with sample variances.
    *  @exception NullPointerException if `m` is `null`.
    *  @exception IllegalArgumentException if `m.rows()` does not
    * correspond to  \ref #rows , or `m.columns()` does not correspond to
    * #columns.
    */
   public void variance (DoubleMatrix2D m) {
      if (m.rows () != rows ())
         throw new IllegalArgumentException (
               "Invalid number of rows in the given matrix: required "
                     + rows () + " but found " + m.rows ());
      if (m.columns () != columns ())
         throw new IllegalArgumentException (
               "Invalid number of columns in the given matrix: required "
                     + columns () + " but found " + m.columns ());
      for (int r = 0; r < rows (); r++)
         for (int c = 0; c < columns (); c++) {
            Tally tally = get (r, c);
            m.setQuick (r, c, tally != null && tally.numberObs () >= 2 ? tally
                  .variance () : Double.NaN);
         }
   }

   /**
    * For each tally in the matrix, computes the standard deviation, and
    * stores it into the given matrix. This is equivalent to calling
    * #variance(DoubleMatrix2D) and performing a square root on every
    * element of the filled matrix.
    *  @param m            the matrix to be filled with standard
    *                      deviations.
    *  @exception NullPointerException if `m` is `null`.
    *  @exception IllegalArgumentException if `m.rows()` does not
    * correspond to  umontreal.ssj.stat.matrix.MatrixOfStatProbes.rows, or
    * `m.columns()` does not correspond to
    * umontreal.ssj.stat.matrix.MatrixOfStatProbes.columns.
    */
   public void standardDeviation (DoubleMatrix2D m) {
      if (m.rows () != rows ())
         throw new IllegalArgumentException (
               "Invalid number of rows in the given matrix: required "
                     + rows () + " but found " + m.rows ());
      if (m.columns () != columns ())
         throw new IllegalArgumentException (
               "Invalid number of columns in the given matrix: required "
                     + columns () + " but found " + m.columns ());
      for (int r = 0; r < rows (); r++)
         for (int c = 0; c < columns (); c++) {
            Tally tally = get (r, c);
            m.setQuick (r, c, tally != null && tally.numberObs () >= 2 ? tally
                  .standardDeviation () : Double.NaN);
         }
   }

   /**
    * Clones this object.
    */
   public MatrixOfTallies<E> clone() {
      return (MatrixOfTallies<E>)super.clone ();
   }
}

