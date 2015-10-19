/*
 * Class:        MatrixOfFunctionOfMultipleMeansTallies
 * Description:  Statistical collectors for functions of multiple means.
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
   import umontreal.ssj.stat.FunctionOfMultipleMeansTally;
   import umontreal.ssj.util.MultivariateFunction;
   import cern.colt.matrix.DoubleMatrix2D;
   import cern.colt.matrix.DoubleMatrix3D;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.list.ListOfFunctionOfMultipleMeansTallies;
import umontreal.ssj.stat.list.ListOfTallies;
import umontreal.ssj.util.PrintfFormat;

/**
 * Represents a matrix of statistical collectors for functions of multiple
 * means. Each element of such a matrix is an instance of
 * @ref umontreal.ssj.stat.FunctionOfMultipleMeansTally, and observations can
 * be added with the  #add(double[][][]) method. This class defines a factory
 * method to construct a matrix of tallies computing the same function and
 * sharing the same dimension.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class MatrixOfFunctionOfMultipleMeansTallies
             <E extends FunctionOfMultipleMeansTally>
              extends MatrixOfStatProbes<E> {
   double[][][] temp;

   /**
    * Constructs a new unnamed matrix of function of multiple means
    * tallies with `numRows` rows, and `numColumns` columns, and filled
    * with `null` references.
    *  @param numRows      the number of rows in the matrix.
    *  @param numColumns   the number of columns in the matrix.
    *  @exception NegativeArraySizeException if `numRows` or `numColumns`
    * are negative.
    */
   public MatrixOfFunctionOfMultipleMeansTallies (int numRows, int numColumns) {
      super (numRows, numColumns);
   }

   /**
    * Constructs a new empty matrix of function of multiple means tallies
    * with name `name`, `numRows` rows, and `numColumns` columns, and
    * filled with `null` references.
    *  @param name         the global name of the matrix.
    *  @param numRows      the number of rows in the matrix.
    *  @param numColumns   the number of columns in the matrix.
    *  @exception NegativeArraySizeException if `numRows` or `numColumns`
    * are negative.
    */
   public MatrixOfFunctionOfMultipleMeansTallies (String name, int numRows,
                                                  int numColumns) {
      super (name, numRows, numColumns);
   }

   /**
    * This factory method constructs and returns a matrix of function of
    * multiple means tallies with `numRows` rows, `numColumns` columns,
    * and filled with instances of<br>
    * @ref umontreal.ssj.stat.FunctionOfMultipleMeansTally. Each tally
    * computes the multivariate function `func`, in `d` dimensions.
    *  @param func         the multivariate function computed by the
    *                      tallies.
    *  @param d            the dimension of the tallies.
    *  @param numRows      the number of rows in the matrix.
    *  @param numColumns   the number of columns in the matrix.
    *  @return the created matrix.
    */
   public static MatrixOfFunctionOfMultipleMeansTallies
                 <FunctionOfMultipleMeansTally>
      create (MultivariateFunction func, int d, int numRows, int numColumns) {
      MatrixOfFunctionOfMultipleMeansTallies<FunctionOfMultipleMeansTally> matrix =
          new MatrixOfFunctionOfMultipleMeansTallies<FunctionOfMultipleMeansTally>
                   (numRows, numColumns);
      for (int r = 0; r < numRows; r++)
         for (int c = 0; c < numColumns; c++)
            matrix.set (r, c, new FunctionOfMultipleMeansTally (func, d));
      return matrix;
   }

   /**
    * For each function of multiple means tally with row index&nbsp;`r`
    * and column index&nbsp;`c`, adds the vector of observations `x[r][c]`
    * if collecting is turned ON.
    *  @param x            the vectors of observations to be added.
    *  @exception IllegalArgumentException if the dimensions of `x`,
    * `x[r]`, or `x[r][c]` are incorrect.
    */
   public void add (double[][][] x) {
      if (x.length != rows())
         throw new IllegalArgumentException
            ("The length of the given array must " +
             "correspond to the number of rows in the matrix");
      if (collect)
         for (int i = 0; i < x.length; i++) {
            if (x[i].length != columns())
               throw new IllegalArgumentException
                  ("The length of the array x[" + i +
                   "] must correspond to the number of columns in the matrix");
            for (int j = 0; j < x[i].length; j++)
               get (i, j).add (x[i][j]);
         }
   }

   /**
    * Equivalent to  {@link #add(double[][][]) add(x.toArray())}, without
    * copying the elements of `x` into a temporary 3D array. Each slice of
    * the 3D matrix `x` corresponds to a row of the matrix of probes. Each
    * row of the 3D matrix corresponds to a column of the matrix of
    * probes. Each column of the 3D matrix corresponds to a component of a
    * vector of observations added to a function of multiple means tally.
    * This can be used only when all the function of multiple means
    * tallies in the array share the same dimension.
    *  @param x            the matrix of observations.
    *  @exception IllegalArgumentException if the dimensions of the matrix
    * are incorrect.
    */
   public void add (DoubleMatrix3D x) {
      if (x.rows() != rows())
         throw new IllegalArgumentException
            ("The number of rows of the given matrix must " +
             "correspond to the number of rows in the matrix of probes");
      if (x.columns() != columns())
         throw new IllegalArgumentException
            ("The number of columns of the given matrix must " +
             "correspond to the number of columns in the matrix of probes");
      if (collect)
         for (int i = 0; i < rows(); i++)
            for (int j = 0; j < x.columns(); j++)
               get (i, j).add
                  (x.viewSlice (i).viewRow (j).toArray());
   }

  /**
   * For each element (<tt>r</tt>, <tt>c</tt>) of this matrix of tallies,
   * adds the vector of observations `x[0].get (r, c)`, …, `x[d-1].get (r,
   * c)`. This method can be used only when all tallies in this matrix
   * share the same dimension. It creates a transposed 3D array compatible
   * with  #add(double[][][]), and calls the latter method with this array.
   * For example, let `mt` be a matrix of tallies with dimensions
   * @f$r\times c@f$ whose <tt>l</tt>&nbsp;@f$=rc@f$ elements correspond to
   * ratios. If a program can generate two matrices `num` and `den` of
   * dimensions @f$r\times c@f$ representing observations for the numerator
   * and denominator of the ratios, respectively, the observations can be
   * added using `mt.add (num, den)` instead of creating an intermediate 3D
   * array of observations.
   *  @param x            the array of matrices of observations.
   *  @exception IllegalArgumentException if the length of `x` does not
   * correspond to the dimension of the tally, or the dimensions of the
   * matrices `x[j]` are not equal for all <tt>j=0,…,d-1</tt>.
   */
  public void addSameDimension (DoubleMatrix2D... x) {
     final int nr = rows();
     final int nc = columns();
     final int d = getDimension();

     if (x.length != d)
        throw new IllegalArgumentException
        ("The length of the given array must be " + d +
              " while its actual length is " + x.length);

     if (x.length == 0 || nr == 0 || nc == 0)
        return;

     if (nr != x[0].rows() || nc != x[0].columns())
        throw new IllegalArgumentException (
              "The given matrices must have the same dimensions");

     for (int i = 0; i < x.length - 1; i++)
        if (x[i].rows() != x[i+1].rows() ||
              x[i].columns() != x[i+1].columns())
           throw new IllegalArgumentException (
           "The given arrays must have the same length");

      if ((temp == null) || (temp.length != nr) || (temp[0].length != nc) || (temp[0][0].length != d))
         temp = new double[nr][nc][d];

      for (int r = 0; r < nr; r++)
         for (int c = 0; c < nc; c++)
            for (int j = 0; j < d; j++)
               temp[r][c][j] = x[j].getQuick (r, c);
     add (temp);
  }

   /**
    * Assuming that each tally in this matrix has the same dimension,
    * returns the dimension of tally (0, 0), or 0 if the matrix has no row
    * or column.
    *  @return the dimension.
    */
   public int getDimension() {
      if (rows () == 0 || columns () == 0)
         return 0;
      FunctionOfMultipleMeansTally t0 = get (0, 0);
      return t0 == null ? 0 : t0.getDimension ();
   }

   /**
    * Assuming that each tally in this matrix contains the same number of
    * observations, returns the number of observations in tally (0, 0), or
    * 0 if the matrix has no row or column.
    *  @return the number of observations.
    */
   public int numberObs() {
      if (rows () == 0 || columns () == 0)
         return 0;
      FunctionOfMultipleMeansTally t0 = get (0, 0);
      return t0 == null ? 0 : t0.numberObs ();
   }

   /**
    * Tests that every tally in this matrix contains the same number of
    * observations. This returns `true` if and only if the matrix has no
    * row and no column, or if all tallies have the same number of
    * observations.
    *  @return the success indicator of the test.
    */
   public boolean areAllNumberObsEqual() {
      final int nr = rows ();
      final int nc = columns ();
      int n = numberObs ();
      for (int r = 0; r < nr; r++)
         for (int c = 0; c < nc; c++) {
            FunctionOfMultipleMeansTally t = get (r, c);
            if (t.numberObs () != n)
               return false;
         }
      return true;
   }

   /**
    * Computes the average for each function of multiple means tally in
    * the matrix. If the tally @f$(r, c)@f$ has no observation, the
    * `Double.NaN` value is stored in the given matrix, at row&nbsp;`r`,
    * column&nbsp;`c`.
    */
   public void average (DoubleMatrix2D m) {
      super.average (m);
      for (int i = 0; i < m.rows(); i++)
         for (int j = 0; j < m.columns(); j++)
            if (!Double.isNaN (m.getQuick (i, j)) &&
                get (i, j).numberObs() == 0)
               m.setQuick (i, j, Double.NaN);
   }

   /**
    * For each tally in the matrix, computes the sample variance, and
    * stores it into the given matrix. If, for some tally @f$(r, c)@f$,
    * there are not enough observations for estimating the variance,
    * `Double.NaN` is stored in the corresponding element of the given
    * matrix `m`.
    *  @param m            the matrix to be filled with sample variances.
    *  @exception NullPointerException if `m` is `null`.
    *  @exception IllegalArgumentException if `m.rows()` does not
    * correspond to  umontreal.ssj.stat.matrix.MatrixOfStatProbes.rows, or
    * `m.columns()` does not correspond to
    * umontreal.ssj.stat.matrix.MatrixOfStatProbes.columns.
    */
   public void variance (DoubleMatrix2D m) {
      if (m.rows() != rows())
         throw new IllegalArgumentException
            ("Invalid number of rows in the given matrix: required " + rows() +
             " but found " + m.rows());
      if (m.columns() != columns())
         throw new IllegalArgumentException
            ("Invalid number of columns in the given matrix: required " + columns() +
             " but found " + m.columns());
      for (int r = 0; r < rows(); r++)
         for (int c = 0; c < columns(); c++) {
            FunctionOfMultipleMeansTally tally = get (r, c);
            m.setQuick (r, c, tally != null && tally.numberObs() >= 2 ?
                   tally.variance() : Double.NaN);
         }
   }

   /**
    * For each tally in the matrix, computes the standard deviation, and
    * stores it into the matrix `m`. This is equivalent to calling
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
      if (m.rows() != rows())
         throw new IllegalArgumentException
            ("Invalid number of rows in the given matrix: required " + rows() +
             " but found " + m.rows());
      if (m.columns() != columns())
         throw new IllegalArgumentException
            ("Invalid number of columns in the given matrix: required " + columns() +
             " but found " + m.columns());
      for (int r = 0; r < rows(); r++)
         for (int c = 0; c < columns(); c++) {
            FunctionOfMultipleMeansTally tally = get (r, c);
            m.setQuick (r, c, tally != null && tally.numberObs() >= 2 ?
                   tally.standardDeviation() : Double.NaN);
         }
   }

   /**
    * Clones this object. Makes a shallow copy of this matrix, i.e., does
    * not clone all the probes in the matrix.
    */
   public MatrixOfFunctionOfMultipleMeansTallies<E> clone() {
      MatrixOfFunctionOfMultipleMeansTallies<E> clone = (MatrixOfFunctionOfMultipleMeansTallies<E>)super.clone();
      if (temp != null) {
         clone.temp = new double[temp.length][temp[0].length][temp[0].length];
         for (int i = 0; i < temp.length; i++)
            for (int j = 0; j < temp[i].length; i++)
               for (int k = 0; k < temp[i][j].length; k++)
                  clone.temp[i][j][k] = temp[i][j][k];
      }

      return clone;
   }
}