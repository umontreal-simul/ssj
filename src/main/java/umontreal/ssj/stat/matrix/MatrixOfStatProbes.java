/*
 * Class:        MatrixOfStatProbes
 * Description:  Matrix of statistical probes
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

import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.stat.StatProbe;
import cern.colt.matrix.DoubleMatrix2D;
import java.util.List;
import java.util.ArrayList;
import java.util.AbstractList;
import java.util.Iterator;
import java.util.ConcurrentModificationException;
import java.util.NoSuchElementException;
import java.util.Collections;
import java.util.RandomAccess;

/**
 * Represents a matrix of statistical probes that can be managed
 * simultaneously. Each element of this matrix is a
 * @ref umontreal.ssj.stat.StatProbe instance which can be obtained and
 * manipulated. Alternatively, several methods are provided to manipulate all
 * probes in the matrix simultaneously.
 *
 * Each matrix of probes can have a global name describing the contents of
 * its elements, and local names for each cell. For example, a matrix of
 * statistical probes for the waiting times can have the global name
 * "<tt>Waiting times</tt>" while the first cell has local name "<tt>type 1,
 * period 1</tt>".
 *
 * Facilities are provided to fill matrices of sums, averages, etc. obtained
 * from the individual statistical probes.  DoubleMatrix2D is used instead of
 * 2D arrays because it more efficiently stores the values, and it supports
 * computations on the elements.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class MatrixOfStatProbes<E extends StatProbe>
                                implements Cloneable, Iterable<E> {
   private List<MatrixOfObservationListener> listeners = new ArrayList<MatrixOfObservationListener>();
   protected boolean collect = true;
   protected boolean broadcast = false;
   protected String name;
   private E[] probes;
   private int numRows;
   private int numColumns;
   private int modCount = 0;

   private static enum ListType { ROW, COLUMN }

   /**
    * Constructs a new unnamed matrix of statistical probes with `numRows`
    * rows, and `numColumns` columns, and filled with `null` references.
    *  @param numRows      the number of rows in the matrix.
    *  @param numColumns   the number of columns in the matrix.
    *  @exception NegativeArraySizeException if `numRows` or `numColumns`
    * are negative.
    */
   public MatrixOfStatProbes (int numRows, int numColumns) {
      createProbes (numRows, numColumns);
   }

   /**
    * Constructs a new matrix of statistical probes with name `name`,
    * `numRows` rows, and `numColumns` columns, and filled with `null`
    * references.
    *  @param name         the global name of the matrix.
    *  @param numRows      the number of rows in the matrix.
    *  @param numColumns   the number of columns in the matrix.
    *  @exception NegativeArraySizeException if `numRows` or `numColumns`
    * are negative.
    */
   public MatrixOfStatProbes (String name, int numRows, int numColumns) {
      this.name = name;
      createProbes (numRows, numColumns);
   }

   
   // @SuppressWarnings("unchecked")
   private void createProbes (int numRows, int numColumns) {
      if (numRows < 0)
         throw new NegativeArraySizeException
            ("The number of rows must not be negative");
      if (numColumns < 0)
         throw new NegativeArraySizeException
            ("The number of columns must not be negative");
      this.numRows = numRows;
      this.numColumns = numColumns;
      int length = numRows*numColumns;

      probes = (E[])new StatProbe[length];
   }


   /**
    * Returns the global name of this matrix of statistical probes.
    *  @return the global name of the matrix.
    */
   public String getName() {
      return name;
   }

   /**
    * Sets the global name of this matrix to `name`.
    *  @param name         the new global name of the matrix.
    */
   public void setName (String name) {
      this.name = name;
   }

   /**
    * Returns the number of rows in this matrix.
    *  @return the total number of rows.
    */
   public int rows() {
      return numRows;
   }

   /**
    * Returns the number of columns in this matrix.
    *  @return the total number of columns.
    */
   public int columns() {
      return numColumns;
   }

   /**
    * Sets the number of rows of this matrix of statistical probes to
    * `newRows`, adding or removing cells as necessary. If `newRows` is
    * negative, a `NegativeArraySizeException` is thrown. Otherwise, if
    * `newRows` is equal to  #rows, nothing happens. If the number of rows
    * is reduced, the last  #rows ` - newRows` rows of statistical probes
    * are lost. If the number of rows is increased, the new elements of
    * the matrix are set to `null`.
    *  @param newRows      the new number of rows of the matrix.
    *  @exception IllegalArgumentException if an error occurs during
    * construction of statistical probes.
    *  @exception NegativeArraySizeException if `newRows` is negative.
    */
   public void setRows (int newRows) {
      if (newRows < 0)
         throw new NegativeArraySizeException
            ("The given number of rows is negative");
      if (rows() == newRows)
         return;
      E[] newProbes = (E[])new StatProbe[newRows*columns()];
      int m = Math.min (rows(), newRows);
      System.arraycopy (probes, 0, newProbes, 0, m*columns());
      probes = newProbes;
      numRows = newRows;
      ++modCount;
   }

   /**
    * Similar to  #setRows(int), for setting the number of columns.
    *  @param newColumns   the new number of columns of the matrix.
    *  @exception IllegalArgumentException if an error occurs during
    * construction of statistical probes.
    *  @exception NegativeArraySizeException if `newolumns` is negative.
    */
   public void setColumns (int newColumns) {
      if (newColumns < 0)
         throw new IllegalArgumentException
            ("The given number of columns is negative");
      if (columns() == newColumns)
         return;
      //E[] newProbes = (E[])new StatProbe[probeClass, rows()*newColumns];
      E[] newProbes = (E[])new StatProbe[rows()*newColumns];
      int m = Math.min (columns(), newColumns);
      for (int r = 0; r < rows(); r++)
         System.arraycopy (probes, columns()*r, newProbes, newColumns*r, m);
      probes = newProbes;
      numColumns = newColumns;
      ++modCount;
   }

   /**
    * Returns the statistical probe corresponding to the row&nbsp;`r` and
    * column&nbsp;`c`.
    *  @param r            the row to look at.
    *  @param c            the column to look at.
    *  @return the corresponding statistical probe.
    *
    *  @exception ArrayIndexOutOfBoundsException if `r` or `c` are
    * negative, if `r` is greater than or equal to  #rows, or if `c` is
    * greater than or equal to  #columns.
    */
   public E get (int r, int c) {
      if (r < 0 || r >= numRows)
         throw new ArrayIndexOutOfBoundsException
            ("Row index out of bounds: " + r);
      if (c < 0 || c >= numColumns)
         throw new ArrayIndexOutOfBoundsException
            ("Column index out of bounds: " + c);
      return probes[numColumns*r + c];
   }

   /**
    * Sets the statistical probe corresponding to the row&nbsp;`r` and
    * column&nbsp;`c` to `probe`.
    *  @param r            the row to modify.
    *  @param c            the column to modify.
    *  @param probe        t
    * he new probe.
    *  @exception ArrayIndexOutOfBoundsException if `r` or `c` are
    * negative, if `r` is greater than or equal to  #rows, or if `c` is
    * greater than or equal to  #columns.
    */
   public void set (int r, int c, E probe) {
      if (r < 0 || r >= numRows)
         throw new ArrayIndexOutOfBoundsException
            ("Row index out of bounds: " + r);
      if (c < 0 || c >= numColumns)
         throw new ArrayIndexOutOfBoundsException
            ("Column index out of bounds: " + c);
      probes[numColumns*r + c] = probe;
      ++modCount;
   }

   /**
    * Initializes this matrix of statistical probes by calling
    * `StatProbe.init` on each element.
    */
   public void init() {
      int rows = rows();
      int columns = columns();
      for (int r = 0; r < rows; r++)
         for (int c = 0; c < columns; c++)
            get (r, c).init();
   }

   /**
    * For each probe in the matrix, computes the sum by calling
    * umontreal.ssj.stat.StatProbe.sum, and stores it into the given
    * matrix `m`.
    *  @param m            the matrix to be filled with sums.
    *  @exception NullPointerException if `m` is `null`.
    *  @exception IllegalArgumentException if `m.rows()` does not
    * correspond to  #rows, or `m.columns()` does not correspond to
    * #columns.
    */
   public void sum (DoubleMatrix2D m) {
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
            StatProbe probe = get (r, c);
            m.setQuick (r, c, probe == null ? Double.NaN : probe.sum());
         }
   }

   /**
    * For each statistical probe in the matrix, computes the average by
    * calling  umontreal.ssj.stat.StatProbe.average, and stores it into
    * the given matrix `m`.
    *  @param m            the matrix to be filled with averages.
    *  @exception NullPointerException if `m` is `null`.
    *  @exception IllegalArgumentException if `m.rows()` does not
    * correspond to  #rows, or `m.columns()` does not correspond to
    * #columns.
    */
   public void average (DoubleMatrix2D m) {
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
            StatProbe probe = get (r, c);
            m.setQuick (r, c, probe == null ? Double.NaN : probe.average());
         }
   }

   /**
    * Determines if this matrix of statistical probes is collecting
    * values. The default is `true`.
    *  @return the status of statistical collecting.
    */
   public boolean isCollecting() {
      return collect;
   }

   /**
    * Sets the status of the statistical collecting mechanism to `c`. A
    * `true` value turns statistical collecting ON, a `false` value turns
    * it OFF.
    *  @param c            the status of statistical collecting.
    */
   public void setCollecting (boolean c) {
      collect = c;
   }

   /**
    * Determines if this matrix of statistical probes is broadcasting
    * values to registered observers. The default is `false`.
    *  @return the status of broadcasting.
    */
   public boolean isBroadcasting() {
      return broadcast;
   }

   /**
    * Sets the status of the observation broadcasting mechanism to `b`. A
    * `true` value turns broadcasting ON, a `false` value turns it OFF.
    *  @param b            the status of broadcasting.
    */
   public void setBroadcasting (boolean b) {
      broadcast = b;
   }

   /**
    * Adds the observation listener `l` to the list of observers of this
    * matrix of statistical probes.
    *  @param l            the new observation listener.
    *  @exception NullPointerException if `l` is `null`.
    */
   public void addMatrixOfObservationListener (MatrixOfObservationListener l) {
      if (l == null)
         throw new NullPointerException();
      if (!listeners.contains (l))
         listeners.add (l);
   }

   /**
    * Removes the observation listener `l` from the list of observers of
    * this matrix of statistical probes.
    *  @param l            the observation listener to be deleted.
    */
   public void removeMatrixOfObservationListener
                      (MatrixOfObservationListener l) {
      listeners.remove (l);
   }

   /**
    * Removes all observation listeners from the list of observers of this
    * matrix of statistical probes.
    */
   public void clearMatrixOfObservationListeners() {
      listeners.clear();
   }

   /**
    * Notifies the observation `x` to all registered observers if
    * broadcasting is ON. Otherwise, does nothing.
    */
   public void notifyListeners (DoubleMatrix2D x) {
      if (!broadcast)
         return;
      // We could also use the enhanced for loop here, but this is less efficient.
      final int nl = listeners.size();
      for (int i = 0; i < nl; i++)
         listeners.get (i).newMatrixOfObservations (this, x);
   }

   /**
    * Returns a list representing a view on row `r` of this matrix of
    * statistical probe. The returned list cannot be modified, and becomes
    * invalid if the number of rows in this matrix of statistical probes
    * is changed.
    *  @param r            the row to look at.
    *  @return the list of statistical probes on the row.
    */
   public List<E> viewRow (int r) {
      return new MyList<E> (this, ListType.ROW, r);
   }

   /**
    * Returns a list representing a view on column `c` of this matrix of
    * statistical probe. The returned list cannot be modified, and becomes
    * invalid if the number of columns in this matrix of statistical
    * probes is changed.
    *  @param c            the column to look at.
    *  @return the list of statistical probes on the column.
    */
   public List<E> viewColumn (int c) {
      return new MyList<E> (this, ListType.COLUMN, c);
   }

   /**
    * Formats a report for the row&nbsp;`r` of the statistical probe
    * matrix. The returned string is constructed by getting a view of row
    * `r` and using `StatProbe.report` on this list.
    *  @param r            the row being reported.
    *  @return the report formatted as a string.
    */
   public String rowReport (int r) {
      return StatProbe.report (getName(), viewRow (r));
   }

   /**
    * Formats a report for the column&nbsp;`c` of the statistical probe
    * matrix. The returned string is constructed by getting a view of
    * column `c` and using `StatProbe.report` on this list.
    *  @param c            the column being reported.
    *  @return the report formatted as a string.
    */
   public String columnReport (int c) {
      return StatProbe.report (getName(), viewColumn (c));
   }

   /**
    * Clones this object. This makes a shallow copy of this matrix, i.e.,
    * this does not clone all the probes in the matrix.
    */
   public MatrixOfStatProbes<E> clone() {
      MatrixOfStatProbes<E> sm;
      try {
         sm = (MatrixOfStatProbes<E>)super.clone();
      }
      catch (CloneNotSupportedException cne) {
         throw new IllegalStateException ("CloneNotSupportedException for a class implementing Cloneable");
      }
      if (probes != null)
         sm.probes = (E[])probes.clone();
      return sm;
   }


   private class MyIterator implements Iterator<E> {
      private int index = 0;
      private int expectedModCount = modCount;

      public boolean hasNext() {
         if (modCount != expectedModCount)
            throw new ConcurrentModificationException();
         return (index < probes.length - 1);
      }

      public E next() {
         if (!hasNext())
            throw new NoSuchElementException();
         return probes[index++];
      }

      public void remove() {
         throw new UnsupportedOperationException("Can not remove an element in a matrix");
      }
   }

   public Iterator<E> iterator() {
      return new MyIterator();
   }

   private static class MyList<E extends StatProbe> extends AbstractList<E> implements RandomAccess {
      private MatrixOfStatProbes<E> matrix;
      private ListType type;
      private int index;

      public MyList (MatrixOfStatProbes<E> matrix, ListType type, int index) {
         this.matrix = matrix;
         this.type = type;
         this.index = index;
      }

      public E get (int index) {
         if (type == ListType.ROW) {
            return matrix.get (this.index, index);
         }
         else {
            return matrix.get (index, this.index);
         }
      }

      public int size() {
         if (type == ListType.ROW) {
            return matrix.numColumns;
         }
         else {
            return matrix.numRows;
         }
      }
   }

}

