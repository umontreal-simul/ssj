/*
 * Class:        ListOfTalliesWithCV
 * Description:  List of tallies with control variables
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Eric Buist
 * @since        August 2007

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
package umontreal.ssj.stat.list.lincv;

import cern.colt.list.DoubleArrayList;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

import umontreal.ssj.probdist.StudentDist;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.TallyStore;
import umontreal.ssj.stat.list.ListOfTallies;
import umontreal.ssj.stat.list.ListOfTalliesWithCovariance;
import umontreal.ssj.util.PrintfFormat;

/**
 * Represents a list of tallies with control variables that inherits the
 * functionalities of a list of tallies, and accepts vectors of length
 * @f$p+q@f$. The first @f$p@f$ tallies in the list correspond to
 * @f$\mathbf{X}@f$ while the @f$q@f$ last tallies correspond to
 * @f$\mathbf{C}@f$. Methods are provided for adding observations with
 * control variables, and estimating the @f$\boldsymbol{\beta}^*@f$ matrix
 * from the sample covariances. Convenience methods are also provided for the
 * common cases where @f$q=1@f$ and @f$p=1@f$.
 *
 * During simulation or pilot runs, one uses an `add` method to add
 * observations with the associated values of the control variables. After
 * the vectors of observations (and controls) are collected, one can obtain
 * an estimate of @f$\boldsymbol{\beta}^*@f$, by using  #estimateBeta. One
 * can also set @f$\boldsymbol{\beta}@f$ to an arbitrary matrix.
 *
 * Before obtaining results with control variables, it is important to set
 * the expected values of the controls, by using
 * #setExpectedValues(double[]). After @f$E[\mathbf{C}]@f$ is set, methods
 * are then available to get the controlled average, and the controlled
 * covariance matrix. One can also compute a confidence interval on a
 * component @f$X_i@f$ of @f$\mathbf{X}@f$ taking the control variables into
 * account.
 *
 * The @f$\boldsymbol{\beta}^*@f$ matrix can be estimated using pilot runs.
 * In this context, if the list of tallies is reset, the values of
 * @f$\boldsymbol{\beta}@f$ and @f$E[\mathbf{C}]@f$ are retained.
 * Observations generated from production runs can then be added to the list
 * of tallies, and the controlled average can be computed.
 *
 * The following pseudocode illustrates how @f$q@f$ control variables can be
 * used with a @f$p@f$-dimensional vector, the value of @f$E[\mathbf{C}]@f$
 * being stored in `expControls`.
 *
 * @code
 *
 *    ListOfTalliesWithCV<Tally> list = ListOfTalliesWithCV.createWithTally
 * (p, q);
 *    list.setExpectedValues (expControls);
 *    // Choose between NO PILOT RUNS or PILOT RUNS
 *
 *    NO PILOT RUNS {
 *       // Make production runs and use list.add to add (X, C) vectors.
 *       list.estimateBeta();  // Biased consistent estimator of the beta
 * matrix
 *    }
 *
 *    PILOT RUNS {
 *       // Make pilot runs and use list.add to add (X, C) vectors.
 *       list.estimateBeta();   // Unbiased estimator of the beta matrix
 *       list.init();
 *       // Make production runs and use list.add to add (X, C) vectors.
 *    }
 *
 *    // Call list.averageWithCV, list.covarianceWithCV, etc.
 *
 * @endcode
 *
 * <div class="SSJ-bigskip"></div>
 */
public class ListOfTalliesWithCV<E extends Tally>
             extends ListOfTalliesWithCovariance<E> {
   private static Algebra alg = new Algebra();
   private DoubleMatrix2D beta;
   private double[] exp;
   private double[] tmp;
   private int q;
   
   private DoubleMatrix2D tempPP;
   private DoubleMatrix2D tempQQ;
   private DoubleMatrix2D tempPQ;
   private DoubleMatrix2D tempQP;

   /**
    * Constructs a new empty list of tallies with no control variable.
    * After calling this constructor, one must fill the list with tallies,
    * call `setNumControlVariables`, and `init`.
    */
   public ListOfTalliesWithCV() {
      super();
   }

   /**
    * Constructs a new empty list of tallies with no control variable and
    * name `name`. After calling this constructor, one must fill the list
    * with tallies, call `setNumControlVariables`, and `init`.
    *  @param name         the name of the list.
    */
   public ListOfTalliesWithCV (String name) {
      super (name);
   }

   /**
    * This factory method constructs and returns a list of tallies with
    * `p+q` new instances of  @ref umontreal.ssj.stat.Tally, `q` being the
    * number of control variables.
    *  @param p            the value of @f$p@f$.
    *  @param q            the value of @f$q@f$.
    *  @return the created list.
    */
   public static ListOfTalliesWithCV<Tally> createWithTally (int p, int q) {
      ListOfTalliesWithCV<Tally> list = new ListOfTalliesWithCV<Tally>();
      int size = p + q;
      for (int i = 0; i < size; i++)
         list.add (new Tally());
      list.setNumControlVariables (q);
      list.init();
      return list;
   }

   /**
    * This factory method constructs and returns a list of tallies with
    * `p+q` new instances of  @ref umontreal.ssj.stat.TallyStore, `q`
    * being the number of control variables.
    *  @param p            the value of @f$p@f$.
    *  @param q            the value of @f$q@f$.
    *  @return the created list.
    */
   public static ListOfTalliesWithCV<TallyStore> createWithTallyStore
                                                    (int p, int q) {
      ListOfTalliesWithCV<TallyStore> list = new ListOfTalliesWithCV<TallyStore>();
      int size = p + q;
      for (int i = 0; i < size; i++)
         list.add (new TallyStore());
      list.setNumControlVariables (q);
      list.init();
      return list;
   }


   public void init() {
      super.init();
      internalInit();
   }

   private void internalInit() {
      int p = sizeWithoutCV();
      beta = new DenseDoubleMatrix2D (q, p);
      exp = new double[q];
      tmp = new double[p + q];
      
      tempPP = new DenseDoubleMatrix2D (p, p);
      tempQQ = new DenseDoubleMatrix2D (q, q);
      tempPQ = new DenseDoubleMatrix2D (p, q);
      tempQP = new DenseDoubleMatrix2D (q, p);
      setUnmodifiable();
   }

/**
 * Returns the current matrix @f$\boldsymbol{\beta}@f$. By default, this
 * matrix is filled with 0’s.
 *  @return the @f$\boldsymbol{\beta}@f$ matrix.
 */
public DoubleMatrix2D getBeta() {
      return beta;
   }

   /**
    * Sets the @f$\boldsymbol{\beta}@f$ matrix to `beta`. The given matrix
    * must be @f$q\times p@f$, or an exception is thrown.
    *  @param beta         the new @f$\boldsymbol{\beta}@f$ matrix.
    */
   public void setBeta (DoubleMatrix2D beta) {
      if (beta.rows() != getNumControlVariables())
         throw new IllegalArgumentException (
               "The number of rows in beta must be equal to q");
      if (beta.columns() != sizeWithoutCV())
         throw new IllegalArgumentException (
               "The number of columns in beta must be equal to p");
      this.beta = beta;
   }

   /**
    * Gets the expected value of the @f$i@f$th control variable. By
    * default, this is set to 0.
    *  @param i            the index of the control variable.
    *  @return the queried expectation.
    */
   public double getExpectedValue (int i) {
      return exp[i];
   }

   /**
    * Sets the expected value of the @f$i@f$th control variable to `e`.
    *  @param i            the index of the control variable.
    *  @param e            the new value of the expectation.
    */
   public void setExpectedValue (int i, double e) {
      exp[i] = e;
   }

   /**
    * Returns @f$E[\mathbf{C}]@f$, the expected value of the vector of
    * control variables.
    *  @return the expected values of the control variables.
    */
   public double[] getExpectedValues() {
      return exp;
   }

   /**
    * Sets @f$E[\mathbf{C}]@f$ to `exp`. The length of the given array
    * must be @f$q@f$, or an exception is thrown.
    *  @param exp          the new expected values for control variables.
    */
   public void setExpectedValues (double[] exp) {
      if (exp.length != this.exp.length)
         throw new IllegalArgumentException ("Invalid length of exp");
      this.exp = exp;
   }

   /**
    * Returns the size of this list excluding the control variables. This
    * corresponds to the number of output variables @f$p@f$.
    *  @return the size of this list excluding control variables.
    */
   public int sizeWithoutCV() {
      return size() - q;
   }

   /**
    * Returns the number @f$q@f$ of control variables.
    *  @return the number of control variables.
    */
   public int getNumControlVariables() {
      return q;
   }

   /**
    * Sets the number of control variables to `q`. The new value of `q`
    * must not be negative or exceed `size()`. This method throws an
    * exception if it is called after  #init.
    *  @param q            the new number of control variables.
    */
   public void setNumControlVariables (int q) {
      if (beta != null)
         throw new IllegalArgumentException ("Cannot change the number of control variables");
      if (q < 0 || q >= size())
         throw new IllegalArgumentException ("q is negative or greater than or equal to " + size());
      this.q = q;
   }

   /**
    * Fills `c` with the sample correlation matrix of @f$\mathbf{X}@f$.
    *  @param c            the @f$p\times p@f$ matrix filled with
    *                      correlations.
    */
   public void correlationX (DoubleMatrix2D c) {
      int l = sizeWithoutCV();
      if (c.rows() != l)
         throw new IllegalArgumentException (
               "Invalid number of rows in covariance matrix");
      if (c.columns() != l)
         throw new IllegalArgumentException (
               "Invalid number of columns in covariance matrix");
      for (int i1 = 0; i1 < l; i1++)
         for (int i2 = 0; i2 < l; i2++)
            c.setQuick (i1, i2, correlation (i1, i2));
   }

   /**
    * Fills `c` with the sample covariance matrix of @f$\mathbf{X}@f$.
    *  @param c            the @f$p\times p@f$ matrix filled with
    *                      covariances.
    */
   public void covarianceX (DoubleMatrix2D c) {
      int l = sizeWithoutCV();
      if (c.rows() != l)
         throw new IllegalArgumentException (
               "Invalid number of rows in covariance matrix");
      if (c.columns() != l)
         throw new IllegalArgumentException (
               "Invalid number of columns in covariance matrix");
      for (int i1 = 0; i1 < l; i1++)
         for (int i2 = 0; i2 < l; i2++)
            c.setQuick (i1, i2, covariance (i1, i2));
   }

   /**
    * Fills `c` with the sample correlation matrix of @f$\mathbf{C}@f$.
    *  @param c            the @f$q\times q@f$ matrix filled with
    *                      correlations.
    */
   public void correlationC (DoubleMatrix2D c) {
      int p = sizeWithoutCV();
      int q = getNumControlVariables();
      if (c.rows() != q)
         throw new IllegalArgumentException (
               "Invalid number of rows in covariance matrix");
      if (c.columns() != q)
         throw new IllegalArgumentException (
               "Invalid number of columns in covariance matrix");
      for (int i1 = 0; i1 < q; i1++)
         for (int i2 = 0; i2 < q; i2++)
            c.setQuick (i1, i2, correlation (p + i1, p + i2));
   }

   /**
    * Fills `c` with the sample covariance matrix of @f$\mathbf{C}@f$.
    *  @param c            the @f$q\times q@f$ matrix filled with
    *                      covariances.
    */
   public void covarianceC (DoubleMatrix2D c) {
      int p = sizeWithoutCV();
      int q = getNumControlVariables();
      if (c.rows() != q)
         throw new IllegalArgumentException (
               "Invalid number of rows in covariance matrix");
      if (c.columns() != q)
         throw new IllegalArgumentException (
               "Invalid number of columns in covariance matrix");
      for (int i1 = 0; i1 < q; i1++)
         for (int i2 = 0; i2 < q; i2++)
            c.setQuick (i1, i2, covariance (p + i1, p + i2));
   }

   /**
    * Fills `c` with the sample correlation matrix of @f$\mathbf{C}@f$ and
    * @f$\mathbf{X}@f$.
    *  @param c            the @f$q\times p@f$ matrix filled with
    *                      correlations.
    */
   public void correlationCX (DoubleMatrix2D c) {
      int p = sizeWithoutCV();
      int q = getNumControlVariables();
      if (c.rows() != q)
         throw new IllegalArgumentException (
               "Invalid number of rows in covariance matrix");
      if (c.columns() != p)
         throw new IllegalArgumentException (
               "Invalid number of columns in covariance matrix");
      for (int i1 = 0; i1 < q; i1++)
         for (int i2 = 0; i2 < p; i2++)
            c.setQuick (i1, i2, correlation (p + i1, i2));
   }

   /**
    * Fills `c` with the sample covariance matrix of @f$\mathbf{C}@f$ and
    * @f$\mathbf{X}@f$.
    *  @param c            the @f$q\times p@f$ matrix filled with
    *                      covariances.
    */
   public void covarianceCX (DoubleMatrix2D c) {
      int p = sizeWithoutCV();
      int q = getNumControlVariables();
      if (c.rows() != q)
         throw new IllegalArgumentException (
               "Invalid number of rows in covariance matrix");
      if (c.columns() != p)
         throw new IllegalArgumentException (
               "Invalid number of columns in covariance matrix");
      for (int i1 = 0; i1 < q; i1++)
         for (int i2 = 0; i2 < p; i2++)
            c.setQuick (i1, i2, covariance (p + i1, i2));
   }

   /**
    * Adds a new observation @f$(\mathbf{X}, \mathbf{C})@f$ to this list
    * of tallies. The array `x` contains the value for @f$\mathbf{X}@f$,
    * while `c` contains the value of @f$\mathbf{C}@f$.
    *  @param x            the value of @f$\mathbf{X}@f$.
    *  @param c            the value of @f$\mathbf{C}@f$.
    */
   public void add (double[] x, double[] c) {
      if (x.length != sizeWithoutCV())
         throw new IllegalArgumentException
         ("Invalid length of x");
      if (c.length != getNumControlVariables())
         throw new IllegalArgumentException
         ("Invalid length of c");
      System.arraycopy (x, 0, tmp, 0, x.length);
      System.arraycopy (c, 0, tmp, x.length, c.length);
      add (tmp);
   }

   /**
    * Variant of the  #add(double[],double[]) method that can be used when
    * there is only one output variable.
    *  @param x            the output variable.
    *  @param c            the vector of control variables.
    */
   public void add (double x, double[] c) {
      if (sizeWithoutCV() != 1)
         throw new IllegalArgumentException
         ("Cannot use this method if p != 1");
      if (c.length != getNumControlVariables())
         throw new IllegalArgumentException
         ("Invalid length of c");
      tmp[0] = x;
      System.arraycopy (c, 0, tmp, 1, c.length);
      add (tmp);
   }

   /**
    * Variant of the  #add(double[],double[]) that can be used when
    * @f$p=q=1@f$.
    *  @param x            the output variable.
    *  @param c            the control variable.
    */
   public void add (double x, double c) {
      if (sizeWithoutCV() != 1)
         throw new IllegalArgumentException
         ("Cannot use this method if p != 1");
      if (getNumControlVariables() != 1)
         throw new IllegalArgumentException
         ("Cannot use this method if q != 1");
      tmp[0] = x;
      tmp[1] = c;
      add (tmp);
   }

   /**
    * Returns the average of the @f$i@f$th component of
    * @f$\mathbf{X}_{\mathrm{C}}@f$, denoted @f${X}_{\mathrm{C},i}@f$.
    * This corresponds to
    * @f[
    *   {X}_{\mathrm{C},i}=X_i - (\boldsymbol{\beta}_{\cdot,i})^{\mathsf{t}}(\mathbf{C}- E[\mathbf{C}])
    * @f]
    * where @f$\boldsymbol{\beta}_{\cdot,i}@f$ is the @f$i@f$th column of
    * the @f$\boldsymbol{\beta}@f$ matrix.
    *  @param i            the index of the output variable.
    *  @return the controlled average.
    */
   public double averageWithCV (int i) {
      int p = sizeWithoutCV();
      if (i >= p)
         throw new ArrayIndexOutOfBoundsException (i);
      Tally tally = get (i);
      if (tally == null || tally.numberObs() == 0)
         return Double.NaN;
      else {
         double avg = tally.average();
         for (int j = 0; j < q; j++)
            avg -= beta.getQuick (j, i)
                  * (get (p + j).average() - exp[j]);
         return avg;
      }
   }

   /**
    * Computes the sample covariance of @f$\mathbf{X}_{\mathrm{C}}@f$ by
    * replacing @f$\boldsymbol{\Sigma}_{\mathrm{X}}@f$,
    * @f$\boldsymbol{\Sigma}_{\mathrm{C}}@f$, and
    * @f$\boldsymbol{\Sigma}_{\mathrm{CX}}@f$ with the corresponding
    * matrices of empirical covariances. These are obtained by calling
    * #covarianceX,  #covarianceC, and  #covarianceCX, respectively. This
    * gives an estimate of
    * @f[
    *   \mathrm{Cov}(\mathbf{X}_{\mathrm{C}}) = \boldsymbol{\Sigma}_{\mathrm{X}}+ \boldsymbol{\beta}^{\mathsf{t}}\boldsymbol{\Sigma}_{\mathrm{C}}\boldsymbol{\beta}- 2\boldsymbol{\beta}^{\mathsf{t}}\boldsymbol{\Sigma}_{\mathrm{CX}}.
    * @f]
    * @param covCV        the @f$p\times p@f$ matrix that will be filled
    *                      with covariances.
    */
   public void covarianceWithCV (DoubleMatrix2D covCV) {
      final int p = sizeWithoutCV();
      final int q = getNumControlVariables();
      if (covCV.rows() != p || covCV.columns() != q)
         throw new IllegalArgumentException ("Invalid dimensions of covCV");
      DoubleMatrix2D covX = covCV;
      covarianceX (covX);
      DoubleMatrix2D covC = tempQQ;
      covarianceC (covC);
      DoubleMatrix2D covCX = tempQP;
      covarianceCX (covCX);
      
      // covCV contains covX
      // Add beta^t*sigmaC*beta
      beta.viewDice().zMult (covC, tempPQ).zMult (beta, covCV, 1, 1, false, false);
      // Subtract 2beta^t*sigmaCX
      beta.viewDice().zMult (covCX, covCV, -2, 1, false, false);
   }

   /**
    * Computes the covariance between component `i` and `j` of
    * @f$\mathbf{X}_{\mathrm{C}}@f$. This is given by
    * @f[
    *   \mathrm{Cov}({X}_{\mathrm{C},i}, {X}_{\mathrm{C},j}) =\mathrm{Cov}(X_i, X_j) + (\boldsymbol{\beta}_{\cdot, i})^{\mathsf{t}}\boldsymbol{\Sigma}_{\mathrm{C}}\boldsymbol{\beta}_{\cdot, j} - (\boldsymbol{\Sigma}_{\mathrm{CX}, \cdot, i})^{\mathsf{t}}\boldsymbol{\beta}_{\cdot, j} - (\boldsymbol{\Sigma}_{\mathrm{CX}, \cdot, j})^{\mathsf{t}}\boldsymbol{\beta}_{\cdot, i}.
    * @f]
    * @param i            the index of the first component.
    *  @param j            the index of the second component.
    *  @return the covariance.
    */
   public double covarianceWithCV (int i, int j) {
      final int p = sizeWithoutCV();
      //final int q = getNumControlVariables();
      if (i >= p || j >= p)
         throw new IllegalArgumentException ("i >= p or j >= p");
      double cov = covariance (i, j);
      DoubleMatrix2D covC = tempQQ;
      covarianceC (covC);
      DoubleMatrix2D covCX = tempQP;
      covarianceCX (covCX);
      // Add beta_.i^t*sigmaC*beta.j
      cov += beta.viewColumn (i).zDotProduct (covC.zMult (beta.viewColumn (j), null));
      // Subtract sigmaCX_.i*beta_.j
      cov -= covCX.viewColumn (i).zDotProduct (beta.viewColumn (j));
      // Subtract sigmaCX_.j*beta_.i
      cov -= covCX.viewColumn (j).zDotProduct (beta.viewColumn (i));
      return cov;
   }

   /**
    * Fills the given array with the controlled averages.
    */
   public void averageWithCV (double[] a) {
      int p = sizeWithoutCV();
      if (a.length != p)
         throw new IllegalArgumentException (
               "Invalid length of the given array: given length is " + a.length
                     + ", required length is " + p);
      for (int i = 0; i < p; i++)
         a[i] = averageWithCV (i);
   }

   /**
    * Fills the given array with the averages without control variables.
    *  @param a            the array to be filled with averages.
    */
   public void averageX (double[] a) {
      int l = sizeWithoutCV();
      if (a.length != l)
         throw new IllegalArgumentException (
               "Invalid length of the given array: given length is " + a.length
                     + ", required length is " + l);
      for (int i = 0; i < a.length; i++) {
         Tally tally = get (i);
         a[i] = tally == null ? Double.NaN : tally.average();
         if (tally.numberObs() == 0)
            a[i] = Double.NaN;
      }
   }

   /**
    * Fills the given array with the averages of the control variables.
    *  @param a            the array to be filled with averages.
    */
   public void averageC (double[] a) {
      int p = sizeWithoutCV();
      int l = getNumControlVariables();
      if (a.length != l)
         throw new IllegalArgumentException (
               "Invalid length of the given array: given length is " + a.length
                     + ", required length is " + l);
      for (int i = 0; i < a.length; i++) {
         Tally tally = get (i + p);
         a[i] = tally == null ? Double.NaN : tally.average();
         if (tally.numberObs() == 0)
            a[i] = Double.NaN;
      }
   }

   /**
    * Fills the given array with the square root of the variance of each
    * component of @f$\mathbf{X}_{\mathrm{C}}@f$.
    */
   public void standardDeviationWithCV (double[] std) {
      final int l = sizeWithoutCV();
      if (l != std.length)
         throw new IllegalArgumentException ("Invalid length of given array");
      DoubleMatrix2D covCV = tempPP;
      covarianceWithCV (covCV);
      for (int i = 0; i < std.length; i++)
         std[i] = Math.sqrt (covCV.getQuick (i, i));
   }

   /**
    * Fills the given array with the variance of each component of
    * @f$\mathbf{X}_{\mathrm{C}}@f$.
    */
   public void varianceWithCV (double[] v) {
      final int l = sizeWithoutCV();
      if (l != v.length)
         throw new IllegalArgumentException ("Invalid length of given array");
      DoubleMatrix2D covCV = tempPP;
      covarianceWithCV (covCV);
      for (int i = 0; i < v.length; i++)
         v[i] = covCV.getQuick (i, i);
   }

   /**
    * Computes a confidence interval for the @f$i@f$th component of
    * @f$\mathbf{X}_{\mathrm{C}}@f$. This is the same as
    * umontreal.ssj.stat.Tally.confidenceIntervalStudent(double,double[])
    * except that the variance with control variables, obtained by
    * `covarianceWithCV (i, i)`, is used instead of the ordinary variance.
    *  @param i            the index of the component.
    *  @param level        the level of confidence of the interval.
    *  @param centerAndRadius the array that will be filled with the
    *                         center and radius of the interval.
    */
   public void confidenceIntervalStudentWithCV (int i, double level,
                                                double[] centerAndRadius) {
      // Must return an array object, cannot return 2 doubles directly
      double t;
      Tally tally = get (i);
      int numObs = tally.numberObs();
      if (numObs < 2) throw new RuntimeException (
          "Calling confidenceIntervalStudent with < 2 Observations");
      centerAndRadius[0] =  averageWithCV (i);
      t = StudentDist.inverseF (numObs - 1, 0.5 * (level + 1.0));
      centerAndRadius[1] = t * Math.sqrt (covarianceWithCV (i, i) / (double)numObs);
   }

   /**
    * Estimates the @f$\boldsymbol{\beta}^*@f$ matrix from the
    * observations currently in this list of tallies. This uses
    * #covarianceC(DoubleMatrix2D) and  #covarianceCX(DoubleMatrix2D) to
    * get estimates of @f$\boldsymbol{\Sigma}_{\mathrm{C}}@f$, and
    * @f$\boldsymbol{\Sigma}_{\mathrm{CX}}@f$. The result of the
    * estimation of
    * @f[
    *   \boldsymbol{\beta}^* = \boldsymbol{\Sigma}_{\mathrm{C}}^{-1}\boldsymbol{\Sigma}_{\mathrm{CX}}
    * @f]
    * is stored in the matrix returned by  #getBeta.
    */
   public void estimateBeta() {
      DoubleMatrix2D covC = tempQQ;
      covarianceC (covC);
      DoubleMatrix2D covCX = tempQP;
      covarianceCX (covCX);
      beta = alg.solve (covC, covCX);
      assert beta.rows() == getNumControlVariables() &&
         beta.columns() == sizeWithoutCV();
   }

/**
 * Fills the given list of tallies with controlled observations. The list
 * must have size @f$p@f$, and each element is
 * @f[
 *   \mathbf{X}_{\mathrm{C},i}=\mathbf{X}_i - \boldsymbol{\beta}^{\mathsf{t}}(\mathbf{C}_i - E[\mathbf{C}]),
 * @f]
 * where @f$\mathbf{X}_{\mathrm{C},i}@f$ is the @f$i@f$th controlled
 * observation, and @f$\mathbf{X}_i@f$ is the @f$i@f$th observation. This
 * method uses the matrix returned by  #getBeta, and the expectations
 * returned by  #getExpectedValues to obtain @f$\boldsymbol{\beta}@f$, and
 * @f$E[\mathbf{C}]@f$. Tallies in this list with control variables must be
 * capable of storing the observations.
 *  @param a            the list of tallies to be modified.
 */
/*   public void addControlledObservations (ListOfTallies<? extends Tally> a) {
      int p = sizeWithoutCV();
      if (a.size() != p)
         throw new IllegalArgumentException ("Invalid length of the given array of tallies");
      if (p == 0)
         return;
      int numObs = numberObs();
      double[] tmp = new double[p];
      for (int obs = 0; obs < numObs; obs++) {
         for (int i = 0; i < p; i++) {
            TallyStore tallyStore = (TallyStore)get (i);
            double v = tallyStore.getArray().getQuick (obs);
            for (int j = 0; j < q; j++) {
               tallyStore = (TallyStore)get (j + p);
               double c = tallyStore.getArray().getQuick (obs);
               v -= beta.getQuick (j, i)*(c - exp[j]);
            }
            tmp[i] = v;
         }
         a.add (tmp);
      }
   }*/

/**
 * This is a variant of  #addControlledObservations(ListOfTallies) for the
 * case where @f$p=1@f$. Instead of adding the controlled observations in a
 * list of tallies, the method adds them into a regular tally.
 *  @param ta           the tally filled with controlled observations.
 */
/*   public void addControlledObservations (Tally ta) {
      int p = sizeWithoutCV();
      int q = getNumControlVariables();
      if (p != 1)
         throw new IllegalArgumentException ("This method can only be called with p=1");
      int numObs = get (0).numberObs();
      ta.init();
      TallyStore tallyStore = (TallyStore)get (0);
      DoubleArrayList array = tallyStore.getArray();
      for (int obs = 0; obs < numObs; obs++) {
         double v = array.getQuick (obs);
         for (int j = 0; j < q; j++) {
            tallyStore = (TallyStore)get (j + p);
            double c = tallyStore.getArray().getQuick (obs);
            v -= beta.getQuick (j, 0)*(c - exp[j]);
         }
         ta.add (v);
      }
   }*/

   /**
    * Clones this object. This clones the list of tallies as well as the
    * data structures holding the sums of products, the
    * @f$\boldsymbol{\beta}@f$ matrix, and the @f$\mathbf{C}@f$ vector,
    * but this does not clone the tallies comprising the list. The created
    * clone is modifiable, even though the original list is unmodifiable.
    */
   public ListOfTalliesWithCV<E> clone() {
      ListOfTalliesWithCV<E> l = (ListOfTalliesWithCV<E>)super.clone();
      l.beta = (DoubleMatrix2D) beta.clone();
      l.exp = exp.clone();
      l.tmp = tmp.clone();
      l.q = q;

      l.tempPP = (DoubleMatrix2D) tempPP.clone();
      l.tempQQ = (DoubleMatrix2D) tempQQ.clone();
      l.tempPQ = (DoubleMatrix2D) tempPQ.clone();
      l.tempQP = (DoubleMatrix2D) tempQP.clone();   

      return l;
   }

}