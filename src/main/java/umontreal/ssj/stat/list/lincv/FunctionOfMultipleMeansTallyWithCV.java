/*
 * Class:        FunctionOfMultipleMeansTallyWithCV
 * Description:  function of multiple means tally with control variables
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Éric Buist
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

import cern.colt.matrix.DoubleMatrix2D;

import umontreal.ssj.util.MultivariateFunction;
import umontreal.ssj.stat.FunctionOfMultipleMeansTally;
import umontreal.ssj.stat.Tally;

/**
 * Represents a function of multiple means tally for an estimator with linear
 * control variables. This extends the function of multiple means tally to
 * use a function @f$h(\boldsymbol{\mu},
 * \boldsymbol{\nu})=g(\boldsymbol{\mu})-\boldsymbol{\beta}_{\mathrm{f}}^{\mathsf{t}}(\mathbf{C}-
 * \boldsymbol{\nu})@f$, where
 * @f$\boldsymbol{\beta}_{\mathrm{f}}=(\beta_{\mathrm{f},0},…,\beta_{\mathrm{f},q-1})@f$
 * is a @f$q@f$-dimensional vector. One must provide a definition for
 * @f$g(\boldsymbol{\mu})@f$ and its gradient
 * @f$\nabla g(\boldsymbol{\mu})@f$ through an implementation of the
 * @ref umontreal.ssj.util.MultivariateFunction interface.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class FunctionOfMultipleMeansTallyWithCV extends
      FunctionOfMultipleMeansTally {
   private MultivariateFunction funcNoCV;
   private double[] beta;

   /**
    * Creates a new function of multiple means tally for a function
    * `funcNoCV` of `p` variables, and with `q` control variables. The
    * constructed tally requires vectors of `p+q` values as observations.
    *  @param funcNoCV     the function used.
    *  @param p            the number of dimensions.
    *  @param q            the number of control variables.
    */
   public FunctionOfMultipleMeansTallyWithCV (MultivariateFunction funcNoCV,
                                              int p, int q) {
      // Assumes that the superclass' constructor do not evaluate the function
      super (new LinCVFunction (p + q), ListOfTalliesWithCV.createWithTally (p, q));
      this.funcNoCV = funcNoCV;
      beta = new double[q];
      ((LinCVFunction)getFunction()).initFunctionOfMultipleMeansTallyWithCV (this);
   }

   /**
    * Constructs a new function of multiple means tally with control
    * variables from the list of tallies `l`, and the function `funcNoCV`.
    *  @param funcNoCV     the function being computed.
    *  @param l            the list of tallies used.
    */
   public FunctionOfMultipleMeansTallyWithCV (MultivariateFunction funcNoCV,
                                              ListOfTalliesWithCV<Tally> l) {
      // Assumes that the superclass' constructor do not evaluate the function
      super (new LinCVFunction (l.size()), l);
      this.funcNoCV = funcNoCV;
      beta = new double[l.getNumControlVariables()];
      ((LinCVFunction)getFunction()).initFunctionOfMultipleMeansTallyWithCV (this);
   }

   /**
    * Returns the implementation computing the function
    * @f$g(\boldsymbol{\mu})@f$. This differs from  #getFunction which
    * returns @f$h(\boldsymbol{\mu}, \boldsymbol{\nu})@f$.
    *  @return the implementation computing the function
    * @f$g(\boldsymbol{\mu})@f$.
    */
   public MultivariateFunction getFunctionWithoutCV() {
      return funcNoCV;
   }

   /**
    * Returns the number of control variables being used.
    *  @return the number of control variables.
    */
   public int getNumControlVariables() {
      return beta.length;
   }

   /**
    * Returns the dimension of this tally excluding the control variables.
    * This corresponds to the result of
    * umontreal.ssj.stat.FunctionOfMultipleMeansTally.getDimension minus
    * the number of control variables.
    *  @return the dimension without control variables.
    */
   public int getDimensionWithoutCV() {
      return getDimension() - beta.length;
   }

   /**
    * Returns the value of @f$\beta_{\mathrm{f},i}@f$. This is set to 0 by
    * default.
    *  @param i            the index of the control variable.
    *  @return the value of @f$\beta_{\mathrm{f},i}@f$.
    */
   public double getBeta (int i) {
      return beta[i];
   }

   /**
    * Sets the value of @f$\beta_{\mathrm{f},i}@f$.
    *  @param i            the index of the control variable.
    *  @param b            the value of @f$\beta_{\mathrm{f},i}@f$.
    */
   public void setBeta (int i, double b) {
      beta[i] = b;
   }

   /**
    * Returns the @f$\boldsymbol{\beta}_{\mathrm{f}}@f$ vector. By
    * default, this is set to an array of 0’s.
    *  @return the values of @f$\boldsymbol{\beta}_{\mathrm{f}}@f$.
    */
   public double[] getBeta() {
      return beta;
   }

   /**
    * Sets the value of @f$\boldsymbol{\beta}_{\mathrm{f}}@f$ vector to
    * `beta`. The given array must have length @f$q@f$.
    *  @param beta         the new @f$\boldsymbol{\beta}_{\mathrm{f}}@f$
    *                      vector.
    */
   public void setBeta (double[] beta) {
      if (beta.length != this.beta.length)
         throw new IllegalArgumentException
         ("Invalid length of beta");
      this.beta = beta;
   }

   /**
    * Returns the list of tallies with control variables used by this
    * object.
    */
   public ListOfTalliesWithCV<Tally> getListOfTalliesWithCV() {
      return (ListOfTalliesWithCV<Tally>)super.getListOfTallies();
   }

   /**
    * Gets the expected value of the @f$i@f$th component of
    * @f$\boldsymbol{\nu}@f$. This is set to 0 by default, and is
    * equivalent to calling `getExpectedValue` on the internal list of
    * tallies.
    *  @param i            the component of @f$\boldsymbol{\nu}@f$ to
    *                      query.
    *  @return the queried expectation.
    */
   public double getExpectedValue (int i) {
      return getListOfTalliesWithCV().getExpectedValue (i);
   }

   /**
    * Sets the expected value of the @f$i@f$th component of
    * @f$\boldsymbol{\nu}@f$ to `e`. This is equivalent to calling
    * `setExpectedValue` on the internal list of tallies.
    *  @param i            the component of @f$\boldsymbol{\nu}@f$ to set.
    *  @param e            the new value of the expectation.
    */
   public void setExpectedValue (int i, double e) {
      getListOfTalliesWithCV().setExpectedValue (i, e);
   }

   /**
    * Gets an array containing the vector @f$\boldsymbol{\nu}@f$. This
    * returns an array of 0’s by default, and is equivalent to calling
    * `getExpectedValues` on the internal list of tallies.
    *  @return the vector @f$\boldsymbol{\nu}@f$.
    */
   public double[] getExpectedValues() {
      return getListOfTalliesWithCV().getExpectedValues();
   }

   /**
    * Sets the vector @f$\boldsymbol{\nu}@f$ to to given array `exp`. The
    * length of `exp` must be @f$q@f$. This is equivalent to calling
    * `setExpectedValues` on the internal list of tallies.
    *  @param exp          the new @f$\boldsymbol{\nu}@f$ vector.
    */
   public void setExpectedValues (double[] exp) {
      getListOfTalliesWithCV().setExpectedValues (exp);
   }

   /**
    * Uses the sample averages and covariances obtained from the internal
    * list of tallies to estimate the
    * @f$\boldsymbol{\beta}_{\mathrm{f}}@f$ vector minimizing the variance
    * of @f$h(\bar{\mathbf{X}}_n, \mathbf{C})@f$. The asymptotically
    * optimal vector is
    * @f[
    *   \boldsymbol{\beta}_{\mathrm{f}}^* = \boldsymbol{\Sigma}_{\mathrm{C}}^{-1}\boldsymbol{\Sigma}_{\mathrm{CX}}\nabla g(\boldsymbol{\mu}).
    * @f]
    * The method estimates this vector by replacing every quantity by its
    * estimator. It first uses
    * umontreal.ssj.stat.list.lincv.ListOfTalliesWithCV.estimateBeta from
    * the internal list of tallies to get a @f$\boldsymbol{\beta}@f$
    * matrix, and calls  #estimateBetaFromMatrix(DoubleMatrix2D) to get
    * the @f$\boldsymbol{\beta}_{\mathrm{f}}@f$ vector.
    */
   public void estimateBeta() {
      ListOfTalliesWithCV a = getListOfTalliesWithCV();
      a.estimateBeta();
      estimateBetaFromMatrix (a.getBeta());
   }

   /**
    * Multiples the given @f$q\times p@f$ matrix by the gradient
    * @f$\nabla g(\bar{\mathbf{X}}_n)@f$ to get an estimate of the
    * @f$\boldsymbol{\beta}_{\mathrm{f}}^*@f$ vector minimizing the
    * variance. The matrix is usually obtained by
    * umontreal.ssj.stat.list.lincv.ListOfTalliesWithCV.estimateBeta. from
    * the internal list of tallies.
    *  @param mbeta        the @f$\boldsymbol{\beta}@f$ matrix.
    */
   public void estimateBetaFromMatrix (DoubleMatrix2D mbeta) {
      ListOfTalliesWithCV a = getListOfTalliesWithCV();
      int p = getDimension() - beta.length;
      int q = beta.length;
      double[] avg = new double[p];
      double[] gradient = new double[p];
      for (int i = 0; i < p; i++)
         avg[i] = a.get (i).average();
      for (int i = 0; i < p; i++)
         gradient[i] = funcNoCV.evaluateGradient (i, avg);
      for (int i = 0; i < q; i++) {
         beta[i] = 0;
         for (int j = 0; j < p; j++)
            beta[i] += mbeta.getQuick (i, j) * gradient[j];
      }
   }

   /**
    * Clones this object and the function which is stored inside. This
    * clones the internal list of tallies as well as each tally in this
    * list.
    */
   public FunctionOfMultipleMeansTallyWithCV clone() {
      FunctionOfMultipleMeansTallyWithCV mta = (FunctionOfMultipleMeansTallyWithCV) super
            .clone();
      mta.beta = beta.clone();

      LinCVFunction fct = new LinCVFunction (getDimension());
      mta.func = fct;
      fct.initFunctionOfMultipleMeansTallyWithCV (mta);

      return mta;
   }


   private static class LinCVFunction implements MultivariateFunction {
      private FunctionOfMultipleMeansTallyWithCV fcv;
      private int pplusq;
      private double[] tmp;

      public LinCVFunction (int pplusq) {
         this.pplusq = pplusq;
      }
     
      public void initFunctionOfMultipleMeansTallyWithCV
      (FunctionOfMultipleMeansTallyWithCV fcv) {
         this.fcv = fcv;
         tmp = new double[fcv.getListOfTalliesWithCV().sizeWithoutCV()];
      }

       public int getDimension() {
          return pplusq;
       }

      public double evaluate (double... x) {
         MultivariateFunction funcNoCV = fcv.getFunctionWithoutCV();

         if (x.length != getDimension())
            throw new IllegalArgumentException ("x has length " + x.length
                  + ", which differs from the dimension " + getDimension());
         int d = getDimension() - fcv.beta.length;
         System.arraycopy (x, 0, tmp, 0, d);
         double gValue = funcNoCV.evaluate (tmp);
         ListOfTalliesWithCV<?> a = fcv.getListOfTalliesWithCV();
         for (int i = 0; i < fcv.beta.length; i++)
            gValue -= fcv.beta[i] * (x[d + i] - a.getExpectedValue (i));
         return gValue;
      }

      public double evaluateGradient (int i, double... x) {
         MultivariateFunction funcNoCV = fcv.getFunctionWithoutCV();

         if (x.length != getDimension())
            throw new IllegalArgumentException ("x has length " + x.length
               + ", which differs from the dimension " + getDimension());
         int d = getDimension() - fcv.beta.length;
         if (i < d) {
            System.arraycopy (x, 0, tmp, 0, d);
            return funcNoCV.evaluateGradient (i, tmp);
         }
         return -fcv.beta[i - d];
      }
   }
}