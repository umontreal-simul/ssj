/*
 * Class:        FunctionOfMultipleMeansTally
 * Description:  statistical collector for a function of multiple means
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author
 * @since

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
package umontreal.ssj.stat;
   import umontreal.ssj.util.MultivariateFunction;
   import umontreal.ssj.stat.list.ListOfTalliesWithCovariance;
import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.probdist.NormalDist;
import cern.colt.matrix.DoubleMatrix1D;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Represents a statistical collector for estimating a function of multiple
 * means with a confidence interval based on the delta theorem
 * @cite tSER80a&thinsp;. Let @f$\bar{\mathbf{X}}_n=(\bar{X}_{n, 0}, …,
 * \bar{X}_{n, d-1})@f$ be a random vector computed by averaging vectors of
 * observations:
 * @f[
 *   \bar{\mathbf{X}}_n=(1/n)\sum_{i=0}^{n-1}\mathbf{X}_i.
 * @f]
 * Then let @f$\{\bar{\mathbf{X}}_n, n \ge0\}@f$ be a sequence of vectors
 * converging to a vector @f$\boldsymbol{\mu}@f$ when @f$n\to\infty@f$.
 * Then, if @f$g(\bar{\mathbf{X}}_n)@f$ is a continuous function, it
 * converges to @f$g(\boldsymbol{\mu})@f$ as @f$n\to\infty@f$.
 *
 * This class collects @f$\mathbf{X}@f$ vectors in order to compute
 * @f$g(\bar{\mathbf{X}}_n)@f$, to estimate @f$\nu=g(\boldsymbol{\mu})@f$
 * with a confidence interval. The function @f$g(\boldsymbol{\mu})@f$ as well
 * as its gradient @f$\nabla g(\boldsymbol{\mu})@f$ are defined using an
 * implementation of  @ref umontreal.ssj.util.MultivariateFunction.
 *
 * This class defines the methods  {@link #add() add(double[])} for adding
 * vectors of observations to the tally,  {@link #average() average()} for
 * estimating @f$\nu@f$, and  {@link #confidenceIntervalDelta()
 * confidenceIntervalDelta(double, double[])} for applying the delta theorem
 * to compute a confidence interval on @f$\nu@f$. It uses an internal
 * @ref umontreal.ssj.stat.list.ListOfTalliesWithCovariance instance to
 * manage the tallies and covariance estimation.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class FunctionOfMultipleMeansTally extends StatProbe
                                          implements Cloneable {
   private ListOfTalliesWithCovariance<Tally> ta;
   protected MultivariateFunction func;
   private double[] temp;
   private double[] delta;

   private Logger log = Logger.getLogger ("umontreal.ssj.stat");

   private static enum CIType {CI_NONE, CI_DELTA};

   protected CIType confidenceInterval = CIType.CI_NONE;
   protected double level = 0.95;

   /**
    * Constructs a function of multiple means tally with dimension `d`,
    * and function `func`.
    *  @param func         the function being computed.
    *  @param d            the dimension of the tally.
    *  @exception NegativeArraySizeException if `d` is negative.
    */
   public FunctionOfMultipleMeansTally (MultivariateFunction func, int d) {
      ta = ListOfTalliesWithCovariance.createWithTally (d);
      ta.setUnmodifiable();
      this.func = func;
      internalInit();
   }

   /**
    * Constructs a function of multiple means tally with name `name`,
    * dimension `d`, and function `func`. The given name is also used as a
    * global name for the internal list of tallies.
    *  @param func         the function being computed.
    *  @param name         the name of the tally.
    *  @param d            the dimension of the tally.
    *  @exception NegativeArraySizeException if `d` is negative.
    */
   public FunctionOfMultipleMeansTally (MultivariateFunction func,
                                        String name, int d) {
      ta = ListOfTalliesWithCovariance.createWithTally (d);
      ta.setName (name);
      ta.setUnmodifiable();
      this.func = func;
      this.name = name;
      internalInit();
   }

   /**
    * Constructs a function of multiple means tally using the function
    * `func` and the list of tallies `ta` for observation management and
    * covariance estimation.
    *  @param func         the function being computed.
    *  @param ta           the list of tallies to be used.
    *  @exception NullPointerException if `ta` is `null`.
    */
   public FunctionOfMultipleMeansTally (MultivariateFunction func,
                                        ListOfTalliesWithCovariance<Tally> ta) {
      if (ta == null)
         throw new NullPointerException
            ("The list of tallies cannot be null");
      this.ta = ta;
      this.func = func;
      ta.setUnmodifiable();
      name = ta.getName();
      internalInit();
   }


   public void setName (String name) {
      super.setName (name);
      ta.setName (name);
   }

/**
 * Returns the (unmodifiable) list of tallies internally used by this object.
 * *Note*: for the sample covariances to be consistent, each tally in the
 * returned list must always have the same number of observations. Therefore,
 * it is not recommended to manually add observations to the individual
 * tallies, by using  umontreal.ssj.stat.Tally.add(double).
 *  @return the internal list of tallies.
 */
public ListOfTalliesWithCovariance<Tally> getListOfTallies() {
      return ta;
   }

   /**
    * Returns the function of multiple means used by this tally.
    *  @return the function of multiple means.
    */
   public MultivariateFunction getFunction() {
      return func;
   }

   /**
    * Returns the dimension of this tally, i.e., the size of any vector of
    * observations.
    *  @return the dimension of the tally.
    */
   public int getDimension() {
      return ta.size();
   }

   /**
    * Adds a new vector <tt>x = </tt>@f$(X_0, …, X_{d-1})@f$ of
    * observations to this tally. If statistical collecting is turned ON,
    * this adds the vector to the internal list of tallies. This method
    * takes a variable number of arguments to allow one to pass the
    * components of @f$\mathbf{X}@f$ without creating an array when the
    * dimension is fixed. However, one can also give an array to this
    * method. For example, if the computed function is a ratio, for
    * estimating a ratio of two expectations, @f$d=2@f$ and the call
    * `add(x, y)` is equivalent to <tt>add(new double[] \{x, y\})</tt>. On
    * the other hand, if the function is a sum or a product of @f$d@f$
    * terms, the dimension is not constant in general, e.g., it can depend
    * on an input parameter. In this case, one must pass an array to this
    * method, and could also call the `add` method of the internal list of
    * tallies directly.
    *  @param x            the vector of observations being added.
    *  @exception NullPointerException if `x` is `null`.
    *  @exception IllegalArgumentException if the length of `x` does not
    * correspond to the dimension of this tally.
    */
   public void add (double... x) {
      if (x.length != ta.size())
         throw new IllegalArgumentException
            ("Incompatible length of vectors of observations: given length is " + x.length +
             ", but required length is " + ta.size());
      if (collect)
         ta.add (x);
   }

   /**
    * Returns the number of vectors of observations given to this probe
    * since its last initialization.
    *  @return the number of collected vectors of observations.
    */
   public int numberObs() {
      return ta.numberObs();
   }

   /**
    * Computes @f$\bar{\nu}_n=g(\bar{\mathbf{X}}_n)@f$, an estimate of
    * the function of means @f$\nu@f$. Note that if @f$g(\mathbf{X})@f$
    * is non-linear, the @f$\bar{\nu}_n@f$ estimator is biased but
    * consistent: @f$g(\bar{\mathbf{X}}_n)\to g(\boldsymbol{\mu})@f$ if
    * @f$n\to\infty@f$.
    *  @return the estimation of @f$\nu@f$.
    */
   public double average() {
      if (temp.length != ta.size())
         // The dimension of the internal tally could (rarely) change.
         temp = new double[ta.size()];
      ta.average (temp);
      return func.evaluate (temp);
   }

   /**
    * Estimates @f$n\mathrm{Var}(g(\bar{\mathbf{X}}_n))@f$ where @f$n@f$
    * is the number of vectors of observations given to this collector
    * since the last initialization. Assuming that, as @f$n\to\infty@f$,
    * @f$\sqrt{n}(\bar{\mathbf{X}}_n - \boldsymbol{\mu})@f$ converges to a
    * random vector @f$\mathbf{Y}@f$ with mean @f$\mathbf{0}@f$ and
    * covariance matrix @f$\boldsymbol{\Sigma}@f$
    * (@f$\boldsymbol{\Sigma}@f$ is also the covariance matrix of
    * @f$\mathbf{X}@f$), the delta theorem @cite tSER80a&thinsp; shows
    * that @f$\sqrt{n}(g(\bar{\mathbf{X}}_n) -
    * g(\boldsymbol{\mu}))\Rightarrow(\nabla g(\boldsymbol{\mu}))^{\mathsf{t}}\mathbf{Y}@f$.
    * As a result,
    * @f$n\mathrm{Var}[g(\bar{\mathbf{X}}_n)]=n\mathrm{Var}[g(\bar{\mathbf{X}}_n) -
    * g(\boldsymbol{\mu})]=\mathrm{Var}[\sqrt{n}(g(\bar{\mathbf{X}}_n) -
    * g(\boldsymbol{\mu}))]\to\mathrm{Var}[(\nabla g(\boldsymbol{\mu}))^{\mathsf{t}}\mathbf{Y}]=(\nabla g(\boldsymbol{\mu}))^{\mathsf{t}}\boldsymbol{\Sigma}\nabla g(\boldsymbol{\mu})=\sigma^2@f$.
    * Here, @f$\nabla g(\boldsymbol{\mu})@f$ is the gradient of
    * @f$g(\boldsymbol{\mu})@f$. This method computes
    * @f[
    *   S_n^2 = (\nabla g(\bar{\mathbf{X}}_n))^{\mathsf{t}}\mathbf{S}_n\nabla g(\bar{\mathbf{X}}_n),
    * @f]
    * where @f$\mathbf{S}_n@f$ is the matrix of empirical covariances of
    * @f$\mathbf{X}@f$. When @f$n\to\infty@f$, @f$S_n^2\to\sigma^2@f$
    * if @f$\mathbf{S}_n\to\boldsymbol{\Sigma}@f$, and
    * @f$\nabla g(\bar{\mathbf{X}}_n)\to\nabla g(\boldsymbol{\mu})@f$.
    * Therefore, @f$S_n^2@f$ is a biased but consistent estimator of
    * @f$n\mathrm{Var}(g(\bar{\mathbf{X}}_n))@f$.
    *  @return the estimate of the variance.
    */
   public double variance() {
      if (ta.numberObs() < 2) {
         //System.out.println
         //   ("******* FunctionOfMultipleMeansTally with name " + getName()
         //         + ":   Calling variance() with "
         //    + ta.numberObs() + " observation");
         log.logp (Level.WARNING, "FunctionOfMultipleMeansTally", "variance",
            "FunctionOfMultipleMeansTally " + name + ":   calling variance() with " + ta.numberObs() +
             " observation");
         return Double.NaN;
      }
      if (temp.length != ta.size()) {
         // The dimension of the internal tally could (rarely) change.
         temp = new double[ta.size()];
         delta = new double[ta.size()];
      }
      ta.average (temp);
      for (int i = 0; i < delta.length; i++)
         delta[i] = func.evaluateGradient (i, temp);
      for (int i = 0; i < temp.length; i++) {
         temp[i] = 0;
         for (int j = 0; j < delta.length; j++)
            temp[i] += delta[j]*ta.covariance (j, i);
      }
      double sigma = 0;
      for (int i = 0; i < temp.length; i++)
         sigma += temp[i]*delta[i];
      // sigma could be negative but near 0, e.g., -1e-13, so round it
      // to 0 if it is negative.
      return sigma < 0 ? 0 : sigma;
   }

   /**
    * Returns the square root of  #variance.
    *  @return the standard deviation.
    */
   public double standardDeviation() {
      return Math.sqrt (variance());
   }

   /**
    * Computes a confidence interval with confidence level `level` on
    * @f$\nu=g(\boldsymbol{\mu})@f$, using the delta theorem
    * @cite tSER80a&thinsp;. Returns, in elements 0 and 1 of the array
    * object `centerAndRadius[]`, the center and half-length (radius) of a
    * confidence interval on the true function of means @f$\nu@f$, with
    * confidence level `level`, assuming that the vectors of observations
    * given to this collector are independent and identically distributed
    * (i.i.d.) copies of @f$\mathbf{X}@f$, and that
    * @f$\bar{\mathbf{X}}_n@f$ has the @f$d@f$-dimensional normal
    * distribution, which is approximately true if @f$n@f$ is large
    * enough. With this assumption, as @f$n\to\infty@f$,
    * @f[
    *   \sqrt{n}(\bar{\mathbf{X}}_n - \boldsymbol{\mu})
    * @f]
    * converges to a random vector @f$\mathbf{Y}@f$ following the
    * @f$d@f$-dimensional normal distribution with mean @f$\mathbf{0}@f$
    * and covariance matrix @f$\boldsymbol{\Sigma}@f$. According to the
    * delta theorem,
    * @f[
    *   \sqrt{n}(g(\bar{\mathbf{X}}_n) - g(\boldsymbol{\mu}))\Rightarrow(\nabla g(\boldsymbol{\mu}))^{\mathsf{t}}\mathbf{Y},
    * @f]
    * which follows the normal distribution with mean 0 and variance
    * @f[
    *   \sigma^2=(\nabla g(\boldsymbol{\mu}))^{\mathsf{t}}\boldsymbol{\Sigma}\nabla g(\boldsymbol{\mu}).
    * @f]
    * This variance is estimated by @f$S_n^2@f$, which replaces
    * @f$\boldsymbol{\mu}@f$ by @f$\bar{\mathbf{X}}_n@f$, and
    * @f$\boldsymbol{\Sigma}@f$ by a matrix of sample covariances computed
    * by  @ref umontreal.ssj.stat.list.ListOfTallies.
    *
    * The confidence interval is computed based on the statistic
    * @f[
    *   Z = {\sqrt{n}(g(\bar{\mathbf{X}}_n) - g(\boldsymbol{\mu}))\over{S_n}}
    * @f]
    * where @f$n@f$ is the number of vectors of observations given to this
    * collector since its last initialization, @f$g(\bar{\mathbf{X}}_n)
    * =@f$&nbsp; #average is the average of these observations, and @f$S_n
    * =@f$&nbsp; #standardDeviation is the standard deviation. Under the
    * previous assumptions, @f$Z@f$ has the normal distribution with
    * mean&nbsp;0 and variance&nbsp;1. The confidence interval given by
    * this method is valid *only if* these assumptions are approximately
    * verified.
    *  @param level        desired probability that the (random)
    *                      confidence interval covers the true function of
    *                      means (a constant).
    *  @param centerAndRadius array of size 2 in which are returned the
    *                         center and radius of the confidence
    *                         interval, respectively.
    */
   public void confidenceIntervalDelta (double level,
                                        double[] centerAndRadius) {
      // Must return an array object, cannot return 2 doubles directly.
      centerAndRadius[0] =  average();
      double z = NormalDist.inverseF01 (0.5 * (level + 1.0));
      centerAndRadius[1] = z * Math.sqrt (variance() / (double)numberObs());
   }

   /**
    * Similar to  #confidenceIntervalDelta(double,double[]), but returns
    * the confidence interval in a formatted string of the form "<tt>95%
    * confidence interval for function of means: (32.431,  32.487)</tt>",
    * using @f$d@f$ decimal digits of accuracy.
    *  @param level        desired probability that the confidence
    *                      interval covers the true function of means.
    *  @param d            the number of decimal digits of accuracy.
    *  @return a confidence interval formatted as a string.
    */
   public String formatCIDelta (double level, int d) {
      PrintfFormat str = new PrintfFormat();
      double ci[] = new double[2];
      confidenceIntervalDelta (level, ci);

      str.append ("  " + (100*level) + "%");
      str.append (" confidence interval for function of means: (");
      str.append (7 + d, d, d - 1, ci[0] - ci[1]).append (',');
      str.append (7 + d, d, d - 1, ci[0] + ci[1]).append (" )" +
                  PrintfFormat.NEWLINE);
      return str.toString();
   }

   /**
    * Same as  {@link #formatCIDelta(double,int) formatCIDelta(level, 3)}.
    */
   public String formatCIDelta (double level) {
      return formatCIDelta (level, 3);
   }

   /**
    * Returns a string containing a formatted report on this probe. The
    * report contains the function of averages, the standard deviation of
    * @f$g(\mathbf{X})@f$, as well as the number of vectors of
    * observations.
    *  @return a statistical report formated as a string.
    */
   public String report() {
      PrintfFormat str = new PrintfFormat();
      str.append ("REPORT on Tally stat. collector ==> " + name);
      str.append (PrintfFormat.NEWLINE + "    func. of averages    standard dev.  ");
      str.append ("num. obs." + PrintfFormat.NEWLINE);
      str.append (20, 3, 2, (double)average());
      str.append (13, 3, 2, standardDeviation());
      str.append (13, (int)numberObs()).append(PrintfFormat.NEWLINE);

      if (confidenceInterval == CIType.CI_DELTA) {
         str.append (formatCIDelta (level));
      }

      return str.toString();
   }

   /**
    * Returns a string containing a formatted report on this probe (as in
    * #report ), followed by a confidence interval (as in
    * #formatCIDelta(double,int) ).
    *  @param level        desired probability that the confidence
    *                      interval covers the true mean.
    *  @param d            the number of decimal digits of accuracy.
    *  @return statistical report with a confidence interval, formatted as
    * a string.
    */
   public String reportAndCIDelta (double level, int d) {
      CIType oldCIType = confidenceInterval;
      double oldLevel = this.level;

      try {
         confidenceInterval = CIType.CI_DELTA;
         this.level = level;
         return report();
      }
      finally {
         confidenceInterval = oldCIType;
         this.level = oldLevel;
      }
   }

   /**
    * Same as  {@link #reportAndCIDelta(double,int)
    * reportAndCIDelta(level, 3)}.
    *  @param level        desired probability that the confidence
    *                      interval covers the true mean.
    *  @return statistical report with a confidence interval, formatted as
    * a string.
    */
   public String reportAndCIDelta (double level) {
      return reportAndCIDelta (level, 3);
   }


   public String shortReportHeader() {
      PrintfFormat pf = new PrintfFormat();
      pf.append (-20, " func. of averages").append ("   ");
      pf.append (-9, "std. dev.");
      if (showNobs)
         pf.append ("   ").append (-5, "nobs.");
      if (confidenceInterval != CIType.CI_NONE)
         pf.append ("   ").append (-25, "conf. int.");

      return pf.toString();
   }

/**
 * Formats and returns a short statistical report for this function of
 * multiple means tally. The returned single-line report contains the
 * function of averages @f$g(\bar{\mathbf{X}}_n)@f$, and the standard
 * deviation @f$S_n@f$, in that order, separated by three spaces. If the
 * number of observations is included in short reports, a column containing
 * the number of observations is added.
 *  @return the string containing the report.
 */
public String shortReport() {
      PrintfFormat pf = new PrintfFormat();
      pf.append (20, 3, 2, average()).append ("   ");
      if (numberObs() >= 2)
         pf.append (9, 3, 2, standardDeviation());
      else
         pf.append (9, "---");
      if (showNobs)
         pf.append ("   ").append (5, numberObs());

      if (confidenceInterval == CIType.CI_DELTA) {
         double[] ci = new double[2];
         confidenceIntervalDelta (level, ci);
         pf.append ("   ").append ((100*level) + "% (");
         pf.append (9, 3, 2, ci[0] - ci[1]).append (',');
         pf.append (9, 3, 2, ci[0] + ci[1]).append (")");
      }

      return pf.toString();
   }

   /**
    * Returns the level of confidence for the intervals on the mean
    * displayed in reports. The default confidence level is 0.95.
    *  @return desired probability that the (random) confidence interval
    * covers the true mean (a constant)
    */
   public double getConfidenceLevel() {
      return level;
   }

   /**
    * Sets the level of confidence for the intervals on the mean displayed
    * in reports.
    *  @param level        desired probability that the (random)
    *                      confidence interval covers the true mean (a
    *                      constant)
    */
   public void setConfidenceLevel (double level) {
      if (level < 0.0)
         throw new IllegalArgumentException("level < 0");
      if (level >= 1.0)
         throw new IllegalArgumentException("level >= 1");
      this.level = level;
   }

   /**
    * Indicates that no confidence interval needs to be printed in reports
    * formatted by  #report, and  #shortReport. This restores the default
    * behavior of the reporting system.
    */
   public void setConfidenceIntervalNone() {
      confidenceInterval = CIType.CI_NONE;
   }

   /**
    * Indicates that a confidence interval on the true mean, based on the
    * delta and central limit theorems, needs to be included in reports
    * formatted by  #report, and  #shortReport. The confidence interval is
    * formatted using  #formatCIDelta.
    */
   public void setConfidenceIntervalDelta() {
      confidenceInterval = CIType.CI_DELTA;
   }

   /**
    * Determines if the number of observations must be displayed in
    * reports. By default, the number of observations is displayed.
    *  @param showNumObs   the value of the indicator.
    */
   public void setShowNumberObs (boolean showNumObs) {
      showNobs = showNumObs;
   }

   /**
    * Returns `Double.NaN`.
    */
   public double min() {
      return Double.NaN;
   }

   /**
    * Returns `Double.NaN`.
    */
   public double max() {
      return Double.NaN;
   }

   /**
    * Returns `Double.NaN`.
    */
   public double sum() {
      return Double.NaN;
   }


   private void internalInit() {
      if ((func.getDimension() != ta.size()) && (func.getDimension() != -1))
         throw new IllegalArgumentException (
            "The dimension of the function must be equal to d, or equal to -1");
      temp = new double[ta.size()];
      delta = new double[ta.size()];
   }

   public void init() {
      ta.init();
      internalInit();
   }

/**
 * Clones this object. This clones the internal list of tallies as well as
 * each tally in this list.
 */
public FunctionOfMultipleMeansTally clone() {
      FunctionOfMultipleMeansTally mta;
      try {
         mta = (FunctionOfMultipleMeansTally)super.clone();
      }
      catch (CloneNotSupportedException cne) {
         throw new IllegalStateException
            ("CloneNotSupportedException for a class implementing Cloneable");
      }

      mta.ta = new ListOfTalliesWithCovariance<Tally>();
      for (int i = 0; i < ta.size(); i++)
      {
         mta.ta.add (ta.get(i).clone());
      }

      mta.temp = new double[ta.size()];
      mta.delta = new double[ta.size()];
      return mta;
   }
}