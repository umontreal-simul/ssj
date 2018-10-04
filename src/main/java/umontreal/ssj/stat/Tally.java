/*
 * Class:        Tally
 * Description:  statistical collector
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001--2018  Pierre L'Ecuyer and Universite de Montreal
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
package umontreal.ssj.stat;
import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.probdist.StudentDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.ChiSquareDist;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * A subclass of  @ref StatProbe. This type of statistical collector takes a
 * sequence of real-valued observations @f$X_1,X_2,X_3,…@f$ and can return
 * the average, the variance, a confidence interval for the theoretical mean,
 * etc. Each call to  #add provides a new observation. When the broadcasting
 * to observers is activated, the method  #add will also pass this new
 * information to its registered observers. This type of collector does not
 * memorize the individual observations, but only their number, sum, sum of
 * squares, maximum, and minimum. The subclass  @ref TallyStore offers a
 * collector that memorizes the observations.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class Tally extends StatProbe implements Cloneable {
   protected int numObs;
   //private double sumSquares;
   private double curAverage;  // The average of the first numObs observations
   private double curSum2;     // The sum (xi - average)^2 of the first numObs
                               // observations.
   private Logger log = Logger.getLogger ("umontreal.ssj.stat");

   private static enum CIType {CI_NONE, CI_NORMAL, CI_STUDENT};

   protected CIType confidenceInterval = CIType.CI_NONE;
   protected double level = 0.95;
   protected int digits = 3;

   /**
    * Constructs a new unnamed `Tally` statistical probe.
    */
   public Tally() {
      super();
      init();
   }

   /**
    * Constructs a new `Tally` statistical probe with name `name`.
    *  @param name         name of the tally
    */
   public Tally (String name) {
      super();
      this.name = name;
      init();
   }


   /**
    * Set the name of this `Tally` to `name`.
    *  @param name         name of the tally
    */
   public void setName (String name) {
      this.name = name;
   }

   public void init() {
       maxValue = Double.NEGATIVE_INFINITY;
       minValue = Double.POSITIVE_INFINITY;
       sumValue = 0.0;
       // sumSquares = 0.0;
       curAverage = 0.0;
       curSum2 = 0.0;
       numObs = 0;
   }

	/**
	 * Gives a new observation `x` to the statistical collector. If broadcasting to observers is
	 * activated for this object, this method also transmits the new information to the registered
	 * observers by invoking the method #notifyListeners.
	 * 
	 * @param x
	 *            observation being added to this Tally object
	 */
	public void add(double x) {
      if (collect) {
         if (x < minValue) minValue = x;
         if (x > maxValue) maxValue = x;
         numObs++;
         // Algorithme dans Knuth ed. 3, p. 232; voir Wikipedia
         // http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#cite_note-1
         double y = x - curAverage;
         curAverage += y / numObs;
         curSum2 += y*(x - curAverage);
         // On pourrait utiliser l'algorithme correcteur de Kahan pour une
         // meilleure precision.
         // (voir http://en.wikipedia.org/wiki/Kahan_summation_algorithm)
      }
      notifyListeners (x);
   }

	/**
	 * Adds the first `number` observations from the array `x` to this probe.
	 */
	public void add(double[] x, int number) {
		if (collect)
			for (int i = 0; i < number; i++)
				add(x[i]);
	}

	/**
	 * Returns the number of observations given to this probe since its last initialization.
	 * 
	 * @return the number of collected observations
	 */
   public int numberObs() {
      return numObs;
   }
   
   @Override
   public double sum() {
      return numObs * curAverage;
   }

   /**
    * Returns the average value of the observations since the last
    * initialization.
    */
   public double average() {
      if (numObs < 1) {
         //System.err.println (
         //    "******* Tally " + name + ":   calling average() with " + numObs +
         //    " Observation");
         log.logp (Level.WARNING, "Tally", "average",
            "Tally " + name + ":   calling average() with " + numObs +
             " observation");
         return Double.NaN;
      }
      return curAverage;
   }

   /**
    * Returns the sample variance of the observations since the last
    * initialization. This returns `Double.NaN` if the tally contains less
    * than two observations.
    *  @return the variance of the observations
    */
   public double variance() {
      // throws NumberObservationException {
      // if (numObs < 2) throw NumberObservationException;
      if (numObs < 2) {
         //System.err.println (
         //    "******* Tally " + name + ":   calling variance() with " + numObs +
         //    " Observation");
         log.logp (Level.WARNING, "Tally", "variance",
            "Tally " + name + ":   calling variance() with " + numObs +
             " observation");
         return Double.NaN;
      }
      return curSum2 / (numObs-1);
   }

   /**
    * Returns the sample standard deviation of the observations since the
    * last initialization. This returns `Double.NaN` if the tally contains
    * less than two observations.
    *  @return the standard deviation of the observations
    */
   public double standardDeviation() {
      return Math.sqrt (variance());
   }

   /**
    * Computes a confidence interval on the mean. Returns, in elements 0
    * and 1 of the array object `centerAndRadius[]`, the center and
    * half-length (radius) of a confidence interval on the true mean of
    * the random variable @f$X@f$, with confidence level `level`, assuming
    * that the @f$n@f$ observations given to this collector are
    * independent and identically distributed (i.i.d.) copies of @f$X@f$,
    * and that @f$n@f$ is large enough for the central limit theorem to
    * hold. This confidence interval is computed based on the statistic
    * @f[
    *   Z = {\bar{X}_n - \mu\over{S_{n,x}/\sqrt{n}}}
    * @f]
    * where @f$n@f$ is the number of observations given to this collector
    * since its last initialization, @f$\bar{X}_n =@f$ `average()` is the
    * average of these observations, @f$S_{n,x} =@f$ `standardDeviation()`
    * is the empirical standard deviation. Under the assumption that the
    * observations of @f$X@f$ are i.i.d. and @f$n@f$ is large, @f$Z@f$ has
    * the standard normal distribution. The confidence interval given by
    * this method is valid *only if* this assumption is approximately
    * verified.
    *  @param level        desired probability that the (random)
    *                      confidence interval covers the true mean (a
    *                      constant)
    *  @param centerAndRadius array of size 2 in which are returned the
    *                         center and radius of the confidence
    *                         interval, respectively
    */
   public void confidenceIntervalNormal (double level,
                                         double[] centerAndRadius) {
      // Must return an array object, cannot return 2 doubles directly
      double z;
      if (numObs < 2) throw new RuntimeException (
          "Tally " + name +
          ": Calling confidenceIntervalStudent with < 2 Observations");
      centerAndRadius[0] =  average();
      z = NormalDist.inverseF01 (0.5 * (level + 1.0));
      centerAndRadius[1] = z * Math.sqrt (variance() / (double)numObs);
   }

   /**
    * Computes a confidence interval on the mean. Returns, in elements 0
    * and 1 of the array object `centerAndRadius[]`, the center and
    * half-length (radius) of a confidence interval on the true mean of
    * the random variable @f$X@f$, with confidence level `level`, assuming
    * that the observations given to this collector are independent and
    * identically distributed (i.i.d.) copies of @f$X@f$, and that @f$X@f$
    * has the normal distribution. This confidence interval is computed
    * based on the statistic
    * @f[
    *   T = {\bar{X}_n - \mu\over{S_{n,x}/\sqrt{n}}}
    * @f]
    * where @f$n@f$ is the number of observations given to this collector
    * since its last initialization, @f$\bar{X}_n =@f$ `average()` is the
    * average of these observations, @f$S_{n,x} =@f$ `standardDeviation()`
    * is the empirical standard deviation. Under the assumption that the
    * observations of @f$X@f$ are i.i.d. and normally distributed, @f$T@f$
    * has the Student distribution with @f$n-1@f$ degrees of freedom. The
    * confidence interval given by this method is valid *only if* this
    * assumption is approximately verified, or if @f$n@f$ is large enough
    * so that @f$\bar{X}_n@f$ is approximately normally distributed.
    *  @param level        desired probability that the (random)
    *                      confidence interval covers the true mean (a
    *                      constant)
    *  @param centerAndRadius array of size 2 in which are returned the
    *                         center and radius of the confidence
    *                         interval, respectively
    */
   public void confidenceIntervalStudent (double level,
                                          double[] centerAndRadius) {
      // Must return an array object, cannot return 2 doubles directly
      double t;
      if (numObs < 2) throw new RuntimeException (
          "Tally " + name +
          ": Calling confidenceIntervalStudent with < 2 Observations");
      centerAndRadius[0] =  average();
      t = StudentDist.inverseF (numObs - 1, 0.5 * (level + 1.0));
      centerAndRadius[1] = t * Math.sqrt (variance() / (double)numObs);
   }

   /**
    * Similar to  #confidenceIntervalNormal. Returns the confidence
    * interval in a formatted string of the form <br><center>"<tt>95%
    * confidence interval for mean (normal): (32.431,
    * 32.487)</tt>",</center> using @f$d@f$ fractional decimal digits.
    *  @param level        desired probability that the (random)
    *                      confidence interval covers the true mean (a
    *                      constant)
    *  @param d            number of fractional decimal digits
    *  @return a confidence interval formatted as a string
    */
   public String formatCINormal (double level, int d) {
      PrintfFormat str = new PrintfFormat();
      double ci[] = new double[2];
      confidenceIntervalNormal (level, ci);
      str.append ("  " + (100*level) + "%");
      str.append (" conf. interval for the mean (normal approx.): (");
      str.append (7 + d, d-1, d, ci[0] - ci[1]).append (',');
      str.append (7 + d, d-1, d, ci[0] + ci[1]).append (" )" + PrintfFormat.NEWLINE);
      return str.toString();
}

   /**
    * Equivalent to `formatCINormal (level, 3)`.
    *  @param level        desired probability that the (random)
    *                      confidence interval covers the true mean (a
    *                      constant)
    *  @return a confidence interval formatted as a string
    */
   public String formatCINormal (double level) {
      return formatCINormal (level, 3);
   }

   /**
    * Similar to  #confidenceIntervalStudent. Returns the confidence
    * interval in a formatted string of the form<br><center>"<tt>95%
    * confidence interval for mean (student): (32.431,
    * 32.487)</tt>",</center> using @f$d@f$ fractional decimal digits.
    *  @param level        desired probability that the (random)
    *                      confidence interval covers the true mean (a
    *                      constant)
    *  @param d            number of fractional decimal digits
    *  @return a confidence interval formatted as a string
    */
   public String formatCIStudent (double level, int d) {
      PrintfFormat str = new PrintfFormat();
      double ci[] = new double[2];
      confidenceIntervalStudent (level, ci);
      str.append ("  " + (100*level) + "%");
      str.append (" conf. interval for the mean (Student approx.): (");
      str.append (7 + d, d, d-1, ci[0] - ci[1]).append (',');
      str.append (7 + d, d, d-1, ci[0] + ci[1]).append (" )" + PrintfFormat.NEWLINE);
      return str.toString();
   }

   /**
    * Equivalent to `formatCIStudent (level, 3)`.
    *  @param level        desired probability that the (random)
    *                      confidence interval covers the true mean (a
    *                      constant)
    *  @return a confidence interval formatted as a string
    */
   public String formatCIStudent (double level) {
      return formatCIStudent (level, 3);
   }

   /**
    * Computes a confidence interval on the variance. Returns, in elements
    * 0 and 1 of array `interval`, the left and right boundaries
    * @f$[I_1,I_2]@f$ of a confidence interval on the true variance
    * @f$\sigma^2@f$ of the random variable @f$X@f$, with confidence
    * level `level`, assuming that the observations given to this
    * collector are independent and identically distributed (i.i.d.)
    * copies of @f$X@f$, and that @f$X@f$ has the normal distribution.
    * This confidence interval is computed based on the statistic @f$
    * \chi^2_{n-1} = (n-1)S^2_n/\sigma^2 @f$ where @f$n@f$ is the number
    * of observations given to this collector since its last
    * initialization, and @f$S^2_n =@f$ `variance()` is the empirical
    * variance of these observations. Under the assumption that the
    * observations of @f$X@f$ are i.i.d. and normally distributed,
    * @f$\chi^2_{n-1}@f$ has the chi-square distribution with @f$n-1@f$
    * degrees of freedom. Given the `level` @f$ = 1 - \alpha@f$, one has
    * @f$P[\chi^2_{n-1} < x_1] = P[\chi^2_{n-1} > x_2] = \alpha/2@f$
    * and @f$[I_1,I_2] = [(n-1)S^2_n/x_2,\; (n-1)S^2_n/x_1]@f$.
    *  @param level        desired probability that the (random)
    *                      confidence interval covers the true mean (a
    *                      constant)
    *  @param interval     array of size 2 in which are returned the left
    *                      and right boundaries of the confidence
    *                      interval, respectively
    */
   public void confidenceIntervalVarianceChi2 (double level,
                                               double[] interval) {
      // Must return an array object, cannot return 2 doubles directly
      if (numObs < 2) throw new RuntimeException (
          "Tally " + name +
          ":   calling confidenceIntervalVarianceChi2 with < 2 observations");
      double w = (numObs - 1)*variance();
      double x2 = ChiSquareDist.inverseF (numObs - 1, 0.5 * (1.0 + level));
      double x1 = ChiSquareDist.inverseF (numObs - 1, 0.5 * (1.0 - level));
      interval[0] = w / x2;
      interval[1] = w / x1;
   }

   /**
    * Similar to  #confidenceIntervalVarianceChi2. Returns the confidence
    * interval in a formatted string of the form <br><center>"<tt>95.0%
    * confidence interval for variance (chi2): ( 510.642, 519.673
    * )</tt>",</center> using @f$d@f$ fractional decimal digits.
    *  @param level        desired probability that the (random)
    *                      confidence interval covers the true variance
    *  @param d            number of fractional decimal digits
    *  @return a confidence interval formatted as a string
    */
   public String formatCIVarianceChi2 (double level, int d) {
      PrintfFormat str = new PrintfFormat();
      double ci[] = new double[2];
      confidenceIntervalVarianceChi2 (level, ci);
      str.append ("  " + (100*level) + "%");
      str.append (" conf. interval for the variance (chi2 approx.): (");
      str.append (7 + d, d, d-1, ci[0]).append (',');
      str.append (7 + d, d, d-1, ci[1]).append (" )" + PrintfFormat.NEWLINE);
      return str.toString();
}

   /**
    * Returns a formatted string that contains a report on this probe.
    *  @return a statistical report formatted as a string
    */
   public String report() {
      return report(level, digits);
   }

   /**
    * Returns a formatted string that contains a report on this probe with
    * a confidence interval level `level` using @f$d@f$ fractional decimal
    * digits.
    *  @param level        desired probability that the confidence
    *                      interval covers the true mean
    *  @param d            number of fractional decimal digits
    *  @return a statistical report formatted as a string
    */
   public String report (double level, int d) {
      PrintfFormat str = new PrintfFormat();
      str.append ("REPORT on Tally stat. collector ==> " + name);
      str.append (PrintfFormat.NEWLINE + "    num. obs.      min          max        average     variance    standard dev." + PrintfFormat.NEWLINE);
      str.append (7 + d, (int)numObs);   str.append (" ");
      str.append (9 + d, d, d-1, (double)minValue);   str.append (" ");
      str.append (9 + d, d, d-1, (double)maxValue);   str.append (" ");
      str.append (9 + d, d, d-1, (double)average());  str.append (" ");
      str.append (9 + d, d, d-1, (double)variance());  str.append (" ");
      str.append (9 + d, d, d-1, standardDeviation());
      str.append (PrintfFormat.NEWLINE);

      switch (confidenceInterval) {
         case CI_NORMAL:
            str.append (formatCINormal (level, d));
            break;
         case CI_STUDENT:
            str.append (formatCIStudent (level, d));
            break;
      }

      return str.toString();
   }


   public String shortReportHeader() {
      PrintfFormat pf = new PrintfFormat();
      if (showNobs)
         pf.append (-8, "num obs.").append ("  ");
      pf.append (-8, "   min").append ("   ");
      pf.append (-8, "   max").append ("   ");
      pf.append (-8, "   average").append ("   ");
      pf.append (-8, "   variance").append ("   ");
      pf.append (-8, "   std. dev.");
      if (confidenceInterval != CIType.CI_NONE)
         pf.append ("   ").append (-12, "conf. int.");

      return pf.toString();
   }

/**
 * Formats and returns a short statistical report for this tally. The
 * returned single-line report contains the minimum value, the maximum value,
 * the average, the variance, and the standard deviation, in that order, separated by three
 * spaces. If the number of observations is shown in the short report, a
 * column containing the number of observations in this tally is added.
 *  @return the string containing the report.
 */
public String shortReport() {
      PrintfFormat pf = new PrintfFormat();
      if (showNobs)
         pf.append (-8, numberObs());
      pf.append (9, 3, 2, min()).append ("   ");
      pf.append (9, 3, 2, max()).append ("   ");
      pf.append (10, 3, 2, average()).append ("   ");
      if (numberObs() >= 2) {
         pf.append (10, 3, 2, variance()).append ("   ");
         pf.append (11, 3, 2, standardDeviation());
      }
      else
         pf.append (21, "---");

      if (confidenceInterval != CIType.CI_NONE) {
         double[] ci = new double[2];
         switch (confidenceInterval) {
         case CI_NORMAL:
            confidenceIntervalNormal (level, ci);
            break;
         case CI_STUDENT:
            confidenceIntervalStudent (level, ci);
            break;
         }
         pf.append ("   ").append ((100*level) + "% (");
         pf.append (9, 3, 2, ci[0] - ci[1]).append (',');
         pf.append (9, 3, 2, ci[0] + ci[1]).append (")");
      }
      return pf.toString();
   }

    /**
     * Returns a formatted string that contains a report on this probe
     * (as in  #report ), followed by a confidence interval (as in
     * #formatCIStudent ), using @f$d@f$ fractional decimal digits.
     *  @param level        desired probability that the (random)
     *                      confidence interval covers the true mean (a
     *                      constant)
     *  @param d            number of fractional decimal digits
     *  @return a statistical report with a confidence interval,
     * formatted as a string
     */
    public String reportAndCIStudent (double level, int d) {
      CIType oldCIType = confidenceInterval;

      try {
         confidenceInterval = CIType.CI_STUDENT;
         return report(level, d);
      } finally {
         confidenceInterval = oldCIType;
      }
  }

   /**
    * Same as  {@link #reportAndCIStudent() reportAndCIStudent(level, 3)}.
    *  @param level        desired probability that the (random)
    *                      confidence interval covers the true mean (a
    *                      constant)
    *  @return a statistical report with a confidence interval, formatted
    * as a string
    */
   public String reportAndCIStudent (double level) {
      return reportAndCIStudent (level, 3);
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
    * central limit theorem, needs to be included in reports formatted by
    * #report and  #shortReport. The confidence interval is formatted
    * using  #formatCINormal.
    */
   public void setConfidenceIntervalNormal() {
      confidenceInterval = CIType.CI_NORMAL;
   }

   /**
    * Indicates that a confidence interval on the true mean, based on the
    * normality assumption, needs to be included in reports formatted by
    * #report and  #shortReport. The confidence interval is formatted
    * using  #formatCIStudent.
    */
   public void setConfidenceIntervalStudent() {
      confidenceInterval = CIType.CI_STUDENT;
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
    * Clones this object.
    */
   public Tally clone() {
      try {
         return (Tally)super.clone();
      } catch (CloneNotSupportedException e) {
         throw new IllegalStateException ("This Tally cannot be cloned");
      }
   }

}