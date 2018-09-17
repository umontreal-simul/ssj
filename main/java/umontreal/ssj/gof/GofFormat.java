/*
 * Class:        GofFormat
 * Description:
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
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
package umontreal.ssj.gof;
   import cern.colt.list.*;
import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.probdist.*;
import java.io.PrintWriter;

/**
 * This class contains methods used to format results of GOF test statistics,
 * or to apply a series of tests simultaneously and format the results. It is
 * in fact a translation from C to Java of a set of functions that were
 * specially written for the implementation of TestU01, a software package
 * for testing uniform random number generators @cite iLEC01t&thinsp;.
 *
 * Strictly speaking, applying several tests simultaneously makes the
 * @f$p@f$-values "invalid" in the sense that the probability of having *at
 * least one* @f$p@f$-value less than 0.01, say, is larger than 0.01. One
 * must therefore be careful with the interpretation of these @f$p@f$-values
 * (one could use, e.g., the Bonferroni inequality @cite sLAW00a&thinsp;).
 * Applying simultaneous tests is convenient in some situations, such as in
 * screening experiments for detecting statistical deficiencies in random
 * number generators. In that context, rejection of the null hypothesis
 * typically occurs with extremely small @f$p@f$-values (e.g., less than
 * @f$10^{-15}@f$), and the interpretation is quite obvious in this case.
 *
 * The class also provides tools to plot an empirical or theoretical
 * distribution function, by creating a data file that contains a graphic
 * plot in a format compatible with the software specified by the environment
 * variable  #graphSoft. NOTE: see also the more recent package
 * @ref umontreal.ssj.charts.
 *
 * Note: This class uses the Colt library.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class GofFormat {
   private GofFormat() {}

   /**
    * @name Plotting distribution functions
    * @{
    */

   /**
    * Data file format used for plotting functions with Gnuplot.
    */
   public static final int GNUPLOT = 0;

   /**
    * Data file format used for creating graphics with Mathematica.
    */
   public static final int MATHEMATICA = 1;

   /**
    * Environment variable that selects the type of software to be used
    * for plotting the graphs of functions. The data files produced by
    * #graphFunc and  #graphDistUnif will be in a format suitable for this
    * selected software. The default value is `GNUPLOT`. To display a
    * graphic in file `f` using `gnuplot`, for example, one can use the
    * command "<tt>plot f with steps, x with lines</tt>" in `gnuplot`.
    *  `graphSoft` can take the values  #GNUPLOT or  #MATHEMATICA.
    */
   public static int graphSoft = GNUPLOT;

   private static String formatMath2 (double x, double y)    {
      // Writes the pair (x, y) in file f, in a format understood
      // by Mathematica
      StringBuffer sb = new StringBuffer();
      String S;

      sb.append ("   { ");
      if ((x != 0.0) && (x < 0.1 || x > 1.0)) {
         S = PrintfFormat.E (16, 7, x);
         int exppos = S.indexOf ('E');
         if (exppos != -1)
            S = S.substring (0, exppos) + "*10^(" +
                             S.substring (exppos+1) + ")";
      }
      else
         S = PrintfFormat.g (16, 8, x);

      sb.append (S + ",     ");

      if (y != 0.0 && (y < 0.1 || y > 1.0)) {
         S = PrintfFormat.E (16, 7, y);
         int exppos = S.indexOf ('E');
         if (exppos != -1)
            S = S.substring (0, exppos) + "*10^(" +
                             S.substring (exppos+1) + ")";
      }
      else
        S = PrintfFormat.g (16, 8, y);

      sb.append (S + " }");
      return sb.toString();
   }


   private static String graphFunc (ContinuousDistribution dist, double a,
                                    double b, int m, int mono, String desc) {
// Renommer drawCDF en fixant mono = 1 et éliminant mono.
      int i;
      double yprec, y, x, h;
      StringBuffer sb = new StringBuffer();
      String openComment = "";
      String closeComment = "";
      String openGraph = "";
      String closeGraph = "";
      if (mono != 1 && mono != -1)
         throw new IllegalArgumentException ("mono must be 1 or -1");
      switch (graphSoft) {
      case GNUPLOT:
        openComment = "# ";
        closeComment = "";
        openGraph = "";
        closeGraph = PrintfFormat.NEWLINE;
        break;
      case MATHEMATICA:
        openComment = "(* ";
        closeComment = " *)";
        openGraph = "points = { " + PrintfFormat.NEWLINE;
        closeGraph = "}" + PrintfFormat.NEWLINE;
        break;
      }

      sb.append (openComment + "----------------------------------" +
                   closeComment  + PrintfFormat.NEWLINE);
      sb.append (openComment + PrintfFormat.s (-70, desc)
                 + closeComment  + PrintfFormat.NEWLINE +
                   PrintfFormat.NEWLINE);

      sb.append (openGraph);
      h = (b - a) / m;
      if (mono == 1)
         yprec = -Double.MAX_VALUE;
      else if (mono == -1)
         yprec = Double.MAX_VALUE;
      else
         yprec = 0.0;

      for (i = 0; i <= m; i++) {
         x = a + i*h;
         y = mono == 1 ? dist.cdf (x) : dist.barF (x);
         switch (graphSoft) {
         case MATHEMATICA:
            sb.append (formatMath2 (x, y));
            if (i < m)
               sb.append (',');
            break;
         default: // Default and GNUPLOT
            sb.append (PrintfFormat.g (20, 14, x) +  "      " +
                       PrintfFormat.g (20, 14, y));
         }

         switch (mono) {
         case 1:
            if (y < yprec)
               sb.append ("    " + openComment +
                    "  DECREASING" + closeComment);
            break;
         case -1:
            if (y > yprec)
               sb.append ("    " + openComment +
                    "  INCREASING" + closeComment);
            break;
         default:
            break;
         }
         sb.append (PrintfFormat.NEWLINE);
         yprec = y;
      }
      sb.append (closeGraph);
      return sb.toString();
   }

   /**
    * Formats data to plot the graph of the distribution function @f$F@f$
    * over the interval @f$[a,b]@f$, and returns the result as a  String.
    * The method `dist.cdf(x)` returns the value of @f$F@f$ at @f$x@f$.
    * The  String `desc` gives a short caption for the graphic plot. The
    * method computes the @f$m+1@f$ points @f$(x_i,  F (x_i))@f$, where
    * @f$x_i = a + i (b-a)/m@f$ for @f$i=0,1,…,m@f$, and formats these
    * points into a `String` in a format suitable for the software
    * specified by  #graphSoft. NOTE: see also the more recent class
    * @ref umontreal.ssj.charts.ContinuousDistChart.
    *  @param dist         continuous distribution function to plot
    *  @param a            lower bound of the interval to plot
    *  @param b            upper bound of the interval to plot
    *  @param m            number of points in the plot minus one
    *  @param desc         short caption describing the plot
    *  @return a string representation of the plot data
    */
   public static String drawCdf (ContinuousDistribution dist, double a,
                                 double b, int m, String desc) {
      return graphFunc (dist, a, b, m, 1, desc);
   }

   /**
    * Formats data to plot the graph of the density @f$f(x)@f$ over the
    * interval @f$[a,b]@f$, and returns the result as a  String. The
    * method `dist.density(x)` returns the value of @f$f(x)@f$ at @f$x@f$.
    * The  String `desc` gives a short caption for the graphic plot. The
    * method computes the @f$m+1@f$ points @f$(x_i,  f(x_i))@f$, where
    * @f$x_i = a + i (b-a)/m@f$ for @f$i=0,1,…,m@f$, and formats these
    * points into a `String` in a format suitable for the software
    * specified by  #graphSoft. NOTE: see also the more recent class
    * @ref umontreal.ssj.charts.ContinuousDistChart.
    *  @param dist         continuous density function to plot
    *  @param a            lower bound of the interval to plot
    *  @param b            upper bound of the interval to plot
    *  @param m            number of points in the plot minus one
    *  @param desc         short caption describing the plot
    *  @return a string representation of the plot data
    */
   public static String drawDensity (ContinuousDistribution dist, double a,
                                     double b, int m, String desc) {
      int i;
      double y, x, h;
      StringBuffer sb = new StringBuffer();
      String openComment = "";
      String closeComment = "";
      String openGraph = "";
      String closeGraph = "";

      switch (graphSoft) {
      case GNUPLOT:
        openComment = "# ";
        closeComment = "";
        openGraph = "";
        closeGraph = PrintfFormat.NEWLINE;
        break;
      case MATHEMATICA:
        openComment = "(* ";
        closeComment = " *)";
        openGraph = "points = { " + PrintfFormat.NEWLINE;
        closeGraph = "}" + PrintfFormat.NEWLINE;
        break;
      }

      sb.append (openComment + "----------------------------------" +
                   closeComment  + PrintfFormat.NEWLINE);
      sb.append (openComment + PrintfFormat.s (-70, desc)
                    + closeComment  + PrintfFormat.NEWLINE +
                      PrintfFormat.NEWLINE);

      sb.append (openGraph);
      h = (b - a) / m;

      for (i = 0; i <= m; i++) {
         x = a + i*h;
         y = dist.density (x);

         switch (graphSoft) {
         case MATHEMATICA:
            sb.append (formatMath2 (x, y));
            if (i < m)
               sb.append (',');
            break;
         default: // Default and GNUPLOT
            sb.append (PrintfFormat.g (16, 8, x) +  "      " +
                       PrintfFormat.g (16, 8, y));
         }
         sb.append (PrintfFormat.NEWLINE);
      }
      sb.append (closeGraph);
      return sb.toString();
   }

   /**
    * Formats data to plot the empirical distribution of
    * @f$U_{(1)},…,U_{(N)}@f$, which are assumed to be in `data[0...N-1]`,
    * and to compare it with the uniform distribution. The @f$U_{(i)}@f$
    * must be sorted. The two endpoints @f$(0, 0)@f$ and @f$(1, 1)@f$ are
    * always included in the plot. The string `desc` gives a short caption
    * for the graphic plot. The data is printed in a format suitable for
    * the software specified by  #graphSoft. NOTE: see also the more
    * recent class  @ref umontreal.ssj.charts.EmpiricalChart.
    *  @param data         array of observations to plot
    *  @param desc         short caption describing the plot
    *  @return a string representation of the plot data
    */
   public static String graphDistUnif (DoubleArrayList data, String desc) {
      double[] u = data.elements();
      int n = data.size();
      int i;
      double unSurN = 1.0/n;
      StringBuffer sb = new StringBuffer();

      switch (graphSoft) {
      case GNUPLOT:
         sb.append ("#----------------------------------" +
                     PrintfFormat.NEWLINE);
         sb.append ("# " + PrintfFormat.s (-70, desc) +
                     PrintfFormat.NEWLINE + PrintfFormat.NEWLINE);
         sb.append (PrintfFormat.g (16, 8, 0.0) + "  " +
                    PrintfFormat.g (16, 8, 0.0) + PrintfFormat.NEWLINE);
         for (i = 0; i < n; i++)
            sb.append (PrintfFormat.g (16, 8, u[i]) + "  " +
                       PrintfFormat.g (16, 8, (i + 1)*unSurN) +
                       PrintfFormat.NEWLINE);

         sb.append (PrintfFormat.g (16, 8, 1.0) + "  " +
                    PrintfFormat.g (16, 8, 1.0) + PrintfFormat.NEWLINE +
                    PrintfFormat.NEWLINE);
         break;
      case MATHEMATICA:
         sb.append ("(*----------------------------------*)" +
                     PrintfFormat.NEWLINE);
         sb.append ("(* " + PrintfFormat.s (-70, desc)  +
                     PrintfFormat.NEWLINE + " *)" +
                     PrintfFormat.NEWLINE + PrintfFormat.NEWLINE +
                     "points = { " + PrintfFormat.NEWLINE);

         sb.append (formatMath2 (0.0, 0.0) + "," + PrintfFormat.NEWLINE);
         for (i = 0; i < n; i++)
            sb.append (formatMath2 (u[i], (i + 1)*unSurN) + "," +
                       PrintfFormat.NEWLINE);
         sb.append (formatMath2 (1.0, 1.0) + PrintfFormat.NEWLINE);
         break;
      default:
         throw new IllegalArgumentException ("graphSoft unknown");
      }
      return sb.toString();
   }

   /**
    * @}
    */

   /**
    * @name Computing and printing \(p\)-values for EDF test statistics
    * @{
    */

   /**
    * Environment variable used in  #formatp0 to determine which
    * @f$p@f$-values are too close to 0 or 1 to be printed explicitly. If
    * `EPSILONP` @f$= \epsilon@f$, then any @f$p@f$-value less than
    * @f$\epsilon@f$ or larger than @f$1-\epsilon@f$ is *not* written
    * explicitly; the program simply writes "<tt>eps</tt>" or
    * "<tt>1-eps</tt>". The default value is @f$10^{-15}@f$.
    */
   public static double EPSILONP = 1.0E-15;

   /**
    * Environment variable used in  #formatp1 to determine which
    * @f$p@f$-values should be marked as suspect when printing test
    * results. If `SUSPECTP` @f$= \alpha@f$, then any @f$p@f$-value less
    * than @f$\alpha@f$ or larger than @f$1-\alpha@f$ is considered
    * suspect and is "singled out" by `formatp1`. The default value is
    * 0.01.
    */
   public static double SUSPECTP = 0.01;

   /**
    * Returns the @f$p@f$-value @f$p@f$ of a test, in the format
    * "@f$1-p@f$" if @f$p@f$ is close to 1, and @f$p@f$ otherwise. Uses
    * the environment variable  #EPSILONP and replaces @f$p@f$ by
    * @f$\epsilon@f$ when it is too small.
    *  @param p            the @f$p@f$-value to be formated
    *  @return the string representation of the @f$p@f$-value
    */
   public static String formatp0 (double p) {
      // Formats the p-value of a test, without a descriptor
      if ((p >= 0.01) && (p <= 0.99))
         return PrintfFormat.format (8, 2, 1, p);
      else if (p < EPSILONP)
         return "   eps  ";
      else if (p < 0.01)
         return PrintfFormat.format (8, 2, 2, p);
      else if (p >= 1.0 - EPSILONP)
         return " 1 - eps ";
      else
         return " 1 - " + PrintfFormat.g (8, 2, 1.0 - p);
   }

   /**
    * Returns the string "<tt>p-value of test : </tt>", then calls
    * #formatp0 to print @f$p@f$, and adds the marker "<tt>****</tt>" if
    * @f$p@f$ is considered suspect (uses the environment variable
    * `SUSPECTP` for this).
    *  @param p            the @f$p@f$-value to be formated
    *  @return the string representation of the p-value of test
    */
   public static String formatp1 (double p) {
      // Prints the p-value of a test, with a descriptor.
      StringBuffer sb = new StringBuffer();
      sb.append ("p-value of test                       :" + formatp0 (p));
      if (p < SUSPECTP || p > 1.0 - SUSPECTP)
         sb.append ("    *****");

      sb.append (PrintfFormat.NEWLINE + PrintfFormat.NEWLINE);
      return sb.toString();
   }

   /**
    * Returns `x` on a single line, then go to the next line and calls
    * #formatp1.
    *  @param x            value of the statistic for which the p-value is
    *                      formated
    *  @param p            the @f$p@f$-value to be formated
    *  @return the string representation of the p-value of test
    */
   public static String formatp2 (double x, double p) {
      // Prints the statistic x and its p-value p.
      return PrintfFormat.format (8, 2, 1, x) + PrintfFormat.NEWLINE +
             formatp1 (p);
   }

   /**
    * Formats the test statistic `x` for a test named `testName` with
    * @f$p@f$-value `p`. The first line of the returned string contains
    * the name of the test and the statistic whereas the second line
    * contains its p-value. The formated values of `x` and `p` are
    * aligned.
    *  @param testName     name of the test that was performed
    *  @param x            value of the test statistic
    *  @param p            @f$p@f$-value of the test
    *  @return the string representation of the test result
    */
   public static String formatp3 (String testName, double x, double p) {
      final String SLT = "p-value of test";
      int l = Math.max (SLT.length(), testName.length());
      PrintfFormat pf = new PrintfFormat();
      pf.append (-l, testName).append (" : ").append (8, 2, 1, x).append
               (PrintfFormat.NEWLINE);
      pf.append (-l, SLT).append (" : ").append (formatp0 (p));
      if (p < SUSPECTP || p > 1.0 - SUSPECTP)
         pf.append ("    *****");
      pf.append (PrintfFormat.NEWLINE + PrintfFormat.NEWLINE);
      return pf.toString();
   }

   /**
    * Computes the @f$p@f$-value of the chi-square statistic `chi2` for a
    * test with `k` intervals. Uses @f$d@f$ decimal digits of precision in
    * the calculations. The result of the test is returned as a string.
    * The @f$p@f$-value is computed using  GofStat.pDisc.
    *  @param k            number of subintervals for the chi-square test
    *  @param chi2         chi-square statistic
    *  @return the string representation of the test result and
    * @f$p@f$-value
    */
   public static String formatChi2 (int k, int d, double chi2) {
      StringBuffer sb = new StringBuffer();
      sb.append ("Chi2 statistic                        : " +
                  PrintfFormat.format (8, 2, 1, chi2));
      sb.append (PrintfFormat.NEWLINE +
                 "p-value                               : " +
                 formatp0 (GofStat.pDisc
                          (ChiSquareDist.cdf (k - 1, d, chi2),
                           ChiSquareDist.barF (k - 1, d, chi2))));
      sb.append (PrintfFormat.NEWLINE + PrintfFormat.NEWLINE);
      return sb.toString();
   }

   /**
    * Computes the @f$p@f$-values of the three Kolmogorov-Smirnov
    * statistics @f$D_N^+@f$, @f$D_N^-@f$, and @f$D_N@f$, whose values are
    * in `dp, dm, d`, respectively, assuming a sample of size `n`. Then
    * formats these statistics and their @f$p@f$-values using  #formatp2
    * for each one.
    *  @param n            sample size
    *  @param dp           value of the @f$D_N^+@f$ statistic
    *  @param dm           value of the @f$D_N^-@f$ statistic
    *  @param d            value of the @f$D_N@f$ statistic
    *  @return the string representation of the Kolmogorov-Smirnov
    * statistics and their p-values
    */
   public static String formatKS (int n, double dp,
                                  double dm, double d) {
      // Prints the results of a Kolmogorov-Smirnov test
      return "Kolmogorov-Smirnov+ statistic = D+    :" +
             formatp2 (dp, KolmogorovSmirnovPlusDist.barF (n, dp)) +
             "Kolmogorov-Smirnov- statistic = D-    :" +
             formatp2 (dm, KolmogorovSmirnovPlusDist.barF (n, dm)) +
             "Kolmogorov-Smirnov statistic = D      :" +
             formatp2 (d, KolmogorovSmirnovDistQuick.barF (n, d)) +
                       PrintfFormat.NEWLINE + PrintfFormat.NEWLINE;
   }

   /**
    * Computes the KS test statistics to compare the empirical
    * distribution of the observations in `data` with the theoretical
    * distribution `dist` and formats the results. See also method
    * {@link umontreal.ssj.gof.GofStat.kolmogorovSmirnov()
    * kolmogorovSmirnov(double[],ContinuousDistribution,double[],double[])}.
    *  @param data         array of observations to be tested
    *  @param dist         assumed distribution of the observations
    *  @return the string representation of the Kolmogorov-Smirnov
    * statistics and their p-values
    */
   public static String formatKS (DoubleArrayList data,
                                  ContinuousDistribution dist) {

      double[] v = data.elements();
      int n = data.size();

      DoubleArrayList dataUnif = GofStat.unifTransform (data, dist);
      dataUnif.quickSortFromTo (0, dataUnif.size() - 1);
      double[] ret = GofStat.kolmogorovSmirnov (dataUnif);
      return formatKS (n, ret[0], ret[1], ret[2]);
   }

   /**
    * Similar to  #formatKS(int,double,double,double), but for the KS
    * statistic @f$D_N^+(a)@f$ defined in (
    * {@link REF_gof_FDist_eq_KSPlusJumpOne KSPlusJumpOne}
    * ). Writes a header, computes the @f$p@f$-value and calls  #formatp2.
    *  @param n            sample size
    *  @param a            size of the jump
    *  @param dp           value of @f$D_N^+(a)@f$
    *  @return the string representation of the Kolmogorov-Smirnov
    * statistic and its p-value
    */
   public static String formatKSJumpOne (int n, double a, double dp) {
      double d = 1.0 - FDist.kolmogorovSmirnovPlusJumpOne (n, a, dp);

      return PrintfFormat.NEWLINE +
             "Kolmogorov-Smirnov+ statistic = D+    : " +
             PrintfFormat.g (8, 2, dp) + PrintfFormat.NEWLINE +
             formatp1 (d) + PrintfFormat.NEWLINE;
   }

   /**
    * Similar to  #formatKS(DoubleArrayList,ContinuousDistribution), but
    * for @f$D_N^+(a)@f$ defined in (
    * {@link REF_gof_FDist_eq_KSPlusJumpOne KSPlusJumpOne}
    * ).
    *  @param data         array of observations to be tested
    *  @param dist         assumed distribution of the data
    *  @param a            size of the jump
    *  @return string representation of the Kolmogorov-Smirnov statistic
    * and its p-value
    */
   public static String formatKSJumpOne (DoubleArrayList data,
                                         ContinuousDistribution dist,
                                         double a) {

      double[] v = data.elements();
      int n = data.size();
      DoubleArrayList dataUnif = GofStat.unifTransform (data, dist);
      dataUnif.quickSortFromTo (0, dataUnif.size() - 1);
      double[] ret =  GofStat.kolmogorovSmirnovJumpOne (dataUnif, a);
      return formatKSJumpOne (n, a, ret[0]);
   }

   /**
    * @}
    */

   /**
    * @name Applying several tests at once and printing results
    *
    * Higher-level tools for applying several EDF goodness-of-fit tests
    * simultaneously are offered here. The environment variable
    * `activeTests` specifies which tests in this list are to be performed
    * when asking for several simultaneous tests via the functions
    * `activeTests`, `formatActiveTests`, etc.
    *
    * @{
    */

   /**
    * Kolmogorov-Smirnov+ test
    */
   public  static final int KSP = 0;

   /**
    * Kolmogorov-Smirnov@f$-@f$ test
    */
   public static final int KSM = 1;

   /**
    * Kolmogorov-Smirnov test
    */
   public static final int KS = 2;

   /**
    * Anderson-Darling test
    */
   public static final int AD = 3;

   /**
    * Cramér-von Mises test
    */
   public static final int CM = 4;

   /**
    * Watson G test
    */
   public static final int WG = 5;

   /**
    * Watson U test
    */
   public static final int WU = 6;

   /**
    * Mean
    */
   public static final int MEAN = 7;

   /**
    * Correlation
    */
   public static final int COR = 8;

   /**
    * Total number of test types
    */
   public static final int NTESTTYPES = 9;

   /**
    * Name of each `testType` test. Could be used for printing the test
    * results, for example.
    */
   public static final String[] TESTNAMES = {
    "KolmogorovSmirnovPlus", "KolmogorovSmirnovMinus",
    "KolmogorovSmirnov", "Anderson-Darling",
    "CramerVon-Mises", "Watson G", "Watson U",
    "Mean", "Correlation"
   };

   /**
    * The set of EDF tests that are to be performed when calling the
    * methods  #activeTests,  #formatActiveTests, etc. By default, this
    * set contains `KSP`, `KSM`, and `AD`. Note: `MEAN` and `COR` are
    * *always excluded* from this set of active tests.
    *  The valid indices for this array are  #KSP,  #KSM,  #KS,  #AD,
    * #CM,  #WG,  #WU,  #MEAN, and  #COR.
    */
   public static boolean[] activeTests = null;
   private static void initActiveTests() {
      activeTests = new boolean[NTESTTYPES];
      for (int i = 0; i < activeTests.length; i++)
        activeTests[i] = false;
      activeTests[KSP] = activeTests[KSM] = true;
      activeTests[AD] = activeTests[MEAN] = activeTests[COR] = true;
   }
   static {
      initActiveTests();
   }

   /**
    * Computes all EDF test statistics enumerated above (except
    * <tt>COR</tt>) to compare the empirical distribution of
    * @f$U_{(0)},…,U_{(N-1)}@f$ with the uniform distribution, assuming
    * that these sorted observations are in `sortedData`. If @f$N > 1@f$,
    * returns `sVal` with the values of the KS statistics @f$D_N^+@f$,
    * @f$D_N^-@f$ and @f$D_N@f$, of the Cramér-von Mises statistic
    * @f$W_N^2@f$, Watson’s @f$G_N@f$ and @f$U_N^2@f$, Anderson-Darling’s
    * @f$A_N^2@f$, and the average of the @f$U_i@f$’s, respectively. If
    * @f$N = 1@f$, only puts @f$1 - {}@f$<tt>sortedData.get (0)</tt> in
    * `sVal[KSP]`. Calling this method is more efficient than computing
    * these statistics separately by calling the corresponding methods in
    * @ref GofStat.
    *  @param sortedData   array of sorted observations
    *  @param sVal         array that will be filled with the results of
    *                      the tests
    */
   public static void tests (DoubleArrayList sortedData, double[] sVal) {
      double[] u = sortedData.elements();
      int n = sortedData.size();
      int i;
      double a2 = 0.0, w2, dm = 0.0, dp = 0.0, w;
      double u1, ui, d2, d1;
      double sumZ;
      double unSurN;

      if (n <= 0)
        throw new IllegalArgumentException ("n <= 0");
      if (sVal.length != NTESTTYPES)
        throw new IllegalArgumentException ("sVal must " +
                              "be of size NTESTTYPES.");

      // We assume that u is already sorted.
      if (n == 1) {
         sVal[KSP] = 1.0 - u[0];
         sVal[MEAN] = u[0];
         return;
      }
      unSurN = 1.0 / n;
      w2 = unSurN / 12.0;
      sumZ = 0.0;
      for (i = 0; i < n; i++) {
         // Statistics KS
         d1 = u[i] - i*unSurN;
         d2 = (i + 1)*unSurN - u[i];
         if (d1 > dm)
            dm = d1;
         if (d2 > dp)
            dp = d2;
         // Watson U and G
         sumZ += u[i];
         w = u[i] - (i + 0.5)*unSurN;
         w2 += w*w;
         // Anderson-Darling
         ui = u[i];
         u1 = 1.0 - ui;
         if (ui < GofStat.EPSILONAD)
            ui = GofStat.EPSILONAD;
         else if (u1 < GofStat.EPSILONAD)
            u1 = GofStat.EPSILONAD;
         a2 += (2*i + 1) * Math.log (ui) + (1 + 2*(n - i - 1))*Math.log (u1);
      }
      if (dm > dp)
         sVal[KS] = dm;
      else
         sVal[KS] = dp;
      sVal[KSM] = dm;
      sVal[KSP] = dp;
      sumZ = sumZ * unSurN - 0.5;
      sVal[CM] = w2;
      sVal[WG] = Math.sqrt ((double) n) * (dp + sumZ);
      sVal[WU] = w2 - sumZ * sumZ * n;
      sVal[AD] = -n - a2 * unSurN;
      sVal[MEAN] = sumZ + 0.5;  // Nouveau ...
   }

   /**
    * The observations @f$V@f$ are in `data`, not necessarily sorted, and
    * their empirical distribution is compared with the continuous
    * distribution `dist`.  If @f$N = 1@f$, only puts `data.get (0)` in
    * `sVal[MEAN]`, and @f$1 - {}@f$<tt>dist.cdf (data.get (0))</tt> in
    * `sVal[KSP]`.
    *  @param data         array of observations to test
    *  @param dist         assumed distribution of the observations
    *  @param sVal         array that will be filled with the results of
    *                      the tests
    */
   public static void tests (DoubleArrayList data,
                             ContinuousDistribution dist, double[] sVal) {

      double[] v = data.elements();
      int n = data.size();

      if (n <= 0)
        throw new IllegalArgumentException ("n <= 0");

      DoubleArrayList sortedData = GofStat.unifTransform (data, dist);
      sortedData.quickSortFromTo (0, sortedData.size()-1);
      tests (sortedData, sVal);
      if (n == 1)
         sVal[MEAN] = v[0];     // On veut v[0], pas u[0].
   }

   /**
    * Computes the EDF test statistics by calling
    * #tests(DoubleArrayList,double[]), then computes the @f$p@f$-values
    * of those that currently belong to `activeTests`, and return these
    * quantities in `sVal` and `pVal`, respectively. Assumes that
    * @f$U_{(0)},…,U_{(N-1)}@f$ are in `sortedData` and that we want to
    * compare their empirical distribution with the uniform distribution.
    * If @f$N = 1@f$, only puts @f$1 - {}@f$<tt>sortedData.get (0)</tt> in
    * `sVal[KSP], pVal[KSP]`, and `pVal[MEAN]`.
    *  @param sortedData   array of sorted observations
    *  @param sVal         array that will be filled with the results of
    *                      the tests
    *  @param pVal         array that will be filled with the
    *                      @f$p@f$-values
    */
   public static void activeTests (DoubleArrayList sortedData,
                                   double[] sVal, double[] pVal) {

      double[] u = sortedData.elements();
      int n = sortedData.size();

      if (n <= 0)
        throw new IllegalArgumentException ("n <= 0");

      if (sVal.length != NTESTTYPES || pVal.length != NTESTTYPES)
        throw new IllegalArgumentException ("sVal and pVal must " +
              "be of length NTESTTYPES.");

      if (n == 1) {
         sVal[KSP] = 1.0 - u[0];
         pVal[KSP] = 1.0 - u[0];
         pVal[MEAN] = pVal[KSP];
         return;
      }
      // We assume that u is already sorted.
      tests (sortedData, sVal);

      if (activeTests.length != NTESTTYPES) {
        initActiveTests();
        System.err.println ("activeTests was invalid, it was reinitialized.");
      }

      if (activeTests[KSP])
         pVal[KSP] = KolmogorovSmirnovPlusDist.barF (n, sVal[KSP]);

      if (activeTests[KSM])
         pVal[KSM] = KolmogorovSmirnovPlusDist.barF (n, sVal[KSM]);

      if (activeTests[KS])
         pVal[KS] = KolmogorovSmirnovDistQuick.barF (n, sVal[KS]);

      if (activeTests[AD])
         pVal[AD] = AndersonDarlingDistQuick.barF (n, sVal[AD]);

      if (activeTests[CM])
         pVal[CM] = CramerVonMisesDist.barF (n, sVal[CM]);

      if (activeTests[WG])
         pVal[WG] = WatsonGDist.barF (n, sVal[WG]);

      if (activeTests[WU])
         pVal[WU] = WatsonUDist.barF (n, sVal[WU]);
   }

   /**
    * The observations are in `data`, not necessarily sorted, and we want
    * to compare their empirical distribution with the distribution
    * `dist`. If @f$N = 1@f$, only puts `data.get(0)` in `sVal[MEAN]`, and
    * @f$1 - {}@f$<tt>dist.cdf (data.get (0))</tt> in `sVal[KSP],
    * pVal[KSP]`, and `pVal[MEAN]`.
    *  @param data         array of observations to test
    *  @param dist         assumed distribution of the observations
    *  @param sVal         array that will be filled with the results of
    *                      the tests
    *  @param pVal         array that will be filled with the
    *                      @f$p@f$-values
    */
   public static void activeTests (DoubleArrayList data,
                                   ContinuousDistribution dist,
                                   double[] sVal, double[] pVal) {
      double[] v = data.elements();
      int n = data.size();

      if (n <= 0)
        throw new IllegalArgumentException ("n <= 0");

      DoubleArrayList sortedData = GofStat.unifTransform (data, dist);
      sortedData.quickSortFromTo (0, sortedData.size() - 1);

      activeTests (sortedData, sVal, pVal);
      if (n == 1)
         sVal[MEAN] = v[0];
   }

   /**
    * Gets the @f$p@f$-values of the *active* EDF test statistics, which
    * are in `activeTests`. It is assumed that the values of these
    * statistics and their @f$p@f$-values are *already computed*, in
    * `sVal` and `pVal`, and that the sample size is `n`. These statistics
    * and @f$p@f$-values are formated using  #formatp2 for each one. If
    * `n=1`, prints only `pVal[KSP]` using  #formatp1.
    *  @param n            sample size
    *  @param sVal         array containing the results of the tests
    *  @param pVal         array containing the @f$p@f$-values
    *  @return the results formated as a string
    */
   public static String formatActiveTests (int n, double[] sVal,
                                           double[] pVal) {

      if (activeTests.length != NTESTTYPES) {
        initActiveTests();
        System.err.println ("activeTests was invalid, it was reinitialized.");
      }
      if (sVal.length != NTESTTYPES || pVal.length != NTESTTYPES)
        throw new IllegalArgumentException ("The length of " +
           "sVal and pVal must be NTESTTYPES.");
      if (n == 1)
         return formatp1 (pVal[KSP]);

      StringBuffer sb = new StringBuffer (PrintfFormat.NEWLINE);
      if (activeTests[KSP])
         sb.append ("Kolmogorov-Smirnov+ statistic = D+    :" +
           formatp2 (sVal[KSP], pVal[KSP]));
      if (activeTests[KSM])
         sb.append ("Kolmogorov-Smirnov- statistic = D-    :" +
           formatp2 (sVal[KSM], pVal[KSM]));
      if (activeTests[KS])
         sb.append ("Kolmogorov-Smirnov statistic  = D     :" +
           formatp2 (sVal[KS], pVal[KS]));
      if (activeTests[AD])
         sb.append ("Anderson-Darling statistic = A2       :" +
           formatp2 (sVal[AD], pVal[AD]));
      if (activeTests[CM])
         sb.append ("Cramer-von Mises statistic = W2       :" +
           formatp2 (sVal[CM], pVal[CM]));
      if (activeTests[WG])
         sb.append ("Watson statistic = G                  :" +
           formatp2 (sVal[WG], pVal[WG]));
      if (activeTests[WU])
         sb.append ("Watson statistic = U2                 :" +
           formatp2 (sVal[WU], pVal[WU]));
      sb.append (PrintfFormat.NEWLINE);
      return sb.toString();
   }

   /**
    * Repeats the following `k` times: Applies the
    * GofStat.iterateSpacings transformation to the
    * @f$U_{(0)},…,U_{(N-1)}@f$, assuming that these observations are in
    * `sortedData`, then computes the EDF test statistics and calls
    * #activeTests(DoubleArrayList,double[],double[]) after each
    * transformation. The function returns the *original* array
    * `sortedData` (the transformations are applied on a copy of
    * <tt>sortedData</tt>). If `printval = true`, stores all the values
    * into the returned  String after each iteration. If `graph = true`,
    * calls  #graphDistUnif after each iteration to print to stream `f`
    * the data for plotting the distribution function of the @f$U_i@f$.
    *  @param sortedData   array containing the sorted observations
    *  @param k            number of times the tests are applied
    *  @param printval     if `true`, stores all the values of the
    *                      observations at each iteration
    *  @param graph        if `true`, the distribution of the @f$U_i@f$
    *                      will be plotted after each iteration
    *  @param f            stream where the plots are written to
    *  @return a string representation of the test results
    */
   public static String iterSpacingsTests (DoubleArrayList sortedData, int k,
                                           boolean printval, boolean graph,
                                           PrintWriter f) {

      int n = sortedData.size();

      DoubleArrayList sortedDataCopy = (DoubleArrayList)sortedData.clone();
      DoubleArrayList diffArrayList = new DoubleArrayList(sortedData.size()+2);

      int j;
      int i;
      double[] sVal = new double[NTESTTYPES], pVal = new double[NTESTTYPES];

      StringBuffer sb = new StringBuffer (PrintfFormat.NEWLINE);

      for (j = 1; j <= k; j++) {
         sb.append ("-----------------------------------" +
                     PrintfFormat.NEWLINE +
                     "EDF Tests after \"iterateSpacings\", level : " +
                     PrintfFormat.d (2, j) + PrintfFormat.NEWLINE);

         GofStat.diff (sortedDataCopy, diffArrayList, 0, n - 1, 0.0, 1.0);
         GofStat.iterateSpacings (sortedDataCopy, diffArrayList);
         sortedDataCopy.quickSortFromTo (0, sortedDataCopy.size() - 1);
         activeTests (sortedDataCopy, sVal, pVal);

         sb.append (formatActiveTests (n, sVal, pVal));
         String desc = "Values of Uniforms after iterateSpacings, level " +
             PrintfFormat.d (2, j);
         if (printval) {
          sb.append (desc + PrintfFormat.NEWLINE +
                     "------------------------" + PrintfFormat.NEWLINE);
          sb.append (sortedDataCopy + PrintfFormat.NEWLINE);
         }
         if (graph && f != null)
          f.print (graphDistUnif (sortedDataCopy, desc));
         else if (graph && f == null)
          sb.append (graphDistUnif (sortedDataCopy, desc));
       }
       return sb.toString();
   }

   /**
    * Similar to  #iterSpacingsTests, but with the  GofStat.powerRatios
    * transformation.
    *  @param sortedData   array containing the sorted observations
    *  @param k            number of times the tests are applied
    *  @param printval     if `true`, stores all the values of the
    *                      observations at each iteration
    *  @param graph        if `true`, the distribution of the @f$U_i@f$
    *                      will be plotted after each iteration
    *  @param f            stream where the plots are written to
    *  @return a string representation of the test results
    */
   public static String iterPowRatioTests (DoubleArrayList sortedData, int k,
                                           boolean printval, boolean graph,
                                           PrintWriter f) {

      int n = sortedData.size();
      DoubleArrayList sortedDataCopy = (DoubleArrayList)sortedData.clone();

      int i;
      int j;
      double[] sVal = new double[NTESTTYPES], pVal = new double[NTESTTYPES];

      StringBuffer sb = new StringBuffer (PrintfFormat.NEWLINE);

      for (i = 1; i <= k; i++) {
         GofStat.powerRatios (sortedDataCopy);
         sb.append ("-----------------------------------" +
                    PrintfFormat.NEWLINE +
                    "EDF Tests after \"powerRatios\", level : " +
                    PrintfFormat.d (2, i) + PrintfFormat.NEWLINE);

         sortedDataCopy.quickSortFromTo (0, sortedDataCopy.size() - 1);

         activeTests (sortedDataCopy, sVal, pVal);
         sb.append (formatActiveTests (n, sVal, pVal));
         String desc = "Values of Uniforms after PowerRatios, level " +
             PrintfFormat.d (2, i);
         if (printval) {
           sb.append (desc + PrintfFormat.NEWLINE +
                      "--------------------------" + PrintfFormat.NEWLINE);
           sb.append (sortedDataCopy + PrintfFormat.NEWLINE);
         }
         if (graph && f != null)
            f.print (graphDistUnif (sortedDataCopy, desc));
         else if (graph && f == null)
            sb.append (graphDistUnif (sortedDataCopy, desc));
      }
      return sb.toString();
   }
}

/**
 * @}
 */