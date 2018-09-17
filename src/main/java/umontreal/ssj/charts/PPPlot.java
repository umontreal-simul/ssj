/*
 * Class:        PPPlot
 * Description:  pp-plot
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        May 2011
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
package umontreal.ssj.charts;
   import umontreal.ssj.probdist.ContinuousDistribution;
   import java.util.Arrays;

/**
 * This class implements *PP-plot* (or probability-probability plot) objects
 * that compare two probability distributions. The data is given as a list of
 * @f$x@f$-coordinates @f$(x_1, x_2, …, x_n)@f$, and one is given a reference
 * continuous probability distribution @f$F(x)@f$. One first sorts the
 * @f$x_i@f$ in ascending order, then noted @f$x_{(i)}@f$, and plots the
 * points @f$(i/n, F(x_{(i)}))@f$, @f$i= 1, 2, …, n@f$, to see if the data
 * @f$x_i@f$ comes from the reference distribution @f$F(x)@f$. The graph of
 * the straight line @f$y=x@f$ is also plotted for comparison.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class PPPlot extends XYLineChart {
   private double[][] U;        // u_i = cdf(x_i)
   private double[][] Lin;      // line y = x

   private void initLinear ()
   {
      // line y = x in [0, 1] by steps of h
      int m = 100;
      double h = 1.0 / m;
      Lin = new double[2][m+1];
      for (int i = 0; i <= m; i++)
         Lin[0][i] = Lin[1][i] = h * i;
   }


   private void initPoints (ContinuousDistribution dist, double[] data,
                            int numPoints)
   {
      int i;
      U = new double[2][numPoints];

      for (i = 0; i < numPoints; i++)
         U[1][i] = dist.cdf(data[i]);
      Arrays.sort(U[1]);
      for (i = 0; i < numPoints; i++)
         U[0][i] = (double) (i+1)/numPoints;
   }

   /**
    * Initializes a new `PPPlot` instance using the points `X`. `title` is
    * a title, `XLabel` is a short description of the @f$x@f$-axis, and
    * `YLabel` a short description of the @f$y@f$-axis. The plot is a
    * PP-plot of the points @f$(i/n, F(x_{(i)})@f$, @f$i= 1, 2, …, n@f$,
    * where @f$x_i = @f$<tt> X[</tt>@f$i@f$<tt>-1]</tt>, @f$x_{(i)}@f$ are
    * the sorted points, and @f$F(x) = @f$<tt>
    * dist.cdf(</tt>@f$x@f$<tt>)</tt>. The points `X` are not sorted.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param dist         Reference distribution
    *  @param X            points.
    */
   public PPPlot (String title, String XLabel, String YLabel,
                  ContinuousDistribution dist, double[] X) {
      this (title, XLabel, YLabel, dist, X, X.length);
   }

   /**
    * Similar to the constructor
    * {@link #PPPlot(String,String,String,ContinuousDistribution,double[])
    * PPPlot(title, XLabel, YLabel, dist, X)} above, except that only *the
    * first* `numPoints` of `X` are plotted.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param dist         Reference distribution
    *  @param X            point set.
    *  @param numPoints    number of points to plot
    */
   public PPPlot (String title, String XLabel, String YLabel,
                  ContinuousDistribution dist, double[] X, int numPoints) {
      super();
      initPoints (dist, X, numPoints);
      initLinear ();
      dataset = new XYListSeriesCollection(U, Lin);
      // --- dashed line for y = x
      ((XYListSeriesCollection)dataset).setDashPattern(1, "dashed");
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `PPPlot` instance. `title` is a title, `XLabel` is
    * a short description of the @f$x@f$-axis, and `YLabel` a short
    * description of the @f$y@f$-axis. The input vectors in `data`
    * represents several sets of @f$x@f$-points. @f$r@f$ determine the set
    * of points to be plotted in the PP-plot, that is, one will plot only
    * the points `data[r][i]`, for @f$i=0, 1, …, (n-1)@f$ and a given
    * @f$r@f$, where @f$n@f$ is the number of points in set @f$r@f$. The
    * points are assumed to follow the distribution `dist`.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param dist         Reference distribution
    *  @param data         series of point sets.
    *  @param r            set of points to plot
    */
   public PPPlot (String title, String XLabel, String YLabel,
                  ContinuousDistribution dist, double[][] data, int r) {
      this (title, XLabel, YLabel, dist, data[r], data[r].length);
   }
}