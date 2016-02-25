/*
 * Class:        YListChart
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
package umontreal.ssj.charts;

/**
 * This class extends the class  @ref umontreal.ssj.charts.XYLineChart. Each
 * @ref YListChart object is associated with a
 * @ref umontreal.ssj.charts.YListSeriesCollection data set. The data is
 * given as one or more lists of @f$y@f$-coordinates. The @f$x@f$-coordinates
 * are regularly-spaced multiples of the indices of the data points.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class YListChart extends XYLineChart {

   /**
    * Empty constructor.
    */
   public YListChart() {
      super();
      // dataset = new XYListSeriesCollection();
      // init (null, null, null);
   }

   /**
    * Initializes a new `YListChart` instance with set of points `data`.
    * `title` is a title, `XLabel` is a short description of the
    * @f$x@f$-axis, and `YLabel` a short description of the @f$y@f$-axis.
    * The input vectors represents a set of plotting data. More
    * specifically, each vector `data` represents a @f$y@f$-coordinates
    * set. Position in the vector will form the @f$x@f$-coordinates.
    * Indeed, the value <tt>data</tt>@f$[j]@f$ corresponds to the point
    * @f$(j+1, \mathtt{data}[j])@f$ (but rescaled) on the chart.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series of point sets.
    */
   public YListChart (String title, String XLabel, String YLabel,
                      double[]... data) {
      super();
      dataset = new YListSeriesCollection(data);
      init (title, XLabel, YLabel);
   }

   /**
    * Similar to the constructor
    * {@link #YListChart(String,String,String,double[]) YListChart(title,
    * XLabel, YLabel, data)} above. Except that if `flag` is `true`, the
    * points are @f$(j+1, \mbox{\texttt{data}}[j])@f$ for each series; but
    * if `flag` is `false`, the points are @f$((j+1)/n,
    * \mbox{\texttt{data}}[j])@f$, where @f$n@f$ is the number of points
    * of each series in `data`.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param flag         to choose the step between @f$x@f$-coordinates
    *  @param data         series of point sets.
    */
   public YListChart (String title, String XLabel, String YLabel,
                      boolean flag, double[]... data) {
      super();
      dataset = new YListSeriesCollection(flag, data);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `YListChart` instance with a set of points `data`.
    * `title` is a title, `XLabel` is a short description of the
    * @f$x@f$-axis, and `YLabel` a short description of the @f$y@f$-axis.
    * The input vector represents a set of plotting data. Position in the
    * vector gives the @f$x@f$-coordinates of the curve. The value
    * <tt>data</tt>@f$[j]@f$ corresponds to the point @f$(j+1@f$,
    * <tt>data</tt>@f$[j]@f$) (but rescaled on the chart) for the curve.
    * However, only *the first* `numPoints` of `data` will be considered
    * to plot the curve.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         point set.
    *  @param numPoints    number of points to plot
    */
   public YListChart (String title, String XLabel, String YLabel,
                      double[] data, int numPoints) {
      super();
      dataset = new YListSeriesCollection(data, numPoints);
      init (title, XLabel, YLabel);
   }

   /**
    * Similar to the constructor
    * {@link #YListChart(String,String,String,double[],int)
    * YListChart(title, XLabel, YLabel, data, numPoints)} above, but the
    * points are @f$(h(j+1), \mbox{\texttt{data}}[j])@f$.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param h            step between @f$x@f$-coordinates
    *  @param data         point set.
    *  @param numPoints    number of points to plot
    */
   public YListChart (String title, String XLabel, String YLabel,
                      double h, double[] data, int numPoints) {
      super();
      dataset = new YListSeriesCollection(h, data, numPoints);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `YListChart` instance with set of points `data`.
    * `title` is a title, `XLabel` is a short description of the
    * @f$x@f$-axis, and `YLabel` a short description of the @f$y@f$-axis.
    * The input vectors represents a set of plotting data. More
    * specifically, for a @f$n@f$-row matrix `data`, each row
    * <tt>data</tt>@f$[i], i=0,…, n-1@f$, represents a @f$y@f$-coordinate
    * set for a curve. Position in the vector gives the
    * @f$x@f$-coordinates of the curves. Indeed, the value
    * <tt>data</tt>@f$[i][j]@f$ corresponds to the point @f$(j+1@f$,
    * <tt>data</tt>@f$[i][j]@f$) (but rescaled on the chart) for curve
    * @f$i@f$. However, only *the first* `numPoints` of each
    * <tt>data</tt>@f$[i]@f$ will be considered to plot each curve.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series of point sets.
    *  @param numPoints    number of points to plot
    */
   public YListChart (String title, String XLabel, String YLabel,
                      double[][] data, int numPoints) {
      super();
      dataset = new YListSeriesCollection(data, numPoints);
      init (title, XLabel, YLabel);
   }
}