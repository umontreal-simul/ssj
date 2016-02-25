/*
 * Class:        YListSeriesCollection
 * Description:  Lists of y-coordinates of charts
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

import   org.jfree.data.xy.XYSeries;
import   org.jfree.data.xy.XYSeriesCollection;
import   org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;

/**
 * This class extends  @ref umontreal.ssj.charts.XYListSeriesCollection. The
 * data is given as lists of @f$y@f$-coordinates. The @f$x@f$-coordinates are
 * regularly spaced multiples of the indices of the data points.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class YListSeriesCollection extends XYListSeriesCollection {

   private void initYListSeries (double h, double[] data, int numPoints)
   {
      renderer = new XYLineAndShapeRenderer(true, false);
      seriesCollection = new XYSeriesCollection();

      XYSeriesCollection tempSeriesCollection =
         (XYSeriesCollection)seriesCollection;
      XYSeries serie = new XYSeries(" ");
      for (int j = 0; j < numPoints; j++)
         serie.add(h*(j+1), data[j]);
      tempSeriesCollection.addSeries(serie);

      // set default colors
      renderer.setSeriesPaint(0, getDefaultColor(0));

      // set default plot style
      plotStyle = new String[1];
      marksType = new String[1];
      dashPattern = new String[1];
      marksType[0] = " ";
      plotStyle[0] = "smooth";
      dashPattern[0] = "solid";
   }


   private void initYListSeries (boolean flag, double[]... data)
   {
      // if flag = true, h = 1; else h = 1/numPoints
      double h;
      renderer = new XYLineAndShapeRenderer(true, false);
      seriesCollection = new XYSeriesCollection();

      XYSeriesCollection tempSeriesCollection =
         (XYSeriesCollection)seriesCollection;
      for (int i = 0; i < data.length; i ++) {
         XYSeries serie = new XYSeries(" ");
         if (flag)
            h = 1;
         else
            h = 1.0 / data[i].length;
         for (int j = 0; j < data[i].length; j++)
            serie.add(h*(j+1), data[i][j]);
         tempSeriesCollection.addSeries(serie);
      }

      final int s = tempSeriesCollection.getSeriesCount();

      // set default colors
      for(int i = 0; i < s; i++) {
         renderer.setSeriesPaint(i, getDefaultColor(i));
      }

      // set default plot style
      plotStyle = new String[s];
      marksType = new String[s];
      dashPattern = new String[s];
      for (int i = 0; i < s; i++) {
         marksType[i] = " ";
         plotStyle[i] = "smooth";
         dashPattern[i] = "solid";
      }
   //   dashPattern[s-1] = "dashed";     // for the line y = x
   }

   /**
    * Creates a new `YListSeriesCollection` instance with default
    * parameters and given data series. The input vectors represent sets
    * of plotting data. More specifically, each vector `data` represents a
    * @f$y@f$-coordinates set. Position in the vector will form the
    * @f$x@f$-coordinates. Indeed the value <tt>data</tt>@f$[j]@f$
    * corresponds to the point @f$(j+1, \mbox{\texttt{data}}[j])@f$ on the
    * chart.
    *  @param data         series of point sets.
    */
   public YListSeriesCollection (double[]... data) {
      initYListSeries (true, data);
   }

   /**
    * Similar to the constructor  {@link #YListSeriesCollection(double[])
    * YListSeriesCollection(data)} above, except that if `flag` is `true`,
    * the points are @f$(j+1, \mbox{\texttt{data}}[j])@f$ for each series;
    * but if `flag` is `false`, the points are @f$((j+1)/n,
    * \mbox{\texttt{data}}[j])@f$, where @f$n@f$ is the number of points
    * of each series in `data`.
    *  @param flag         to choose the step between @f$x@f$-coordinates
    *  @param data         series of point sets.
    */
   public YListSeriesCollection (boolean flag, double[]... data) {
      initYListSeries (flag, data);
   }

   /**
    * Creates a new `YListSeriesCollection` instance with default
    * parameters and one data series. The vector `data` represents the
    * @f$y@f$-coordinate of the points, and position in the vector
    * represents the @f$x@f$-coordinate. However, only *the first*
    * `numPoints` of `data` will be considered in the series. Thus the
    * coordinates of the points are given by @f$(j,
    * \mbox{\texttt{data}}[j-1])@f$, for @f$j=1,2,…,
    * \mathtt{numPoints}@f$.
    *  @param data         point set.
    *  @param numPoints    Number of points to plot
    */
   public YListSeriesCollection (double[] data, int numPoints) {
      initYListSeries (1, data, numPoints);
   }

   /**
    * Similar to the constructor
    * {@link #YListSeriesCollection(double[],int)
    * YListSeriesCollection(data, numPoints)} above, but the points are
    * @f$(hj, \mbox{\texttt{data}}[j-1])@f$, for @f$j=1,2,…,
    * \mathtt{numPoints}@f$.
    *  @param h            step between @f$x@f$-coordinates
    *  @param data         point set.
    *  @param numPoints    Number of points to plot
    */
   public YListSeriesCollection (double h, double[] data, int numPoints) {
      initYListSeries (h, data, numPoints);
   }

   /**
    * Creates a new `YListSeriesCollection` instance with default
    * parameters and given data series. The matrix `data` represents a set
    * of plotting data. More specifically, each row of `data` represents a
    * @f$y@f$-coordinates set. Position in the vector will form the
    * @f$x@f$-coordinates. Indeed, for each serie @f$i@f$, the value
    * <tt>data</tt>@f$[i][j]@f$ corresponds to the point @f$(j+1,
    * \mbox{\texttt{data}}[j])@f$ on the chart. However, only *the first*
    * `numPoints` of `data` will be considered for each series of points.
    *  @param data         series of point sets.
    *  @param numPoints    Number of points to plot
    */
   public YListSeriesCollection (double[][] data, int numPoints) {
      renderer = new XYLineAndShapeRenderer(true, false);
      seriesCollection = new XYSeriesCollection();

      XYSeriesCollection tempSeriesCollection =
         (XYSeriesCollection)seriesCollection;
      for (int i = 0; i < data.length; i ++) {
         XYSeries serie = new XYSeries(" ");
         for (int j = 0; j < numPoints; j++)
            serie.add(j + 1, data[i][j]);
         tempSeriesCollection.addSeries(serie);
      }

      final int s = tempSeriesCollection.getSeriesCount();

      // set default colors
      for (int i = 0; i < s; i++) {
         renderer.setSeriesPaint(i, getDefaultColor(i));
      }

      // set default plot style
      plotStyle = new String[s];
      marksType = new String[s];
      dashPattern = new String[s];
      for (int i = 0; i < s; i++) {
         marksType[i] = " ";
         plotStyle[i] = "smooth";
         dashPattern[i] = "solid";
      }
   }
}