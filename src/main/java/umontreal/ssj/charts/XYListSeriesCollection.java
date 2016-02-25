/*
 * Class:        XYListSeriesCollection
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

import   umontreal.ssj.functions.MathFunction;
import   umontreal.ssj.functionfit.SmoothingCubicSpline;
import   umontreal.ssj.util.RootFinder;

import   org.jfree.data.xy.*;
import   org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;

import   cern.colt.list.DoubleArrayList;

import   java.util.Locale;
import   java.util.Formatter;
import   java.awt.Color;

/**
 * This class extends  @ref umontreal.ssj.charts.SSJXYSeriesCollection. It
 * stores data used in a `XYLineChart` or in other related charts.
 * `XYListSeriesCollection` provides complementary tools to draw simple
 * curves; for example, one may add or remove plots series and modify plot
 * style. This class is linked with the JFreeChart `XYSeriesCollection` class
 * to store data plots, and linked with the JFreeChart
 * `XYLineAndShapeRenderer` to render the plot. Each series must contain
 * enough points to plot a nice curve. It is recommended to use about 30
 * points. However, some rapidly varying functions may require many more
 * points. This class can be used to draw scatter plots.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class XYListSeriesCollection  extends SSJXYSeriesCollection {
   protected String[] marksType;   // marks on points (+, x, *...)
   protected String[] dashPattern; // line dashing (solid, dotted, densely dotted, loosely dotted,
                                 //               dashed, densely dashed, loosely dashed, only marks)
   protected String[] plotStyle;   // plot style (lines, curves...)
   private boolean autoCompletion = false;

   /**
    * Creates a new `XYListSeriesCollection` instance with an empty
    * dataset.
    */
   public XYListSeriesCollection() {
      renderer = new XYLineAndShapeRenderer(true, false);
     // ((XYLineAndShapeRenderer)renderer).setShapesVisible(false);
      seriesCollection = new XYSeriesCollection();
   }

   /**
    * Creates a new `XYListSeriesCollection` instance with default
    * parameters and given data series. The input parameter `data`
    * represents a set of plotting data.
    *
    * For example, if one @f$n@f$-row matrix `data1` is given as argument,
    * then the first row <tt>data1</tt>@f$[0]@f$ represents the
    * @f$x@f$-coordinate vector, and every other row <tt>data1</tt>@f$[i],
    * i=1,…, n-1@f$, represents a @f$y@f$-coordinate set for the points.
    * Therefore matrix <tt>data1</tt>@f$[i][j]@f$, @f$i=0,…, n-1@f$,
    * corresponds to @f$n-1@f$ curves, all with the same
    * @f$x@f$-coordinates.
    *
    * However, one may want to plot several curves with different
    * @f$x@f$-coordinates. In that case, one should give the curves as
    * matrices with two rows. For examples, if the argument `data` is made
    * of three 2-row matrices `data1`, `data2` and `data3`, then they
    * represents three different curves, <tt>data*</tt>@f$[0]@f$ being the
    * @f$x@f$-coordinates, and <tt>data*</tt>@f$[1]@f$ the
    * @f$y@f$-coordinates of the curves.
    *
    * However, we may also consider the sets of points above not as part
    * of curves, but rather as several list of points.
    *  @param data         series of point sets.
    */
   public XYListSeriesCollection (double[][]... data) {
      renderer = new XYLineAndShapeRenderer(true, false);
   //   ((XYLineAndShapeRenderer)renderer).setShapesVisible(false);
      seriesCollection = new XYSeriesCollection();

      XYSeriesCollection tempSeriesCollection = (XYSeriesCollection)seriesCollection;
      for (int i = 0; i < data.length; i ++) {

         if (data[i].length < 2)
            throw new IllegalArgumentException (
               "Unable to render the plot. data["+ i +"] contains less than two rows");

         for (int j = 0; j < data[i].length-1; j++)
            if (data[i][j].length != data[i][j+1].length)
               throw new IllegalArgumentException(
                  "data["+ i +"][" + j + "] and data["+ i +"]["+ (j+1) +"] must share the same length");

         for (int j = 1; j < data[i].length; j++) {
            XYSeries serie = new XYSeries(" ");
            for (int k = 0; k < data[i][0].length; k++)
               serie.add(data[i][0][k], data[i][j][k]);
            tempSeriesCollection.addSeries(serie);
         }
      }

      // set default colors
      for(int i = 0; i < tempSeriesCollection.getSeriesCount(); i++)
         renderer.setSeriesPaint(i, getDefaultColor(i));

      // set default plot style
      plotStyle = new String[tempSeriesCollection.getSeriesCount()];
      marksType = new String[tempSeriesCollection.getSeriesCount()];
      dashPattern = new String[tempSeriesCollection.getSeriesCount()];
      for(int i = 0; i < tempSeriesCollection.getSeriesCount(); i++) {
         marksType[i] = " ";
         plotStyle[i] = "smooth";
         dashPattern[i] = "solid";
      }
   }

   /**
    * Creates a new `XYListSeriesCollection` instance with default
    * parameters and given points `data`. If `data` is a @f$n@f$-row
    * matrix, then the first row <tt>data</tt>@f$[0]@f$ represents the
    * @f$x@f$-coordinate vector, and every other row <tt>data</tt>@f$[i],
    * i=1,…, n-1@f$, represents a @f$y@f$-coordinate set of points.
    * Therefore, if the points represents curves to be plotted,
    * <tt>data</tt>@f$[i][ ]@f$, @f$i=0,…, n-1@f$, corresponds to
    * @f$n-1@f$ curves, all with the same @f$x@f$-coordinates. Only the
    * first `numPoints` of `data` will be considered for each of the set
    * of points.
    *  @param data         series of point sets.
    *  @param numPoints    Number of points to plot
    */
   public XYListSeriesCollection (double[][] data, int numPoints) {
      renderer = new XYLineAndShapeRenderer(true, false);
     // ((XYLineAndShapeRenderer)renderer).setShapesVisible(false);
      seriesCollection = new XYSeriesCollection();

      XYSeriesCollection tempSeriesCollection = (XYSeriesCollection)seriesCollection;
      if (data.length < 2)
         throw new IllegalArgumentException (
            "Unable to render the plot. data contains less than two rows");

      // n-1 curves: data[0] is x; data[i] is y for each curve
      for (int j = 1; j < data.length; j++) {
         XYSeries serie = new XYSeries(" ");
         for (int k = 0; k < numPoints; k++)
            serie.add(data[0][k], data[j][k]);
         tempSeriesCollection.addSeries(serie);
      }

      // set default colors
      for (int i = 0; i < tempSeriesCollection.getSeriesCount(); i++)
         renderer.setSeriesPaint(i, getDefaultColor(i));

      // set default plot style
      plotStyle = new String[tempSeriesCollection.getSeriesCount()];
      marksType = new String[tempSeriesCollection.getSeriesCount()];
      dashPattern = new String[tempSeriesCollection.getSeriesCount()];
      for (int i = 0; i < tempSeriesCollection.getSeriesCount(); i++) {
         marksType[i] = " ";
         plotStyle[i] = "smooth";
         dashPattern[i] = "solid";
      }
   }

   /**
    * Creates a new `XYListSeriesCollection` instance with default
    * parameters and given data. The input parameter represents a set of
    * data plots, the constructor will count the occurrence number @f$Y@f$
    * of each value @f$X@f$ in the `DoubleArrayList`, and plot the point
    * @f$(X, Y)@f$. Each  DoubleArrayList variable corresponds to a curve
    * on the chart.
    *  @param data         series of point sets.
    */
   public XYListSeriesCollection (DoubleArrayList... data) {
      renderer = new XYLineAndShapeRenderer(true, false);
     // ((XYLineAndShapeRenderer)renderer).setShapesVisible(false);
      seriesCollection = new XYSeriesCollection ();
      XYSeriesCollection tempSeriesCollection = (XYSeriesCollection )seriesCollection;

      XYSeries serie;
      double[] elements;
      int count = 0;
      DoubleArrayList temp;
      for(int i = 0; i < data.length; i++) {
         serie = new XYSeries(" ");

         temp = data[i].copy();  // deep copy
         temp.trimToSize();      // set capacity to the current size
         temp.quickSortFromTo(0, temp.size()-1);   // sort list in increasing order, simplify the next processings
         elements = temp.elements();

         int j = 0;
         int l = 0;
         while(j < elements.length) {
            while(j < elements.length && elements[j] == elements[l]) {
               j++;
               count++;
            }
            serie.add(elements[l], count);
            count = 0;
            l = j;
         }
         tempSeriesCollection.addSeries(serie);
      }

      // set default colors
      for(int i = 0; i < tempSeriesCollection.getSeriesCount(); i++) {
         renderer.setSeriesPaint(i, getDefaultColor(i));
      }

      // set default plot style
      plotStyle = new String[tempSeriesCollection.getSeriesCount()];
      marksType = new String[tempSeriesCollection.getSeriesCount()];
      dashPattern = new String[tempSeriesCollection.getSeriesCount()];
      for(int i = 0; i < tempSeriesCollection.getSeriesCount(); i++) {
         marksType[i] = " ";
         plotStyle[i] = "smooth";
         dashPattern[i] = "solid";
      }
   }

   /**
    * Creates a new `XYListSeriesCollection` instance with default
    * parameters and given data series. The input parameter represents a
    * set of plotting data. Each series of the given collection
    * corresponds to a curve on the plot.
    *  @param data         series of point sets.
    */
   public XYListSeriesCollection (XYSeriesCollection data) {
      renderer = new XYLineAndShapeRenderer(true, false);
    //  ((XYLineAndShapeRenderer)renderer).setShapesVisible(false);
      seriesCollection = data;
      for(int i = 0; i < data.getSeriesCount(); i++) {
         XYSeries serie = data.getSeries(i);
      }

      // set default colors
      for(int i = 0; i < data.getSeriesCount(); i++) {
         renderer.setSeriesPaint(i, getDefaultColor(i));
      }

      // set default plot style
      plotStyle = new String[data.getSeriesCount()];
      marksType = new String[data.getSeriesCount()];
      dashPattern = new String[data.getSeriesCount()];
      for(int i = 0; i < data.getSeriesCount(); i++) {
         marksType[i] = " ";
         plotStyle[i] = "smooth";
         dashPattern[i] = "solid";
      }
   }

   /**
    * @name Data control methods
    * @{
    */

   /**
    * Adds a data series into the series collection. Vector `x` represents
    * the @f$x@f$-coordinates and vector `y` represents the
    * @f$y@f$-coordinates of the series.
    *  @param x            @f$x_i@f$ coordinates.
    *  @param y            @f$y_i@f$ coordinates.
    *  @return Integer that represent the new point set’s position in the
    * JFreeChart `XYSeriesCollection` object.
    */
   public int add (double[] x, double[] y) {
      if (x.length != y.length)
         throw new IllegalArgumentException("x and y must have the same length");
      return add (x, y, x.length);
   }

   /**
    * Adds a data series into the series collection. Vector `x` represents
    * the @f$x@f$-coordinates and vector `y` represents the
    * @f$y@f$-coordinates of the series. Only *the first* `numPoints` of
    * `x` and `y` will be added to the new series.
    *  @param x            @f$x_i@f$ coordinates.
    *  @param y            @f$y_i@f$ coordinates.
    *  @param numPoints    Number of points to add
    *  @return Integer that represent the new point set’s position in the
    * JFreeChart `XYSeriesCollection` object.
    */
   public int add (double[] x, double[] y, int numPoints) {
      XYSeries serie = new XYSeries(" ");
      XYSeriesCollection tempSeriesCollection = (XYSeriesCollection )seriesCollection;
      serie.setNotify(true);
      if ((x.length < numPoints) ||(y.length < numPoints))
         throw new IllegalArgumentException("numPoints > length of x or y");
      for (int i = 0; i < numPoints; i++)
         serie.add(x[i], y[i]);
      tempSeriesCollection.addSeries(serie);

      // color
      int j = tempSeriesCollection.getSeriesCount()-1;
      renderer.setSeriesPaint(j, getDefaultColor(j));

      int co = tempSeriesCollection.getSeriesCount();
      String[] newPlotStyle = new String[co];
      String[] newMarksType = new String[co];
      String[] newDashPattern = new String[co];
      for(j = 0; j < co - 1; j++) {
         newPlotStyle[j] = plotStyle[j];
         newMarksType[j] = marksType[j];
         newDashPattern[j] = dashPattern[j];
      }

      newPlotStyle[j] = "smooth";
      newMarksType[j] = " ";
      newDashPattern[j] = "solid";
      plotStyle = newPlotStyle;
      marksType = newMarksType;
      dashPattern = newDashPattern;

      return tempSeriesCollection.getSeriesCount()-1;
   }

   /**
    * Adds a data series into the series collection. The input format of
    * `data` is described in constructor
    * `XYListSeriesCollection(double[][] data)`.
    *  @param data         input data.
    *  @return Integer that represent the number of point sets added to
    * the current dataset.
    */
   public int add (double[][] data) {
      return add (data, data[0].length);
   }

   /**
    * Adds data series into the series collection. The input format of
    * `data` is described in constructor
    * `XYListSeriesCollection(double[][] data)`. Only *the first*
    * `numPoints` of `data` (the first `numPoints` columns of the matrix)
    * will be added to each new series.
    *  @param data         input data.
    *  @param numPoints    Number of points to add for each new series
    *  @return Integer that represent the number of point sets added to
    * the current dataset.
    */
   public int add (double[][] data, int numPoints) {
      XYSeriesCollection tempSeriesCollection =
          (XYSeriesCollection) seriesCollection;
      int n = tempSeriesCollection.getSeriesCount();

      if (data.length < 2)
         throw new IllegalArgumentException(
            "Unable to render the plot. data contains less than two rows");

      for (int j = 0; j < data.length; j++)
         if (data[j].length < numPoints)
            throw new IllegalArgumentException(
               "data[" + j + "] has not enough points");

      for (int j = 1; j < data.length; j++) {
         XYSeries serie = new XYSeries(" ");
         serie.setNotify(true);
         for (int k = 0; k < numPoints; k++)
            serie.add(data[0][k], data[j][k]);
         tempSeriesCollection.addSeries(serie);
      }

      // color
      for(int j = n; j < tempSeriesCollection.getSeriesCount(); j++)
         renderer.setSeriesPaint(j, getDefaultColor(j));

      String[] newPlotStyle = new String[tempSeriesCollection.getSeriesCount()];
      String[] newMarksType = new String[tempSeriesCollection.getSeriesCount()];
      String[] newDashPattern = new String[tempSeriesCollection.getSeriesCount()];
      for (int j = 0; j < n; j++) {
         newPlotStyle[j] = plotStyle[j];
         newMarksType[j] = marksType[j];
         newDashPattern[j] = dashPattern[j];
      }

      for(int j = n; j < tempSeriesCollection.getSeriesCount(); j++) {
         newPlotStyle[j] = "smooth";
         newMarksType[j] = " ";
         newDashPattern[j] = "solid";
      }
      plotStyle = newPlotStyle;
      marksType = newMarksType;
      dashPattern = newDashPattern;

      return (tempSeriesCollection.getSeriesCount()-n);
   }

   /**
    * Adds a data series into the series collection. The input format of
    * `data` is described in constructor `XYListSeriesCollection
    * (DoubleArrayList... data)`.
    *  @param data         data series.
    *  @return Integer that represent the new point set’s position in the
    * JFreeChart `XYSeriesCollection` object.
    */
   public int add (DoubleArrayList data) {
      XYSeries serie = new XYSeries(" ");
      DoubleArrayList temp = data.copy();  // deep copy
      XYSeriesCollection tempSeriesCollection = (XYSeriesCollection) seriesCollection;

      temp.trimToSize();      // set capacity to the current size
      temp.quickSortFromTo(0, temp.size()-1);   // sort list in increasing order, simplify the next processings
      double[] elements = temp.elements();

      int count = 0;
      int j = 0;
      int l = 0;
      while(j < elements.length) {
         while(j < elements.length && elements[j] == elements[l]) {
            j++;
            count++;
         }
         serie.add(elements[l], count);
         count = 0;
         l = j;
      }
      tempSeriesCollection.addSeries(serie);

      // color
      j = tempSeriesCollection.getSeriesCount()-1;
      renderer.setSeriesPaint(j, getDefaultColor(j));

      String[] newPlotStyle = new String[tempSeriesCollection.getSeriesCount()];
      String[] newMarksType = new String[tempSeriesCollection.getSeriesCount()];
      String[] newDashPattern = new String[tempSeriesCollection.getSeriesCount()];
      for(j = 0; j < tempSeriesCollection.getSeriesCount()-1; j++) {
         newPlotStyle[j] = plotStyle[j];
         newMarksType[j] = marksType[j];
         newDashPattern[j] = dashPattern[j];
      }

      newPlotStyle[j] = "smooth";
      newMarksType[j] = " ";
      newDashPattern[j] = "solid";
      plotStyle = newPlotStyle;
      marksType = newMarksType;
      dashPattern = newDashPattern;

      return tempSeriesCollection.getSeriesCount()-1;
   }

   /**
    * Gets the current name of the selected series.
    *  @param series       series index.
    *  @return current name of the series.
    */
   public String getName (int series) {
      return (String)((XYSeriesCollection )seriesCollection).getSeries(series).getKey();
   }

   /**
    * Sets the name of the selected series.
    *  @param series       series index.
    *  @param name         point set new name.
    */
   public void setName (int series, String name) {
      if(name == null)
         name = " ";
      ((XYSeriesCollection)seriesCollection).getSeries(series).setKey(name);
   }

   /**
    * @}
    */

   /**
    * @name Rendering methods
    * @{
    */

   /**
    * Enables the auto completion option. When this parameter is enabled,
    * straight lines are used to approximate points on the chart bounds if
    * the method isn’t able to display all points, because the user
    * defined bounds are smaller than the most significant data point
    * coordinate, for instance. It does not extrapolate the point sets,
    * but simply estimates point coordinates on the curve at bound
    * positions for a better visual rendering.
    */
   public void enableAutoCompletion() {
      this.autoCompletion = true;
   }

   /**
    * Disables auto completion option. Default status is `disabled`.
    */
   public void disableAutoCompletion() {
      this.autoCompletion = false;
   }

   /**
    * Returns the mark type associated with the <tt>series</tt>th data
    * series.
    *  @param series       series index.
    *  @return mark type.
    */
   public String getMarksType (int series) {
      return marksType[series];
   }

   /**
    * Adds marks on the points of a data series. It is possible to use any
    * of the marks provided by the TikZ package, some of which are
    * "<tt>*</tt>", "<tt>+</tt>" and "<tt>x</tt>". A blank character, used
    * by default, disables marks. The PGF/TikZ documentation provides more
    * information about placing marks on plots.
    *  @param series       series index.
    *  @param marksType    mark type.
    */
   public void setMarksType (int series, String marksType) {
      this.marksType[series] = marksType;
   }

   /**
    * Returns the dash pattern associated with the <tt>series</tt>th data
    * series.
    *  @param series       series index.
    *  @return mark type.
    */
   public String getDashPattern (int series) {
      return dashPattern[series];
   }

   /**
    * Selects dash pattern for a data series. It is possible to use all
    * the dash options provided by the TikZ package: "<tt>solid</tt>",
    * "<tt>dotted</tt>", "<tt>densely dotted</tt>", "<tt>loosely
    * dotted</tt>", "<tt>dashed</tt>", "<tt>densely dashed</tt>",
    * "<tt>loosely dashed</tt>" and "<tt>only marks</tt>". If "<tt>only
    * marks</tt>" is chosen, then method  #setMarksType must be called to
    * choose the marks (which are blank by default).
    *  @param series       series index.
    *  @param dashPattern  dash style.
    */
   public void setDashPattern (int series, String dashPattern) {
      this.dashPattern[series] = dashPattern;
      if (dashPattern.equals("only marks")) {
          ((XYLineAndShapeRenderer) renderer).setSeriesLinesVisible(series, false);
          ((XYLineAndShapeRenderer) renderer).setSeriesShapesVisible(series, true);
      } else {
          ((XYLineAndShapeRenderer) renderer).setSeriesLinesVisible(series, true);
          ((XYLineAndShapeRenderer) renderer).setSeriesShapesVisible(series, false);
      }
   }

   /**
    * Gets the current plot style for the selected series.
    *  @param series       series index.
    *  @return current plot style.
    */
   public String getPlotStyle (int series) {
      return plotStyle[series];
   }

   /**
    * Selects the plot style for a given series. It is possible to use all
    * the plot options provided by the TikZ package. Some of which are:
    * "<tt>sharp plot</tt>", which joins points with straight lines,
    * "<tt>smooth</tt>", which joins points with a smoothing curve,
    * "<tt>only marks</tt>", which does not join points, etc. The PGF/TikZ
    * documentation provides more information about smooth plots, sharp
    * plots and comb plots.
    *  @param series       series index.
    *  @param plotStyle    plot style.
    */
   public void setPlotStyle (int series, String plotStyle) {
      this.plotStyle[series] = plotStyle;
   }
   public String toLatex (double XScale, double YScale,
                          double XShift, double YShift,
                          double xmin, double xmax,
                          double ymin, double ymax) {

      // Calcule les bornes reelles du graphique, en prenant en compte la position des axes
      xmin = Math.min(XShift, xmin);
      xmax = Math.max(XShift, xmax);
      ymin = Math.min(YShift, ymin);
      ymax = Math.max(YShift, ymax);

      Formatter formatter = new Formatter(Locale.US);
      XYSeriesCollection tempSeriesCollection = (XYSeriesCollection) seriesCollection;
      double XEPSILON = (1.0E-4/XScale)+XShift;
      double YEPSILON = (1.0E-4/YScale)+YShift;
      boolean outOfBounds = false;
      MathFunction[] spline = null;
      double[] xBounds = getRangeBounds();
      double[] yBounds = getDomainBounds();
      double x, y;
// Smoothing splines, consulter  ref: QA278.2 G74, QA278.2 T35, QA278.2 E87

//       if(xBounds[0] < xmin || xBounds[1] > xmax || yBounds[0] < ymin || yBounds[1] > ymax) {
//          // on sait qu'il y a des points qui vont sortir du chart
//          // initialisation d'une spline pour chaque courbe
//          spline = new SmoothingCubicSpline[seriesCollection.getSeriesCount()];
//          for(int i = 0; i<seriesCollection.getSeriesCount(); i++)
//             spline[i] = new SmoothingCubicSpline(  (seriesCollection.getSeries(i).toArray())[0],
//                                                    (seriesCollection.getSeries(i).toArray())[1], 0.5);
//       }

      // on sait qu'il y a des points qui vont sortir du chart
      // initialisation d'une spline pour chaque courbe
      if (true) {
         spline = new SmoothingCubicSpline[tempSeriesCollection.getSeriesCount()];
         for (int i = 0; i < tempSeriesCollection.getSeriesCount(); i++)
            spline[i] = new SmoothingCubicSpline((tempSeriesCollection.getSeries(i).toArray())[0],
                                                 (tempSeriesCollection.getSeries(i).toArray())[1], 1);
      } else {
         spline = new AffineFit[tempSeriesCollection.getSeriesCount()];
         for (int i = 0; i < tempSeriesCollection.getSeriesCount(); i++)
            spline[i] = new AffineFit((tempSeriesCollection.getSeries(i).toArray())[0],
                                      (tempSeriesCollection.getSeries(i).toArray())[1]);
      }

      for(int i = tempSeriesCollection.getSeriesCount()-1; i >= 0; i--) {
         XYSeries temp = tempSeriesCollection.getSeries(i);

         if (temp.getItemCount() < 2)
            throw new IllegalArgumentException("Unable to plot series " + i +
                                               ": this series must have two points at least");

         Color color = (Color)renderer.getSeriesPaint(i);
         String colorString = detectXColorClassic(color);
         if( colorString == null) {
            colorString = "color"+i;
            formatter.format( "\\definecolor{%s}{rgb}{%.2f, %.2f, %.2f}%n",
                              colorString, color.getRed()/255.0, color.getGreen()/255.0, color.getBlue()/255.0);
         }

         // Cas particulier pour le premier point, on doit savoir si il est dans le chart ou pas
         if (  temp.getX(0).doubleValue() >= xmin && temp.getX(0).doubleValue() <= xmax &&
               temp.getY(0).doubleValue() >= ymin && temp.getY(0).doubleValue() <= ymax) {
            outOfBounds = false;
            formatter.format( "\\draw [%s, color=%s, mark=%s, style=%s] plot coordinates {%%%n",
                              plotStyle[i], colorString, marksType[i], dashPattern[i]);
         }
         else {
            outOfBounds = true;
            formatter.format("%% ");
         }
         formatter.format("(%.2f,%.4f)",   (temp.getX(0).doubleValue()-XShift)*XScale,
                                              (temp.getY(0).doubleValue()-YShift)*YScale);
         formatter.format(" %%   (%f,  %f)%n", temp.getX(0).doubleValue(), temp.getY(0).doubleValue());

         // Cas general
         for(int j = 1; j < temp.getItemCount(); j++) {
            double[] result;
            if (!outOfBounds) { //on est dans le chart
               result = evalLimitValues(xmin, xmax, ymin, ymax, XEPSILON, YEPSILON, spline[i], temp, j, false);
               // on regarde si on ne sort pas du chart, si c'est le cas on evalue le point en limite
               if (result != null) { // le point courant n'est pas dans le chart, on sort donc du chart
                  outOfBounds = true;
                  if (autoCompletion)
                     formatter.format("(%.2f,%.4f) %%%n", (result[0]-XShift)*XScale, (result[1]-YShift)*YScale);
                  formatter.format("}%%%n%% ");
               }
            }
            else { // le point precedent etait hors du chart
               if (  temp.getX(j).doubleValue() >= xmin && temp.getX(j).doubleValue() <= xmax &&
                     temp.getY(j).doubleValue() >= ymin && temp.getY(j).doubleValue() <= ymax) {
                     // on rentre dans le chart, il faut evaluer le point en limite
                  j = j-1;
                  result = evalLimitValues(xmin, xmax, ymin, ymax, XEPSILON, YEPSILON, spline[i], temp, j, true);
                  // ici result ne peut pas etre null
                  formatter.format( ";%%%n\\draw [%s, color=%s, mark=%s, style=%s] plot coordinates {%%%n",
                                    plotStyle[i], colorString, marksType[i], dashPattern[i]);
                  if (autoCompletion)
                     formatter.format("(%.2f,%.4f) %%%n ", (result[0]-XShift)*XScale, (result[1]-YShift)*YScale);
                  formatter.format("%% ");
                  outOfBounds = false;
               }
               else {
                  formatter.format("%% ");
                  // on les donnees sont toujours hors du chart
               }
            }
            /* on affiche les coordonnees du point quoiqu'il arrive,
            si celui ci est hors du chart alors la balise de commentaire a ete deja place */
            formatter.format("(%.2f,%.4f)",   (temp.getX(j).doubleValue()-XShift)*XScale,
                                              (temp.getY(j).doubleValue()-YShift)*YScale);
            if(j == temp.getItemCount()-1)
               formatter.format("}");
            formatter.format(" %%   (%f,  %f)%n", temp.getX(j).doubleValue(), temp.getY(j).doubleValue());
//            formatter.format(" %%%n");
         }
         formatter.format(" node[right] {%s};%n", (String)temp.getKey());
      }
      return formatter.toString();
   }


   /**
    * Compute x and y to chart limit bounds for extra-bounded points
    *
    * @param   xmin     lower bound for x coordinates
    * @param   xmax     upper bound for x coordinates
    * @param   ymin     lower bound for y coordinates
    * @param   ymax     upper bound for y coordinates
    * @param   XEPSILON increment step size for x coordinates
    * @param   YEPSILON increment step size for y coordinates
    * @param   spline   sline used to approximate points
    * @param   temp     current series
    * @param   numPoint point index in the current series
    * @param   sens     direction of the in chart last point
    *                    true  : point numPoint+1 is in the chart
    *                    false : point numPoint-1 is in the chart
    *
    * @return           x and y coordinates on the chart bounds.
    */
   private static double[] evalLimitValues(double xmin, double xmax, double ymin, double ymax, double XEPSILON, double YEPSILON, MathFunction spline, XYSeries temp, int numPoint, boolean sens) {
      int j = numPoint;
      int k = 0;
      double x, y;
      if(sens)
         k = j+1;
      else
         k = j-1;
      if(temp.getX(j).doubleValue() < xmin) {// Hors du chart mais on etait dans le chart au point precedent
         x = xmin;
         y = spline.evaluate(xmin); // spline puis evaluer en xmin
         while(y<ymin) {
            x += XEPSILON;
            y = spline.evaluate(x);
         }  // evaluer un x>xmin tantque y<ymin, y peut etre superieur a ymin car le point precedent est dans le chart
         while(y > ymax) {
            x += XEPSILON;
            y = spline.evaluate(x);
         }  // evaluer un x en ymax avec x > xmin
      }
      else if(temp.getX(j).doubleValue() > xmax) {
         x = xmax;
         y = spline.evaluate(xmax);
         while(y<ymin) {
            x -= XEPSILON;
            y = spline.evaluate(x);
         }  // evaluer un x<xmax tantque y<ymin
         while(y > ymax) {
            x -= XEPSILON;
            y = spline.evaluate(x);
         }  // evaluer un x<xmax tantque y>ymax
      }
      else if(temp.getY(j).doubleValue() < ymin) {
         y = ymin;
         x = evaluateX(spline, y, temp.getX(j).doubleValue(), temp.getX(k).doubleValue());// spline puis evaluer en ymin avec x ente xValue(ptCourant) et xValue(ptCourant-1)
         while(x < xmin) {
            y += YEPSILON;
            x = evaluateX(spline, y, x, temp.getX(k).doubleValue());
         }
         while(x > xmax) {
           y += YEPSILON;
           x = evaluateX(spline, y, x, temp.getX(k).doubleValue());
         }
      }
      else if(temp.getY(j).doubleValue() > ymax) {
         y = ymax;
         x = evaluateX(spline, y, temp.getX(j).doubleValue(), temp.getX(k).doubleValue());// spline puis evaluer en ymax avec x ente xValue(ptCourant) et xValue(ptCourant-1)
         while(x < xmin) {
            y -= YEPSILON;
            x = evaluateX(spline, y, x, temp.getX(k).doubleValue());
         }
         while(x > xmax) {
            y -= YEPSILON;
            x = evaluateX(spline, y, x, temp.getX(k).doubleValue());
         }
      }
      else
         return null;
      double[] retour = new double[2];
      retour[0] = x;
      retour[1] = y;
      return retour;
   }


   private static double evaluateX (final MathFunction spline, final double y, double xPrincipal, double xAnnexe) {
      final MathFunction xFunction = new MathFunction () {
         public double evaluate (double t) {
            return spline.evaluate(t) - y;
         }
      };
      return RootFinder.brentDekker (xPrincipal, xAnnexe-1.0E-6, xFunction, 1e-6);
   }


   private class AffineFit implements MathFunction{

      double[] x;
      double[] y;

      public AffineFit(double[] x, double[] y) {
         this.x =x;
         this.y =y;
      }

      public double evaluate(double t) {
         int i = 0;
         if (t <= x[0])
            return y[0];
         while (i < x.length && t > x[i])
            i++;
         i--;
         if (i == x.length)
            return x[x.length-1];

         return y[i] + ((t - x[i]) / (x[i+1] - x[i])) * (y[i+1] - y[i]);
      }
   }
}

/**
 * @}
 */