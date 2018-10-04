/*
 * Class:        XYLineChart
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

import   org.jfree.chart.*;
import   org.jfree.chart.axis.*;
import   org.jfree.chart.plot.*;
import   org.jfree.chart.renderer.xy.*;
import   org.jfree.data.xy.*;
// import   org.jfree.data.category.*;
// import   org.jfree.chart.renderer.category.*;
import   java.util.Locale;
import   java.util.Formatter;
import   java.lang.Math;
import   java.awt.*;
import   java.awt.geom.*;
import   cern.colt.list.DoubleArrayList;
import   javax.swing.JFrame;

/**
 * This class provides tools to create and manage curve plots. Using the
 * @ref XYLineChart class is the simplest way to produce curve plots only.
 * Each  @ref XYLineChart object is linked with a
 * @ref umontreal.ssj.charts.XYListSeriesCollection data set.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class XYLineChart extends XYChart {

   protected void init (String title, String XLabel, String YLabel) {
      // create the chart...
      chart = ChartFactory.createXYLineChart (
         title,                    // chart title
         XLabel,                   // x axis label
         YLabel,                   // y axis label
         dataset.getSeriesCollection(), // data
         PlotOrientation.VERTICAL,
         false,                    // include legend
         true,                     // tooltips
         false                     // urls
      );

      if (null != title) {
         if (title.startsWith("cdf") || title.startsWith("prob") || title.startsWith("density"))
            setprobFlag (true);
      }

      ((XYPlot)chart.getPlot()).setRenderer(dataset.getRenderer());
      // Initialize axis variables
      initAxis();
   }

   protected void initAxis(){
      XAxis = new Axis((NumberAxis)((XYPlot)chart.getPlot()).getDomainAxis(),
                        Axis.ORIENTATION_HORIZONTAL);
      YAxis = new Axis((NumberAxis)((XYPlot)chart.getPlot()).getRangeAxis(),
                        Axis.ORIENTATION_VERTICAL);
      setAutoRange(true, true);
   }

   /**
    * Initializes a new `XYLineChart` instance with an empty data set.
    */
   public XYLineChart() {
      super();
      dataset = new XYListSeriesCollection();
      init (null, null, null);
   }

   /**
    * Initializes a new `XYLineChart` instance with sets of points `data`.
    * `title` is a title, `XLabel` is a short description of the
    * @f$x@f$-axis, and `YLabel` a short description of the @f$y@f$-axis.
    * The input parameter `data` represents a set of plotting data.
    *
    * For example, if one @f$n@f$-row matrix `data1` is given as argument
    * `data`, then the first row <tt>data1</tt>@f$[0]@f$ represents the
    * @f$x@f$-coordinate vector, and every other row <tt>data1</tt>@f$[i],
    * i=1,…, n-1@f$, represents a @f$y@f$-coordinate set for a curve.
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
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series of point sets.
    */
   public XYLineChart (String title, String XLabel, String YLabel,
                       double[][]... data) {
      super();
      dataset = new XYListSeriesCollection(data);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `XYLineChart` instance with sets of points `data`.
    * `title` is a title, `XLabel` is a short description of the
    * @f$x@f$-axis, and `YLabel` a short description of the @f$y@f$-axis.
    * If `data` is a @f$n@f$-row matrix, then the first row
    * <tt>data</tt>@f$[0]@f$ represents the @f$x@f$-coordinate vector, and
    * every other row <tt>data</tt>@f$[i], i=1,…, n-1@f$, represents a
    * @f$y@f$-coordinate set of points. Therefore matrix
    * <tt>data</tt>@f$[i][ ]@f$, @f$i=0,…, n-1@f$, corresponds to
    * @f$n-1@f$ curves, all with the same @f$x@f$-coordinates. However,
    * only *the first* `numPoints` of `data` will be considered to plot
    * each curve.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series of point sets.
    *  @param numPoints    Number of points to plot
    */
   public XYLineChart (String title, String XLabel, String YLabel,
                       double[][] data, int numPoints) {
      super();
      dataset = new XYListSeriesCollection(data, numPoints);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `XYLineChart` instance using subsets of `data`.
    * `data[x][.]` will form the @f$x@f$-coordinates and `data[y][.]` will
    * form the @f$y@f$-coordinates of the chart. `title` sets a title,
    * `XLabel` is a short description of the @f$x@f$-axis, and `YLabel` is
    * a short description of the @f$y@f$-axis. Warning: if the new
    * @f$x@f$-axis coordinates are not monotone increasing, then they will
    * automatically be sorted in increasing order so the points will be
    * reordered, but the original `data` is not changed.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series of point sets.
    *  @param x            Index of data forming the @f$x@f$-coordinates
    *  @param y            Index of data forming the @f$y@f$-coordinates
    */
   public XYLineChart (String title, String XLabel, String YLabel,
                       double[][] data, int x, int y) {
      super();
      int len = data[0].length;
      double[][] proj = new double[2][len];
      for (int i = 0; i < len; i++) {
         proj[0][i] = data[x][i];
         proj[1][i] = data[y][i];
      }
      dataset = new XYListSeriesCollection(proj);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `XYLineChart` instance with data `data`. The input
    * parameter `data` represents a set of plotting data. A
    * DoubleArrayList from the Colt library is used to store the data. The
    * description is similar to the constructor  @ref YListChart with
    * `double[]... data`.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series of point sets.
    */
   public XYLineChart (String title, String XLabel, String YLabel,
                       DoubleArrayList... data) {
      super();
      dataset = new XYListSeriesCollection(data);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `XYLineChart` instance with data `data`. The input
    * parameter `data` represents a set of plotting data.
    * @ref org.jfree.data.xy.XYSeriesCollection is a `JFreeChart`
    * container class to store @f$XY@f$ plots.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series collection.
    */
   public XYLineChart (String title, String XLabel, String YLabel,
                       XYSeriesCollection data) {
      super();
      dataset = new XYListSeriesCollection(data);
      init (title, XLabel, YLabel);
   }

   /**
    * Adds a data series into the series collection. Vector `x` represents
    * the @f$x@f$-coordinates and vector `y` represents the
    * @f$y@f$-coordinates of the series. `name` and `plotStyle` are the
    * name and the plot style associated to the series.
    *  @param x            @f$x_i@f$ coordinates.
    *  @param y            @f$y_i@f$ coordinates.
    *  @param name         Name of the series.
    *  @param plotStyle    Plot style of the series.
    *  @return Integer that represent the new point set’s position in the
    * JFreeChart `XYSeriesCollection` object.
    */
   public int add (double[] x, double[] y, String name, String plotStyle) {
      int seriesIndex = add(x,y);
      getSeriesCollection().setName(seriesIndex, name);
      getSeriesCollection().setPlotStyle (seriesIndex, plotStyle);
      return seriesIndex;
   }

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
      int seriesIndex = getSeriesCollection().add(x,y);
      initAxis();
      return seriesIndex;
   }

   /**
    * Adds a data series into the series collection. Vector `x` represents
    * the @f$x@f$-coordinates and vector `y` represents the
    * @f$y@f$-coordinates of the series. Only *the first* `numPoints` of
    * `x` and `y` will be taken into account for the new series.
    *  @param x            @f$x_i@f$ coordinates.
    *  @param y            @f$y_i@f$ coordinates.
    *  @param numPoints    Number of points to add
    *  @return Integer that represent the new point set’s position in the
    * JFreeChart `XYSeriesCollection` object.
    */
   public int add (double[] x, double[] y, int numPoints) {
      int seriesIndex = getSeriesCollection().add(x, y, numPoints);
      initAxis();
      return seriesIndex;
   }

   /**
    * Adds the new collection of data series `data` into the series
    * collection. If `data` is a @f$n@f$-row matrix, then the first row
    * <tt>data</tt>@f$[0]@f$ represents the @f$x@f$-coordinate vector, and
    * every other row <tt>data</tt>@f$[i], i=1,…, n-1@f$, represents a
    * @f$y@f$-coordinate set of points. Therefore matrix
    * <tt>data</tt>@f$[i][ ]@f$, @f$i=0,…, n-1@f$, corresponds to
    * @f$n-1@f$ curves, all with the same @f$x@f$-coordinates.
    *  @param data         series of point sets.
    */
   public int add (double[][] data) {
      int seriesIndex = getSeriesCollection().add(data);
      initAxis();
      return seriesIndex;
   }

   /**
    * Adds the new collection of data series `data` into the series
    * collection. If `data` is a @f$n@f$-row matrix, then the first row
    * <tt>data</tt>@f$[0]@f$ represents the @f$x@f$-coordinate vector, and
    * every other row <tt>data</tt>@f$[i], i=1,…, n-1@f$, represents a
    * @f$y@f$-coordinate set of points. Therefore matrix
    * <tt>data</tt>@f$[i][ ]@f$, @f$i=0,…, n-1@f$, corresponds to
    * @f$n-1@f$ curves, all with the same @f$x@f$-coordinates. However,
    * only *the first* `numPoints` of `data` will be taken into account
    * for the new series.
    *  @param data         series of point sets.
    *  @param numPoints    Number of points to plot
    */
   public int add (double[][] data, int numPoints) {
      int seriesIndex = getSeriesCollection().add(data, numPoints);
      initAxis();
      return seriesIndex;
   }

   /**
    * Returns the chart’s dataset.
    *  @return the chart’s dataset.
    */
   public XYListSeriesCollection getSeriesCollection() {
      return (XYListSeriesCollection)dataset;
   }

   /**
    * Links a new dataset to the current chart.
    *  @param dataset      new dataset.
    */
   public void setSeriesCollection (XYListSeriesCollection dataset) {
      this.dataset = dataset;
   }

   /**
    * Synchronizes @f$X@f$-axis ticks to the @f$s@f$-th series
    * @f$x@f$-values.
    *  @param s            series used to define ticks.
    */
   public void setTicksSynchro (int s) {
      XYSeriesCollection seriesCollection =
          (XYSeriesCollection)this.dataset.getSeriesCollection();
      double[] values = new double[seriesCollection.getItemCount(s)];

      for(int i = 0; i < seriesCollection.getItemCount(s); i++)
         values[i] = seriesCollection.getXValue(s, i);

      XAxis.setLabels(values);
   }

   /**
    * Displays chart on the screen using Swing. This method creates an
    * application containing a chart panel displaying the chart. The
    * created frame is positioned on-screen, and displayed before it is
    * returned. The `width` and the `height` of the chart are measured in
    * pixels.
    *  @param width        frame width in pixels.
    *  @param height       frame height in pixels.
    *  @return frame containing the chart.
    *
    * ;
    */
   public JFrame view (int width, int height) {
      JFrame myFrame;
      if(chart.getTitle() != null)
         myFrame = new JFrame("XYLineChart from SSJ: " + chart.getTitle().getText());
      else
         myFrame = new JFrame("XYLineChart from SSJ");
      ChartPanel chartPanel = new ChartPanel(chart);
      chartPanel.setPreferredSize(new java.awt.Dimension(width, height));
      myFrame.setContentPane(chartPanel);
      myFrame.pack();
      myFrame.setDefaultCloseOperation (JFrame.DISPOSE_ON_CLOSE);
      myFrame.setLocationRelativeTo (null);
      myFrame.setVisible(true);
      return myFrame;
   }

   /**
    * Displays bar chart on the screen using Swing. This method creates an
    * application containing a bar chart panel displaying the chart. The
    * created frame is positioned on-screen, and displayed before it is
    * returned. The `width` and the `height` of the chart are measured in
    * pixels.
    *  @param width        frame width in pixels.
    *  @param height       frame height in pixels.
    *  @return frame containing the bar chart.
    *
    * ;
    */
   public JFrame viewBar (int width, int height) {
      JFrame myFrame;
      if (chart.getTitle() != null)
         myFrame = new JFrame("XYLineChart from SSJ: " + chart.getTitle().getText());
      else
         myFrame = new JFrame("XYLineChart from SSJ");

      XYPlot plot = (XYPlot) chart.getPlot();

      //Create the bar
      plot.setDataset(0, dataset.getSeriesCollection());
      final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(false, true);
      renderer.setSeriesPaint(0,Color.ORANGE);
      renderer.setSeriesShape(0,new Line2D.Double(0, 0, 0 , 1000));
      plot.setRenderer(0, renderer);

      //Create the points
      plot.setDataset(1, dataset.getSeriesCollection());
      final XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer(false, true);
      renderer2.setSeriesPaint(0,Color.ORANGE);
      renderer2.setSeriesShape(0,new Ellipse2D.Double(-2.0,-2.0,4.0,4.0));
      plot.setRenderer(1, renderer2);

      ChartPanel chartPanel = new ChartPanel(chart);
      chartPanel.setPreferredSize(new java.awt.Dimension(width, height));
      myFrame.setContentPane(chartPanel);
      myFrame.pack();
      myFrame.setDefaultCloseOperation (JFrame.DISPOSE_ON_CLOSE);
      myFrame.setLocationRelativeTo (null);
      myFrame.setVisible(true);
      return myFrame;
   }

   /**
    * @name Latex-specific method
    * @{
    */
   public String toLatex (double width, double height) {
      double xunit=0, yunit=0;
      double[] save = new double[4];

      if (dataset.getSeriesCollection().getSeriesCount() == 0)
         throw new IllegalArgumentException("Empty chart");

      //Calcul des parametres d'echelle et de decalage
      double XScale = computeXScale(XAxis.getTwinAxisPosition());
      double YScale = computeYScale(YAxis.getTwinAxisPosition());

      // taille d'une unite en x et en cm dans l'objet "tikzpicture"
      xunit = width / ((Math.max(XAxis.getAxis().getRange().getUpperBound(), XAxis.getTwinAxisPosition()) * XScale) - (Math.min(XAxis.getAxis().getRange().getLowerBound(), XAxis.getTwinAxisPosition()) * XScale));
      // taille d'une unite en y et en cm dans l'objet "tikzpicture"
      yunit = height / ((Math.max(YAxis.getAxis().getRange().getUpperBound(), YAxis.getTwinAxisPosition()) * YScale) - (Math.min(YAxis.getAxis().getRange().getLowerBound(), YAxis.getTwinAxisPosition()) * YScale));

      Formatter formatter = new Formatter(Locale.US);

      /*Entete du document*/
      if (latexDocFlag) {
         formatter.format("\\documentclass[12pt]{article}%n%n");
         formatter.format("\\usepackage{tikz}%n\\usetikzlibrary{plotmarks}%n\\begin{document}%n%n");
      }
      if (chart.getTitle() != null)
         formatter.format("%% PGF/TikZ picture from SSJ: %s%n", chart.getTitle().getText());
      else
         formatter.format("%% PGF/TikZ picture from SSJ %n");
      formatter.format("%% XScale = %s,  YScale = %s,  XShift = %s,  YShift = %s%n", XScale, YScale, XAxis.getTwinAxisPosition(), YAxis.getTwinAxisPosition());
      formatter.format("%% Therefore, thisFileXValue = (originalSeriesXValue+XShift)*XScale%n");
      formatter.format("%%        and thisFileYValue = (originalSeriesYValue+YShift)*YScale%n%n");
      if (chart.getTitle() != null)
         formatter.format("\\begin{figure}%n");
      formatter.format("\\begin{center}%n");
      formatter.format("\\begin{tikzpicture}[x=%scm, y=%scm]%n", xunit, yunit);
      formatter.format("\\footnotesize%n");
      if (grid)
         formatter.format("\\draw[color=lightgray] (%s, %s) grid[xstep = %s, ystep=%s] (%s, %s);%n",
            (Math.min(XAxis.getAxis().getRange().getLowerBound(), XAxis.getTwinAxisPosition())-XAxis.getTwinAxisPosition()) * XScale,
            (Math.min(YAxis.getAxis().getRange().getLowerBound(), YAxis.getTwinAxisPosition())-YAxis.getTwinAxisPosition()) * YScale,
            xstepGrid*XScale, ystepGrid*YScale,
            (Math.max(XAxis.getAxis().getRange().getUpperBound(), XAxis.getTwinAxisPosition())-XAxis.getTwinAxisPosition()) * XScale,
            (Math.max(YAxis.getAxis().getRange().getUpperBound(), YAxis.getTwinAxisPosition())-YAxis.getTwinAxisPosition()) * YScale );
      setTick0Flags();
      formatter.format("%s", XAxis.toLatex(XScale) );
      formatter.format("%s", YAxis.toLatex(YScale) );

      formatter.format("%s", dataset.toLatex(
         XScale, YScale,
         XAxis.getTwinAxisPosition(), YAxis.getTwinAxisPosition(),
         XAxis.getAxis().getLowerBound(), XAxis.getAxis().getUpperBound(),
         YAxis.getAxis().getLowerBound(), YAxis.getAxis().getUpperBound()));

      formatter.format("\\end{tikzpicture}%n");
      formatter.format("\\end{center}%n");
      if (chart.getTitle() != null) {
         formatter.format("\\caption{");
         formatter.format(chart.getTitle().getText());
         formatter.format("}%n\\end{figure}%n");
      }
      if (latexDocFlag)
         formatter.format("\\end{document}%n");
      return formatter.toString();
   }

}

/**
 * @}
 */