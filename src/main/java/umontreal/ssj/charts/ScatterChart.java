/*
 * Class:        ScatterChart
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

import   org.jfree.chart.axis.NumberAxis;
import   org.jfree.chart.ChartFactory;
import   org.jfree.chart.ChartPanel;
import   org.jfree.chart.plot.XYPlot;
import   org.jfree.chart.plot.PlotOrientation;
import   org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import   org.jfree.chart.renderer.xy.XYDotRenderer;
import   org.jfree.data.xy.XYSeriesCollection;
import   java.util.Locale;
import   java.util.Formatter;
import   cern.colt.list.DoubleArrayList;
import   javax.swing.JFrame;

/**
 * This class provides tools to create and manage scatter plots. Using the
 * @ref ScatterChart class is the simplest way to produce scatter plots only.
 * Each  @ref ScatterChart object is linked with a
 * @ref umontreal.ssj.charts.XYListSeriesCollection data set.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class ScatterChart extends XYChart {

   protected void init (String title, String XLabel, String YLabel) {
      // create the chart...
      chart = ChartFactory.createScatterPlot (
         title,                    // chart title
         XLabel,                   // x axis label
         YLabel,                   // y axis label
         dataset.getSeriesCollection(), // data
         PlotOrientation.VERTICAL,
         true,                     // include legend
         true,                     // tooltips
         false                     // urls
      );

      ((XYPlot)chart.getPlot()).setRenderer(dataset.getRenderer());
      // Initialize axis variables
      initAxis();

      int nb = getSeriesCollection().getSeriesCollection().getSeriesCount();
      for (int i = 0 ; i < nb ; i++) {
         getSeriesCollection().setDashPattern(i, "only marks");
         getSeriesCollection().setMarksType(i, "+");
      }
   }

   protected void initAxis(){
      XAxis = new Axis((NumberAxis)((XYPlot)chart.getPlot()).getDomainAxis(),
                        Axis.ORIENTATION_HORIZONTAL);
      YAxis = new Axis((NumberAxis)((XYPlot)chart.getPlot()).getRangeAxis(),
                        Axis.ORIENTATION_VERTICAL);
      setAutoRange(true, true);
   }

   /**
    * Initializes a new `ScatterChart` instance with an empty data set.
    */
   public ScatterChart() {
      super();
      dataset = new XYListSeriesCollection();
      init (null, null, null);
   }

   /**
    * Initializes a new `ScatterChart` instance with data `data`. `title`
    * is a title, `XLabel` is a short description of the @f$x@f$-axis and
    * `YLabel` a short description of the @f$y@f$-axis. The input
    * parameter `data` represents sets of plotting data. For example, if
    * one @f$n@f$-row matrix `data1` is given as argument `data`, then the
    * first row <tt>data1</tt>@f$[0]@f$ represents the @f$x@f$-coordinate
    * vector, and every other row <tt>data1</tt>@f$[i], i=1,…, n-1@f$,
    * represents the @f$y@f$-coordinates of a set of points. Therefore
    * matrix `data1` corresponds to @f$n-1@f$ sets of points, all with the
    * same @f$x@f$-coordinates. However, one may want to plot sets of
    * points with different @f$x@f$-coordinates. In that case, one should
    * give the points as matrices with two rows. For examples, if the
    * argument `data` is made of three 2-row matrices `data1`, `data2` and
    * `data3`, then they represents three different sets of points,
    * <tt>data*</tt>@f$[0]@f$ giving the @f$x@f$-coordinates, and
    * <tt>data*</tt>@f$[1]@f$ the @f$y@f$-coordinates of the points.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series of point sets.
    */
   public ScatterChart (String title, String XLabel, String YLabel,
                        double[][]... data) {
      super();
      dataset = new XYListSeriesCollection(data);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `ScatterChart` instance with sets of points
    * `data`. `title` is a title, `XLabel` is a short description of the
    * @f$x@f$-axis, and `YLabel` a short description of the @f$y@f$-axis.
    * If `data` is a @f$n@f$-row matrix, then the first row
    * <tt>data</tt>@f$[0]@f$ represents the @f$x@f$-coordinate vector, and
    * every other row <tt>data</tt>@f$[i], i=1,…, n-1@f$, represents a
    * @f$y@f$-coordinate vector. Therefore matrix <tt>data</tt>@f$[i][
    * ]@f$, @f$i=0,…, n-1@f$, corresponds to @f$n-1@f$ sets of points, all
    * with the same @f$x@f$-coordinates. However, only *the first*
    * `numPoints` of each set <tt>data</tt>@f$[i]@f$ (i.e. the first
    * `numPoints` columns of each row) will be plotted.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series of point sets.
    *  @param numPoints    Number of points to plot
    */
   public ScatterChart (String title, String XLabel, String YLabel,
                        double[][] data, int numPoints) {
      super();
      dataset = new XYListSeriesCollection(data, numPoints);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `ScatterChart` instance using subsets of `data`.
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
   public ScatterChart (String title, String XLabel, String YLabel,
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
    * Initializes a new `ScatterChart` instance with data `data`. The
    * input parameter `data` represents a set of plotting data. A
    * DoubleArrayList from the Colt library is used to store the data. The
    * description is similar to the above constructor with `double[]...
    * data`.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series of point sets.
    */
   public ScatterChart (String title, String XLabel, String YLabel,
                        DoubleArrayList... data) {
      super();
      dataset = new XYListSeriesCollection(data);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `ScatterChart` instance with data `data`. The
    * input parameter `data` represents a set of plotting data.
    * @ref org.jfree.data.xy.XYSeriesCollection is a `JFreeChart`
    * container class to store @f$XY@f$ plots.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series collection.
    */
   public ScatterChart (String title, String XLabel, String YLabel,
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
      return add (x, y, x.length);
   }

   /**
    * Adds a data series into the series collection. Vector `x` represents
    * the @f$x@f$-coordinates and vector `y` represents the
    * @f$y@f$-coordinates of the series. Only *the first* `numPoints` of
    * `x` and `y` will be taken into account for the new series.
    *  @param x            @f$x_i@f$ coordinates.
    *  @param y            @f$y_i@f$ coordinates.
    *  @param numPoints    Number of points to add.
    *  @return Integer that represent the new point set’s position in the
    * JFreeChart `XYSeriesCollection` object.
    */
   public int add (double[] x, double[] y, int numPoints) {
      int seriesIndex = getSeriesCollection().add(x, y, numPoints);
      initAxis();
      getSeriesCollection().setMarksType(seriesIndex, "+");
      getSeriesCollection().setDashPattern(seriesIndex, "only marks");
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
      if (chart.getTitle () != null)
         myFrame = new JFrame ("ScatterChart from SSJ: " + chart.getTitle ().getText ());
      else
         myFrame = new JFrame ("ScatterChart from SSJ");
      XYPlot plot = chart.getXYPlot ();

/*    // The drawn points are somewhat big, of different shapes, unfilled
      XYLineAndShapeRenderer shape = new XYLineAndShapeRenderer(false, true);
      int nb = getSeriesCollection().getSeriesCollection().getSeriesCount();
      for (int i = 0 ; i < nb ; i++) {
         shape.setSeriesShapesFilled(i, false);
         plot.setRenderer(i, shape);
      }
*/
      // The drawn points are all square, filled
      XYDotRenderer shape = new XYDotRenderer();
      final int dotSize = 3;
      shape.setDotWidth(dotSize);
      shape.setDotHeight(dotSize);
      int nb = getSeriesCollection().getSeriesCollection().getSeriesCount();
      for (int i = 0 ; i < nb ; i++)
         plot.setRenderer(i, shape);

      ChartPanel chartPanel = new ChartPanel (chart);
      chartPanel.setPreferredSize (new java.awt.Dimension(width, height));
      myFrame.setContentPane (chartPanel);
      myFrame.pack();
      myFrame.setDefaultCloseOperation (JFrame.DISPOSE_ON_CLOSE);
      myFrame.setLocationRelativeTo (null);
      myFrame.setVisible (true);
      return myFrame;
   }

   /**
    * @name LaTex-specific method
    * @{
    */
   public String toLatex (double width, double height) {
      double xunit=0, yunit=0;
      double[] save = new double[4];

      if(dataset.getSeriesCollection().getSeriesCount() == 0)
         throw new IllegalArgumentException("Empty chart");

      //Calcul des parametres d'echelle et de decalage
      double XScale = computeXScale(XAxis.getTwinAxisPosition());
      double YScale = computeYScale(YAxis.getTwinAxisPosition());

         //taille d'une unite en x et en cm dans l'objet "tikzpicture"
      xunit = width / ((Math.max(XAxis.getAxis().getRange().getUpperBound(),
                            XAxis.getTwinAxisPosition()) * XScale)
             - (Math.min(XAxis.getAxis().getRange().getLowerBound(),
                    XAxis.getTwinAxisPosition()) * XScale));
         //taille d'une unite en y et en cm dans l'objet "tikzpicture"
      yunit = height / ((Math.max(YAxis.getAxis().getRange().getUpperBound(),
                         YAxis.getTwinAxisPosition()) * YScale)
           - (Math.min(YAxis.getAxis().getRange().getLowerBound(),
                  YAxis.getTwinAxisPosition()) * YScale));

      Formatter formatter = new Formatter(Locale.US);

      /*Entete du document*/
      if (latexDocFlag) {
         formatter.format("\\documentclass[12pt]{article}%n%n");
         formatter.format("\\usepackage{tikz}%n\\usetikzlibrary{plotmarks}%n\\begin{document}%n%n");
      }
      if(chart.getTitle() != null)
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
      if(grid)
         formatter.format("\\draw[color=lightgray] (%s, %s) grid[xstep = %s, ystep=%s] (%s, %s);%n",
            (Math.min(XAxis.getAxis().getRange().getLowerBound(),
              XAxis.getTwinAxisPosition())-XAxis.getTwinAxisPosition()) * XScale,
            (Math.min(YAxis.getAxis().getRange().getLowerBound(),
             YAxis.getTwinAxisPosition())-YAxis.getTwinAxisPosition()) * YScale,
            xstepGrid*XScale, ystepGrid*YScale,
            (Math.max(XAxis.getAxis().getRange().getUpperBound(),
              XAxis.getTwinAxisPosition())-XAxis.getTwinAxisPosition()) * XScale,
            (Math.max(YAxis.getAxis().getRange().getUpperBound(),
              YAxis.getTwinAxisPosition())-YAxis.getTwinAxisPosition()) * YScale );
      setTick0Flags();
      formatter.format("%s", XAxis.toLatex(XScale) );
      formatter.format("%s", YAxis.toLatex(YScale) );

      formatter.format("%s", dataset.toLatex(XScale, YScale,
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