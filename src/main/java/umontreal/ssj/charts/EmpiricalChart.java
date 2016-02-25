/*
 * Class:        EmpiricalChart
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

import   umontreal.ssj.stat.TallyStore;

import   org.jfree.chart.ChartPanel;
import   org.jfree.chart.ChartFactory;
import   org.jfree.chart.axis.NumberAxis;
import   org.jfree.chart.plot.XYPlot;
import   org.jfree.chart.plot.PlotOrientation;
import   org.jfree.data.xy.XYSeriesCollection;
import   org.jfree.data.xy.XYSeries;

import   cern.colt.list.DoubleArrayList;

import   java.util.ListIterator;
import   java.util.Locale;
import   java.util.Formatter;
import   javax.swing.JFrame;

/**
 * This class provides additional tools to create and manage empirical plots.
 * Empirical plots are used to plot empirical distributions. The
 * @ref EmpiricalChart class is the simplest way to produce empirical plots
 * only. Each  @ref EmpiricalChart object is linked with an
 * @ref umontreal.ssj.charts.EmpiricalSeriesCollection data set.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class EmpiricalChart extends XYChart {

   protected void init (String title, String XLabel, String YLabel) {
      // create the chart...
      chart = ChartFactory.createXYLineChart(
         title,                    // chart title
         XLabel,                   // x axis label
         YLabel,                   // y axis label
         dataset.getSeriesCollection(), // data
         PlotOrientation.VERTICAL,
         true,                     // include legend
         true,                     // tool tips
         false                     // urls
      );
      ((XYPlot)chart.getPlot()).setRenderer(dataset.getRenderer());
      // Initialize axis variables
      XAxis = new Axis((NumberAxis)((XYPlot) chart.getPlot()).getDomainAxis(),
                Axis.ORIENTATION_HORIZONTAL);
      YAxis = new Axis((NumberAxis)((XYPlot) chart.getPlot()).getRangeAxis(),
                Axis.ORIENTATION_VERTICAL);
      fixZeroPoint();
   }


   private void fixZeroPoint() {
      // reset the first point (x0, 0) with x0 at the beginning of x-axis
      double xmin = Math.min (XAxis.getAxis().getRange().getLowerBound(),
                                                   XAxis.getTwinAxisPosition());
      XYSeriesCollection col = (XYSeriesCollection)dataset.getSeriesCollection();
      for (int i = 0; i < col.getSeriesCount(); i++) {
         XYSeries ser = col.getSeries (i);
         ser.remove(0);   // remove temporary 0-point
         ser.add(xmin, 0); // replace
      }
   }

   /**
    * Initializes a new `EmpiricalChart` instance with an empty data set.
    */
   public EmpiricalChart() {
      super();
      dataset = new EmpiricalSeriesCollection();
      init (null, null, null);
   }

   /**
    * Initializes a new `EmpiricalChart` instance with data `data`.
    * `title` is a title, `XLabel` is a short description of the
    * @f$x@f$-axis and `YLabel` a short description of the @f$y@f$-axis.
    * The input vectors `data` represents a collection of observation
    * sets. Each vector of `data` represents a @f$x@f$-coordinates set.
    * Therefore <tt>data</tt>@f$[i], i = 0,…,n-1@f$, is used to draw the
    * @f$i@f$-th plot. The values of each observation set
    * <tt>data</tt>@f$[i]@f$ *must be sorted* in increasing order.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series of point sets.
    */
   public EmpiricalChart (String title, String XLabel, String YLabel,
                          double[]... data) {
      super();
      dataset = new EmpiricalSeriesCollection(data);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `EmpiricalChart` instance with a set of points
    * `data`. `title` is a title, `XLabel` is a short description of the
    * @f$x@f$-axis and `YLabel` a short description of the @f$y@f$-axis.
    * Vector `data` represents a @f$x@f$-coordinates set. The values of
    * this observation set *must be sorted* in increasing order. Only *the
    * first* `numPoints` of `data` will be considered to plot.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series of point sets.
    *  @param numPoints    number of points to plot
    */
   public EmpiricalChart (String title, String XLabel, String YLabel,
                          double[] data, int numPoints) {
      super();
      dataset = new EmpiricalSeriesCollection(data, numPoints);
      init (title, XLabel, YLabel);
   }

   /**
    * Similar to the above constructor, but with `DoubleArrayList`. A
    * DoubleArrayList from the Colt library is used to store data. The
    * values of each observation set <tt>data</tt>@f$[i]@f$ *must be
    * sorted* in increasing order.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series of point sets.
    */
   public EmpiricalChart (String title, String XLabel, String YLabel,
                          DoubleArrayList... data) {
      super();
      dataset = new EmpiricalSeriesCollection(data);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `EmpiricalChart` instance with data arrays
    * contained in each  @ref umontreal.ssj.stat.TallyStore object. The
    * input parameter `tallies` represents a collection of observation
    * sets. Therefore, the @f$i@f$-th `tallies` is used to draw the
    * @f$i@f$-th plot.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param tallies      series of observation sets.
    */
   public EmpiricalChart (String title, String XLabel, String YLabel,
                          TallyStore... tallies) {
      super();
      dataset = new EmpiricalSeriesCollection(tallies);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `EmpiricalChart` instance with data `data`. The
    * input parameter `data` represents a set of plotting data.
    * `XYSeriesCollection` is a <tt>JFreeChart</tt>-like container class
    * used to store and manage observation sets.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series collection.
    */
   public EmpiricalChart (String title, String XLabel, String YLabel,
                          XYSeriesCollection data) {
      super();
      dataset = new EmpiricalSeriesCollection(data);
      init (title, XLabel, YLabel);
   }

   /**
    * Returns the chart’s dataset.
    *  @return the chart’s dataset.
    */
   public EmpiricalSeriesCollection getSeriesCollection() {
      return (EmpiricalSeriesCollection)dataset;
   }

   /**
    * Links a new dataset to the current chart.
    *  @param dataset      new dataset.
    */
   public void setSeriesCollection (EmpiricalSeriesCollection dataset) {
      this.dataset = dataset;
   }

   /**
    * Synchronizes @f$x@f$-axis ticks to the @f$s@f$-th series
    * @f$x@f$-values.
    *  @param s            series used to define ticks.
    */
   public void setTicksSynchro (int s) {
      XYSeriesCollection seriesCollection = (XYSeriesCollection)this.dataset.getSeriesCollection();
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
    */
   public JFrame view (int width, int height) {
      JFrame myFrame;
      if(chart.getTitle() != null)
         myFrame = new JFrame("EmpiricalChart from SSJ : " + chart.getTitle().getText());
      else
         myFrame = new JFrame("EmpiricalChart from SSJ");
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
    * @name LaTeX-specific method
    * @{
    */
   public String toLatex (double width, double height) {
      double xunit, yunit;
      double[] save = new double[4];

      if(dataset.getSeriesCollection().getSeriesCount() == 0)
         throw new IllegalArgumentException("Empty chart");

      //Calcul des parametres d'echelle et de decalage
      double XScale = computeXScale(XAxis.getTwinAxisPosition());
      double YScale = computeYScale(YAxis.getTwinAxisPosition());

      xunit = width / ((Math.max(XAxis.getAxis().getRange().getUpperBound(), XAxis.getTwinAxisPosition()) * XScale) - (Math.min(XAxis.getAxis().getRange().getLowerBound(), XAxis.getTwinAxisPosition()) * XScale));
      //taille d'une unite en x et en cm dans l'objet "tikzpicture"
      yunit = height / ((Math.max(YAxis.getAxis().getRange().getUpperBound(), YAxis.getTwinAxisPosition()) * YScale) - (Math.min(YAxis.getAxis().getRange().getLowerBound(), YAxis.getTwinAxisPosition()) * YScale));
      //taille d'une unite en y et en cm dans l'objet "tikzpicture"

      Formatter formatter = new Formatter(Locale.US);

      /*Entete du document*/
      if (latexDocFlag) {
         formatter.format("\\documentclass[12pt]{article}%n%n");
         formatter.format("\\usepackage{tikz}%n\\usetikzlibrary{plotmarks}%n\\begin{document}%n%n");
      }
      if(chart.getTitle() != null)
         formatter.format("%% PGF/TikZ picture from SSJ : %s%n", chart.getTitle().getText());
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
            (Math.max(YAxis.getAxis().getRange().getUpperBound(), YAxis.getTwinAxisPosition())-YAxis.getTwinAxisPosition()) * YScale);
      setTick0Flags();
      formatter.format("%s", XAxis.toLatex(XScale));
      formatter.format("%s", YAxis.toLatex(YScale));

      formatter.format("%s", dataset.toLatex(XScale, YScale, XAxis.getTwinAxisPosition(), YAxis.getTwinAxisPosition(),
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