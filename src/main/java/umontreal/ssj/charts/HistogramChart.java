/*
 * Class:        HistogramChart
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

import   umontreal.ssj.stat.*;

import   org.jfree.chart.ChartFactory;
import   org.jfree.chart.ChartPanel;
import   org.jfree.chart.axis.NumberAxis;
import   org.jfree.chart.plot.XYPlot;
import   org.jfree.chart.plot.PlotOrientation;
import   org.jfree.data.statistics.HistogramBin;

import   cern.colt.list.DoubleArrayList;

import   java.util.ListIterator;
import   java.util.Locale;
import   java.util.Formatter;
import   javax.swing.JFrame;

/**
 * This class provides tools to create and manage histograms. The
 * @ref HistogramChart class is the simplest way to produce histograms. Each
 * @ref HistogramChart object is linked with an
 * @ref umontreal.ssj.charts.HistogramSeriesCollection dataset.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class HistogramChart extends XYChart {

   protected void init (String title, String XLabel, String YLabel) {
      // create the chart...
      chart = ChartFactory.createXYLineChart(
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
      //Initialize axis variables
      XAxis = new Axis( (NumberAxis)((XYPlot) chart.getPlot()).getDomainAxis(),
            Axis.ORIENTATION_HORIZONTAL );
      YAxis = new Axis( (NumberAxis)((XYPlot) chart.getPlot()).getRangeAxis() ,
            Axis.ORIENTATION_VERTICAL );
      setAutoRange(false, true, true, true);
   }

   /**
    * Initializes a new `HistogramChart` instance with an empty data set.
    */
   public HistogramChart () {
      super();
      dataset = new HistogramSeriesCollection();
      init (null, null, null);
   }

   /**
    * Initializes a new `HistogramChart` instance with input `data`.
    * `title` is a title, `XLabel` is a short description of the
    * @f$x@f$-axis, and `YLabel` a short description of the @f$y@f$-axis.
    * The input parameter `data` represents a collection of observation
    * sets. Therefore <tt>data</tt>@f$[i], i = 0,…,n-1@f$, is used to plot
    * the @f$i@f$th histogram.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series of point sets.
    */
   public HistogramChart (String title, String XLabel, String YLabel,
                          double[]... data) {
      super();
      dataset = new HistogramSeriesCollection(data);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `HistogramChart` instance with input `data`.
    * `title` is a title, `XLabel` is a short description of the
    * @f$x@f$-axis, and `YLabel` a short description of the @f$y@f$-axis.
    * The input parameter `data` represents an observation set. Only *the
    * first* `numPoints` of `data` will be considered to plot the
    * histogram.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series of point sets.
    *  @param numPoints    Number of points to plot
    */
   public HistogramChart (String title, String XLabel, String YLabel,
                          double[] data, int numPoints) {
      super();
      double[] datan = new double[numPoints];
      System.arraycopy (data, 0, datan, 0, numPoints);
      dataset = new HistogramSeriesCollection(datan);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `HistogramChart` instance with data `data`. Each
    * `DoubleArrayList` input parameter represents a collection of
    * observation sets.  DoubleArrayList is from the Colt library and is
    * used to store data.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series of observation sets.
    */
   public HistogramChart (String title, String XLabel, String YLabel,
                          DoubleArrayList... data) {
      super();
      dataset = new HistogramSeriesCollection(data);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `HistogramChart` instance with data arrays
    * contained in each  @ref umontreal.ssj.stat.TallyStore object. The
    * input parameter `tallies` represents a collection of observation
    * sets.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param tallies      series of observation sets.
    */
   public HistogramChart (String title, String XLabel, String YLabel,
                          TallyStore... tallies) {
      super();
      dataset = new HistogramSeriesCollection(tallies);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `HistogramChart` instance with data `data`. The
    * input parameter `data` represents a set of plotting data.
    * @ref umontreal.ssj.charts.CustomHistogramDataset is a
    * <tt>JFreeChart</tt>-like container class that stores and manages
    * observation sets.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series collection.
    */
   public HistogramChart (String title, String XLabel, String YLabel,
                          CustomHistogramDataset data) {
      super();
      dataset = new HistogramSeriesCollection(data);
      init (title, XLabel, YLabel);
   }

    /**
     * Initializes a new `HistogramChart` instance with data `count` and
     * `bound`. The adjacent categories (or bins) are specified as
     * non-overlapping intervals: bin[j] contains the values in the
     * interval [<tt>bound[j]</tt>, <tt>bound[j+1]</tt>], and `count[j]`
     * is the number of such values. Thus the length of `bound` must be
     * equal to the length of `count` plus one: the last value of `bound`
     * is the right boundary of the last bin.
     *  @param title        chart title.
     *  @param XLabel       Label on @f$x@f$-axis.
     *  @param YLabel       Label on @f$y@f$-axis.
     *  @param count        the number of observation between each bound.
     *  @param bound        the bounds of the observations
     */
    public HistogramChart (String title, String XLabel, String YLabel,
                           int[] count, double[] bound) {
      super();
      if (bound.length != count.length + 1)
         throw new IllegalArgumentException (
            "bound.length must be equal to count.length + 1");
      final int nb = count.length;
      int sum = 0;
      for (int i = 0 ; i < nb; i++) sum +=count[i];
      double[] data = new double [sum];

      int k = 0;
      double h;
      for (int i = 0 ; i < nb; i++) {
         h = bound[i + 1] - bound[i];
         if (count[i] > 0)
            h /= count[i];
         if (i == nb - 1) {
            for (int j = 0 ; j < count[i] ; j++)
               data[k++] = bound[i + 1] - j*h;
         } else {
            for (int j = 0 ; j < count[i] ; j++)
               data[k++] = bound[i] + j*h;
         }
      }

      dataset = new HistogramSeriesCollection(data, sum);
      init (title, XLabel, YLabel);
      ((HistogramSeriesCollection) dataset).setBins(0, nb);
   }

    /**
     * Initializes a new `HistogramChart` instance with data arrays
     * contained in each  @ref umontreal.ssj.stat.TallyHistogram object.
     * The input parameter `tallies` represents a collection of
     * observation sets. The 2 extra bins at the beginning and at the end
     * of the tallies are not counted nor represented in the chart.
     *  @param title        chart title.
     *  @param XLabel       Label on @f$x@f$-axis.
     *  @param YLabel       Label on @f$y@f$-axis.
     *  @param tallies      series of observation sets.
     */
    public HistogramChart (String title, String XLabel, String YLabel,
                           TallyHistogram... tallies) {
      super();
      dataset = new HistogramSeriesCollection(tallies);
      init (title, XLabel, YLabel);
   }

   public void setAutoRange (boolean right, boolean top)  {
         throw new UnsupportedOperationException(
            "You can't use setAutoRange with HistogramChart class, use setAutoRange().");
   }
   public void setManuelRange (double [] range, boolean right, boolean top) {
         throw new UnsupportedOperationException(
            "You can't use setManuelRange with HistogramChart class, use setManuelRange(range).");
   }

/**
 * Returns the chart’s dataset.
 *  @return the chart’s dataset.
 */
public HistogramSeriesCollection getSeriesCollection() {
      return (HistogramSeriesCollection)dataset;
   }

   /**
    * Links a new dataset to the current chart.
    *  @param dataset      new dataset.
    */
   public void setSeriesCollection (HistogramSeriesCollection dataset) {
      this.dataset = dataset;
   }

   /**
    * Synchronizes @f$x@f$-axis ticks to the @f$s@f$-th histogram bins if
    * the number of bins is not larger than 10; otherwise, choose
    * approximately 10 ticks.
    *  @param s            selects histogram used to define ticks.
    */
   public void setTicksSynchro (int s) {
      if (((CustomHistogramDataset)this.dataset.getSeriesCollection()).getBinWidth(s) == -1){
         DoubleArrayList newTicks = new DoubleArrayList();
         ListIterator binsIter = ((HistogramSeriesCollection)this.dataset).getBins(s).listIterator();

         int i = 0;
         HistogramBin prec = (HistogramBin)binsIter.next();
         double var;
         newTicks.add(prec.getStartBoundary());
         newTicks.add(var = prec.getEndBoundary());
         HistogramBin temp;
         while(binsIter.hasNext()) {
            temp = (HistogramBin)binsIter.next();
            if(temp.getStartBoundary() != var) {
               newTicks.add(var = temp.getStartBoundary());
            } else if(temp.getEndBoundary() != var) {
               newTicks.add(var = temp.getEndBoundary());
            }
         }
         XAxis.setLabels(newTicks.elements());
      }
      else {
         // set a label-tick for each bin, if num bins is <= 10
         int n = ((HistogramSeriesCollection)this.dataset).getBins(s).size();
         if (n > 10) {
            // number of bins is too large, set ~10 labels-ticks for histogram
            n = 10;
            double[] B = ((HistogramSeriesCollection)this.dataset).getDomainBounds();
            double w = (B[1] - B[0]) / n;
            XAxis.setLabels(w);
         } else {
            XAxis.setLabels(((CustomHistogramDataset)this.dataset.getSeriesCollection()).getBinWidth(s));
         }
      }
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
         myFrame = new JFrame("HistogramChart from SSJ: " + chart.getTitle().getText());
      else
         myFrame = new JFrame("HistogramChart from SSJ");
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

      if (dataset.getSeriesCollection().getSeriesCount() == 0)
         throw new IllegalArgumentException("Empty chart");
      if (YAxis.getTwinAxisPosition() < 0)
         YAxis.setTwinAxisPosition(0);

      // Calcul des parametres d'echelle et de decalage
      double XScale = computeXScale(XAxis.getTwinAxisPosition());
      double YScale = computeYScale(YAxis.getTwinAxisPosition());

      // taille d'une unite en x et en cm dans l'objet "tikzpicture"
      xunit = width / ( (Math.max(XAxis.getAxis().getRange().getUpperBound(), XAxis.getTwinAxisPosition()) * XScale) - (Math.min(XAxis.getAxis().getRange().getLowerBound(), XAxis.getTwinAxisPosition()) * XScale) );
      // taille d'une unite en y et en cm dans l'objet "tikzpicture"
     yunit = height / ( (Math.max(YAxis.getAxis().getRange().getUpperBound(), YAxis.getTwinAxisPosition()) * YScale) - (Math.min(YAxis.getAxis().getRange().getLowerBound(), YAxis.getTwinAxisPosition()) * YScale) );

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