/*
 * Class:        BoxChart
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
import   org.jfree.chart.plot.CategoryPlot;
import   org.jfree.chart.renderer.category.BoxAndWhiskerRenderer;
import   org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;

import   java.util.Locale;
import   java.util.Formatter;
import   javax.swing.JFrame;

/**
 * This class provides tools to create and manage box-and-whisker plots. Each
 * @ref BoxChart object is linked with a
 * @ref umontreal.ssj.charts.BoxSeriesCollection data set.
 *
 * A boxplot is a convenient way of viewing sets of numerical data through
 * their summaries: the smallest observation, first quartile (@f$Q_1 =
 * x_{.25}@f$), median (@f$Q_2 = x_{.5}@f$), third quartile (@f$Q_3 =
 * x_{.75}@f$), and largest observation. Sometimes, the mean and the outliers
 * are also plotted.
 *
 * In the charts created by this class, the box has its lower limit at
 * @f$Q_1@f$ and its upper limit at @f$Q_3@f$. The median is indicated by the
 * line inside the box, while the mean is at the center of the filled circle
 * inside the box. Define the interquartile range as (@f$Q_3 - Q_1@f$). Any
 * data observation which is more than @f$1.5(Q_3 - Q_1)@f$ lower than the
 * first quartile or @f$1.5(Q_3 - Q_1)@f$ higher than the third quartile is
 * considered an outlier. The smallest and the largest values that are not
 * outliers are connected to the box with a vertical line or "whisker" which
 * is ended by a horizontal line. Outliers are indicated by hollow circles
 * outside the whiskers. Triangles indicate the existence of very far
 * outliers.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class BoxChart extends CategoryChart {

   protected void init (String title, String XLabel, String YLabel) {
      // create the chart...
      chart = ChartFactory.createBoxAndWhiskerChart (
         title,                         // chart title
         XLabel,                        // x axis label
         YLabel,                        // y axis label
         (DefaultBoxAndWhiskerCategoryDataset)dataset.getSeriesCollection(), // data
         true                          // include legend
      );

      ((CategoryPlot)chart.getPlot()).setRenderer(dataset.getRenderer());
      // Initialize axis variables
      initAxis();
   }

   protected void initAxis(){
      YAxis = new Axis((NumberAxis)((CategoryPlot) chart.getPlot()).getRangeAxis(),
                        Axis.ORIENTATION_VERTICAL);
      setAutoRange();
   }

   /**
    * Initializes a new `BoxChart` instance with an empty data set.
    */
   public BoxChart() {
      super();
      dataset = new BoxSeriesCollection();
      init (null, null, null);
   }

   /**
    * Initializes a new `BoxChart` instance with data `data`. `title` is a
    * title, `XLabel` is a short description of the @f$x@f$-axis, and
    * `YLabel` a short description of the @f$y@f$-axis. The input
    * parameter `data` represents a set of plotting data. Only *the first*
    * `numPoints` of `data` will be considered for the plot.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         point sets.
    *  @param numPoints    Number of points to plot
    */
   public BoxChart (String title, String XLabel, String YLabel,
                    double[] data, int numPoints) {
      super();
      dataset = new BoxSeriesCollection(data, numPoints);
      init (title, XLabel, YLabel);
   }

   /**
    * Initializes a new `BoxChart` instance with data `data`. `title` sets
    * a title, `XLabel` is a short description of the @f$x@f$-axis, and
    * `YLabel` is a short description of the @f$y@f$-axis. The input
    * parameter `data` represents a set of plotting data.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    *  @param data         series of point sets.
    */
   public BoxChart (String title, String XLabel, String YLabel,
                    double[]... data) {
      super();
      dataset = new BoxSeriesCollection(data);
      init (title, XLabel, YLabel);
   }

   /**
    * Adds a data series into the series collection. Vector `data`
    * represents a set of plotting data.
    *  @param data         point sets.
    *  @return Integer that represent the new point set’s position in the
    * JFreeChart `BoxSeriesCollection` object.
    */
   public int add (double[] data) {
      return add(data, data.length);
   }

   /**
    * Adds a data series into the series collection. Vector `data`
    * represents a set of plotting data. Only *the first* `numPoints` of
    * `data` will be taken into account for the new series.
    *  @param data         point set.
    *  @param numPoints    number of points to add.
    *  @return Integer that represent the new point set’s position in the
    * JFreeChart `BoxSeriesCollection` object.
    */
   public int add (double[] data, int numPoints) {
      int seriesIndex = getSeriesCollection().add(data, numPoints);
      initAxis();
      return seriesIndex;
   }

   /**
    * Returns the chart’s dataset.
    *  @return the chart’s dataset.
    */
   public BoxSeriesCollection getSeriesCollection() {
      return (BoxSeriesCollection)dataset;
   }

   /**
    * Links a new dataset to the current chart.
    *  @param dataset      new dataset.
    */
   public void setSeriesCollection (BoxSeriesCollection dataset) {
      this.dataset = dataset;
   }

   /**
    * Sets `fill` to `true`, if the boxes are to be filled.
    *  @param fill         true if the boxes are filled
    */
   public void setFillBox (boolean fill) {
      ((BoxAndWhiskerRenderer)dataset.getRenderer()).setFillBox(fill);
   }

   /**
    * Displays chart on the screen using Swing. This method creates an
    * application containing a chart panel displaying the chart. The
    * created frame is positioned on-screen, and displayed before it is
    * returned. The circle represents the mean, the dark line inside the
    * box is the median, the box limits are the first and third quartiles,
    * the lower whisker (the lower line outside the box) is the first
    * decile, and the upper whisker is the ninth decile. The outliers, if
    * any, are represented by empty circles, or arrows if outside the
    * range bounds.
    *  @param width        frame width.
    *  @param height       frame height.
    *  @return frame containing the chart.
    *
    * ;
    */
   public JFrame view (int width, int height) {
      JFrame myFrame;
      if(chart.getTitle() != null)
         myFrame = new JFrame ("BoxChart from SSJ : " + chart.getTitle().getText());
      else
         myFrame = new JFrame ("BoxChart from SSJ");
      ChartPanel chartPanel = new ChartPanel (chart);
      chartPanel.setPreferredSize (new java.awt.Dimension(width, height));
      myFrame.setContentPane (chartPanel);
      myFrame.pack ();
      myFrame.setDefaultCloseOperation (JFrame.DISPOSE_ON_CLOSE);
      myFrame.setLocationRelativeTo (null);
      myFrame.setVisible (true);
      return myFrame;
   }

  /**
   * @name Latex-specific method
   * @{
   */

  /**
   * NOT IMPLEMENTED.
   *  @param width        Chart’s width in centimeters.
   *  @param height       Chart’s height in centimeters.
   *  @return LaTeX source code.
   */
  public String toLatex (double width, double height) {
  throw new UnsupportedOperationException(" NOT implemented yet");
/*
      double yunit;
      double[] save = new double[4];

      if(dataset.getSeriesCollection().getColumnCount() == 0)
         throw new IllegalArgumentException("Empty chart");

      //Calcul des parametres d'echelle et de decalage
      double YScale = computeYScale(YAxis.getTwinAxisPosition());


      yunit = height / ( (Math.max(YAxis.getAxis().getRange().getUpperBound(), YAxis.getTwinAxisPosition()) * YScale) - (Math.min(YAxis.getAxis().getRange().getLowerBound(), YAxis.getTwinAxisPosition()) * YScale) );
      //taille d'une unite en y et en cm dans l'objet "tikzpicture"

      Formatter formatter = new Formatter(Locale.US);

      // Entete du document
      formatter.format("\\documentclass[12pt]{article}%n%n");
      formatter.format("\\usepackage{tikz}%n\\usetikzlibrary{plotmarks}%n\\begin{document}%n%n");
      if(chart.getTitle() != null)
         formatter.format("%% PGF/TikZ picture from SSJ : %s%n", chart.getTitle().getText());
      else
         formatter.format("%% PGF/TikZ picture from SSJ %n");
      formatter.format("%%  YScale = %s, YShift = %s%n", YScale,  YAxis.getTwinAxisPosition());
      formatter.format("%%        and thisFileYValue = (originalSeriesYValue+YShift)*YScale%n%n");
      if (chart.getTitle() != null)
         formatter.format("\\begin{figure}%n");
      formatter.format("\\begin{center}%n");
      formatter.format("\\begin{tikzpicture}[y=%scm]%n", yunit);
      formatter.format("\\footnotesize%n");
      if(grid)
         formatter.format("\\draw[color=lightgray] (%s) grid[ystep=%s] (%s);%n",
            (Math.min(YAxis.getAxis().getRange().getLowerBound(), YAxis.getTwinAxisPosition())-YAxis.getTwinAxisPosition()) * YScale,
            ystepGrid*YScale,
            (Math.max(YAxis.getAxis().getRange().getUpperBound(), YAxis.getTwinAxisPosition())-YAxis.getTwinAxisPosition()) * YScale );
      formatter.format("%s", YAxis.toLatex(YScale) );

      formatter.format("%s", dataset.toLatex(YScale, YAxis.getTwinAxisPosition(),      YAxis.getAxis().getLowerBound(), YAxis.getAxis().getUpperBound()));

      formatter.format("\\end{tikzpicture}%n");
      formatter.format("\\end{center}%n");
      if (chart.getTitle() != null) {
         formatter.format("\\caption{");
         formatter.format(chart.getTitle().getText());
         formatter.format("}%n\\end{figure}%n");
      }
      formatter.format("\\end{document}%n");
      return formatter.toString();
*/
   }

}

/**
 * @}
 */