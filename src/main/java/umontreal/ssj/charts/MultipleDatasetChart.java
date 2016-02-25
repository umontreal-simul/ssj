/*
 * Class:        MultipleDatasetChart
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

import   org.jfree.chart.JFreeChart;
import   org.jfree.chart.ChartPanel;
import   org.jfree.chart.ChartFactory;
import   org.jfree.chart.axis.NumberAxis;
import   org.jfree.chart.plot.XYPlot;
import   org.jfree.chart.plot.PlotOrientation;

import   java.util.Locale;
import   java.util.Formatter;
import   java.util.ArrayList;
import   javax.swing.JFrame;

/**
 * Provides tools to plot many datasets on the same chart. This class is
 * mainly used to draw plots with different styles. Class
 * @ref umontreal.ssj.charts.XYChart and its subclasses are to be preferred
 * to draw simple charts with one style. Datasets are stored in an
 * `ArrayList`. The first dataset is called as the *primary dataset*.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class MultipleDatasetChart {

   protected ArrayList<SSJXYSeriesCollection> datasetList;
   protected Axis XAxis;
   protected Axis YAxis;
   protected JFreeChart chart;
   protected boolean latexDocFlag = true;

   protected boolean autoRange = true;
   protected double[] manualRange;

   protected boolean grid = false;
   protected double xstepGrid;
   protected double ystepGrid;

   /**
    * Initializes a new `MultipleDatasetChart`.
    */
   public MultipleDatasetChart() {
      super();

      // create the chart...
      chart = ChartFactory.createXYLineChart(
         null,                     // chart title
         null,                     // x axis label
         null,                     // y axis label
         null, // data
         PlotOrientation.VERTICAL,
         true,                     // include legend
         true,                     // tool tips
         false                     // urls
      );

      datasetList = new ArrayList<SSJXYSeriesCollection>();
      // Initialize axis variables
      XAxis = new Axis((NumberAxis)((XYPlot) chart.getPlot()).getDomainAxis(),
                 Axis.ORIENTATION_HORIZONTAL);
      YAxis = new Axis((NumberAxis)((XYPlot) chart.getPlot()).getRangeAxis(),
                 Axis.ORIENTATION_VERTICAL);
   }

   /**
    * Initializes a new `MultipleDatasetChart` instance. `title` sets a
    * title, `XLabel` is a short description of the @f$x@f$-axis, and
    * `YLabel` is a short description of the @f$y@f$-axis.
    *  @param title        chart title.
    *  @param XLabel       Label on @f$x@f$-axis.
    *  @param YLabel       Label on @f$y@f$-axis.
    */
   public MultipleDatasetChart (String title, String XLabel, String YLabel) {
      // create the chart...
      chart = ChartFactory.createXYLineChart(
         title,                    // chart title
         XLabel,                   // x axis label
         YLabel,                   // y axis label
         null,                     // data
         PlotOrientation.VERTICAL,
         true,                     // include legend
         true,                     // tool tips
         false                     // urls
      );

      datasetList = new ArrayList<SSJXYSeriesCollection>();
      //Initialize axis variables
      XAxis = new Axis((NumberAxis)((XYPlot) chart.getPlot()).getDomainAxis(),
              Axis.ORIENTATION_HORIZONTAL);
      YAxis = new Axis((NumberAxis)((XYPlot) chart.getPlot()).getRangeAxis(),
              Axis.ORIENTATION_VERTICAL);
   }

   /**
    * Returns the `JFreeChart` variable associated with this chart.
    *  @return the associated JFreeChart variable.
    */
   public JFreeChart getJFreeChart() {
      return chart;
   }

   /**
    * Returns the chart’s domain axis (@f$x@f$-axis) object.
    *  @return chart’s domain axis (@f$x@f$-axis) object.
    */
   public Axis getXAxis() {
      return XAxis;
   }

   /**
    * Returns the chart’s range axis (@f$y@f$-axis) object.
    *  @return chart’s range axis (@f$y@f$-axis) object.
    */
   public Axis getYAxis() {
      return YAxis;
   }

   /**
    * Gets the current chart title.
    *  @return Chart title.
    */
   public String getTitle() {
      return chart.getTitle().getText();
   }

   /**
    * Sets a title to the chart. This title will appear on the chart
    * displayed by method  #view.
    *  @param title        chart title.
    */
   public void setTitle (String title) {
      chart.setTitle(title);
   }

   /**
    * Sets chart range to automatic values.
    */
   public void setAutoRange() {
      autoRange = true;
      double[][] temp = new double[2][datasetList.size()];
      for(int i = 0; i<datasetList.size(); i++) {
         temp[0][i] = (datasetList.get(i).getDomainBounds())[0];
         temp[1][i] = (datasetList.get(i).getDomainBounds())[1];
      }
      XAxis.getAxis().setLowerBound(min(temp[0]));
      XAxis.getAxis().setUpperBound(max(temp[1]));

      for(int i = 0; i<datasetList.size(); i++) {
         temp[0][i] = (datasetList.get(i).getRangeBounds())[0];
         temp[1][i] = (datasetList.get(i).getRangeBounds())[1];
      }
      YAxis.getAxis().setLowerBound(min(temp[0]));
      YAxis.getAxis().setUpperBound(max(temp[1]));
   }

   /**
    * Sets new @f$x@f$-axis and @f$y@f$-axis bounds, with format:
    * `axisRange` = [xmin, xmax, ymin, ymax].
    *  @param axisRange    new axis ranges.
    */
   public void setManualRange (double[] axisRange) {
      if(axisRange.length != 4)
         throw new IllegalArgumentException("axisRange must share the format: [xmin, xmax, ymin, ymax]");
      autoRange = false;
      XAxis.getAxis().setLowerBound(axisRange[0]);
      XAxis.getAxis().setUpperBound(axisRange[1]);
      YAxis.getAxis().setLowerBound(axisRange[2]);
      YAxis.getAxis().setUpperBound(axisRange[3]);
   }

   /**
    * Adds a new dataset to the chart at the end of the list and returns
    * its position.
    *  @param dataset      dataset to add.
    *  @return the dataset position in the list.
    */
   public int add (SSJXYSeriesCollection dataset) {
      ((XYPlot)chart.getPlot()).setDataset(datasetList.size(),
                  dataset.getSeriesCollection());
      ((XYPlot)chart.getPlot()).setRenderer(datasetList.size(),
                  dataset.getRenderer());
      datasetList.add(dataset);
      if(datasetList.size() == 1) {
         XAxis.setLabelsAuto();
         YAxis.setLabelsAuto();
      }
      return datasetList.size()-1;
   }

   /**
    * Gets the primary dataset.
    *  @return dataset.
    */
   public SSJXYSeriesCollection get() {
      return datasetList.get(0);
   }

   /**
    * Sets the primary dataset for the plot, replacing the existing
    * dataset if there is one.
    *  @param dataset      the new primary dataset.
    */
   public void set (SSJXYSeriesCollection dataset) {
      ((XYPlot)chart.getPlot()).setDataset(dataset.getSeriesCollection());
      ((XYPlot)chart.getPlot()).setRenderer(dataset.getRenderer());
      datasetList.set(0, dataset);
   }

   /**
    * Gets the element at the specified position in the dataset list.
    *  @param datasetNum   position in the dataset list.
    *  @return dataset.
    */
   public SSJXYSeriesCollection get (int datasetNum) {
      return datasetList.get(datasetNum);
   }

   /**
    * Replaces the element at the specified position in the dataset list
    * with the specified element.
    *  @param datasetNum   position in the dataset list.
    *  @param dataset      dataset list.
    */
   public void set (int datasetNum, SSJXYSeriesCollection dataset) {
      ((XYPlot)chart.getPlot()).setDataset(datasetNum, dataset.getSeriesCollection());
      ((XYPlot)chart.getPlot()).setRenderer(datasetNum, dataset.getRenderer());
      datasetList.add(datasetNum, dataset);
   }

   /**
    * Returns the dataset list.
    *  @return dataset list.
    */
   public ArrayList<SSJXYSeriesCollection> getList() {
      return datasetList;
   }

   /**
    * Displays chart on the screen using Swing. This method creates an
    * application containing a chart panel displaying the chart. The
    * created frame is positioned on-screen, and displayed before it is
    * returned. The `width` and the `height` of the chart are measured in
    * pixels.
    *  @param width        frame width in pixels.
    *  @param height       frame height in pixels.
    */
   public JFrame view (int width, int height) {
      JFrame myFrame;
      if (chart.getTitle() != null)
         myFrame = new JFrame("MultipleDatasetChart from SSJ: " +
                               chart.getTitle().getText());
      else
         myFrame = new JFrame("MultipleDatasetChart from SSJ");
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
    * Puts grid on the background. It is important to note that the grid
    * is always placed in such a way that it contains the axes. Thus, the
    * grid does not always have an intersection at the corner points; this
    * occurs only if the corner points are multiples of the stepping.
    * `xstep` and `ystep` sets the stepping in each direction.
    *  @param xstep        sets the stepping in the @f$x@f$-direction.
    *  @param ystep        sets the stepping in the @f$y@f$-direction.
    */
   public void enableGrid (double xstep, double ystep) {
      this.grid = true;
      this.xstepGrid = xstep;
      this.ystepGrid = ystep;
   }

   /**
    * Disables the background grid.
    */
   public void disableGrid () {
      this.grid = false;
   }

   /**
    * @name LaTeX-specific method
    * @{
    */

   /**
    * Same as in  @ref XYChart.
    *  @param width        Chart’s width in centimeters.
    *  @param height       Chart’s height in centimeters.
    */
   public String toLatex (double width, double height) {
      double xunit, yunit;
      double[] save = new double[4];

      if(datasetList.size() == 0)
         throw new IllegalArgumentException("Empty chart");

      //Calcul des parametres d'echelle et de decalage
      double XScale = computeXScale(XAxis.getTwinAxisPosition());
      double YScale = computeYScale(YAxis.getTwinAxisPosition());

      xunit = width / ( (Math.max(XAxis.getAxis().getRange().getUpperBound(), XAxis.getTwinAxisPosition()) * XScale) - (Math.min(XAxis.getAxis().getRange().getLowerBound(), XAxis.getTwinAxisPosition()) * XScale) );
      //taille d'une unite en x et en cm dans l'objet "tikzpicture"
      yunit = height / ( (Math.max(YAxis.getAxis().getRange().getUpperBound(), YAxis.getTwinAxisPosition()) * YScale) - (Math.min(YAxis.getAxis().getRange().getLowerBound(), YAxis.getTwinAxisPosition()) * YScale) );
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
      if(grid)
         formatter.format("\\draw[color=lightgray] (%s, %s) grid[xstep = %s, ystep=%s] (%s, %s);%n",
            (Math.min(XAxis.getAxis().getRange().getLowerBound(), XAxis.getTwinAxisPosition())-XAxis.getTwinAxisPosition()) * XScale,
            (Math.min(YAxis.getAxis().getRange().getLowerBound(), YAxis.getTwinAxisPosition())-YAxis.getTwinAxisPosition()) * YScale,
            xstepGrid*XScale, ystepGrid*YScale,
            (Math.max(XAxis.getAxis().getRange().getUpperBound(), XAxis.getTwinAxisPosition())-XAxis.getTwinAxisPosition()) * XScale,
            (Math.max(YAxis.getAxis().getRange().getUpperBound(), YAxis.getTwinAxisPosition())-YAxis.getTwinAxisPosition()) * YScale );

      formatter.format("%s", XAxis.toLatex(XScale) );
      formatter.format("%s", YAxis.toLatex(YScale) );

      for(int i = 0; i < datasetList.size(); i++)
         formatter.format("%s", datasetList.get(i).toLatex(XScale, YScale,
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

   /**
    * Same as in  @ref XYChart.
    */
   public void setLatexDocFlag (boolean flag) {
      latexDocFlag = flag;
   }


   protected double computeXScale (double position) {
      double[] bounds = new double[2];
      bounds[0] = XAxis.getAxis().getLowerBound();
      bounds[1] = XAxis.getAxis().getUpperBound();

      if(position<bounds[0])
         bounds[0] = position;
      if(position>bounds[1])
         bounds[1] = position;
      bounds[0] -= position;
      bounds[1] -= position;
      return computeScale(bounds);
   }

   protected double computeYScale (double position) {
      double[] bounds = new double[2];
      bounds[0] = YAxis.getAxis().getLowerBound();
      bounds[1] = YAxis.getAxis().getUpperBound();

      if(position<bounds[0])
         bounds[0] = position;
      if(position>bounds[1])
         bounds[1] = position;
      bounds[0] -= position;
      bounds[1] -= position;
      return computeScale(bounds);
   }


   protected double computeScale (double[] bounds) {
      int tenPowerRatio=0;
      //echelle < 1 si les valeurs sont grandes
      while(bounds[1] > 1000 || bounds[0] < -1000) {
         bounds[1] /= 10;
         bounds[0] /= 10;
         tenPowerRatio++;
      }
      //echelle > 1 si les valeurs sont petites
      while(bounds[1]<100 && bounds[0] > -100) {
         bounds[1] *= 10;
         bounds[0] *= 10;
         tenPowerRatio--;
      }
      return 1/Math.pow(10, tenPowerRatio);
   }

   private static double max (double[] t) {
      double aux = t[0];
      for (int i=1 ; i < t.length ; i++)
         if (t[i] > aux)
            aux = t[i];
      return aux ;
   }

   private static double min (double[] t) {
      double aux = t[0];
      for (int i=1 ; i < t.length ; i++)
         if (t[i] < aux)
            aux = t[i];
      return aux ;
   }
}

/**
 * @}
 */