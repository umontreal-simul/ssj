/*
 * Class:        CategoryChart
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

import org.jfree.chart.JFreeChart;
import javax.swing.JFrame;

/**
 * This class provides tools to create charts from data in a simple way. Its
 * main feature is to produce TikZ/PGF (see WWW link
 * [http://sourceforge.net/projects/pgf/](http://sourceforge.net/projects/pgf/))
 * compatible source code which can be included in LaTeX documents, but it
 * can also produce charts in other formats. One can easily create a new
 * chart, and customize its appearance using methods of this class, with the
 * encapsulated  @ref umontreal.ssj.charts.SSJCategorySeriesCollection object
 * representing the data, and an  @ref umontreal.ssj.charts.Axis object
 * representing the axis. All these classes depend on the `JFreeChart` API
 * (see WWW link
 * [http://www.jfree.org/jfreechart/](http://www.jfree.org/jfreechart/))
 * which provides tools to build charts with Java, to draw them, and export
 * them to files. However, only basic features are used here.
 *
 * Moreover, `CategoryChart` provides methods to plot data using a MATLAB
 * friendly syntax. None of these methods provides new features; they just
 * propose a different syntax to create charts. Therefore some features are
 * unavailable when using these methods only.
 *
 * <div class="SSJ-bigskip"></div>
 */
public abstract class CategoryChart {

   protected Axis YAxis;
   protected SSJCategorySeriesCollection dataset;
   protected JFreeChart chart;
   protected boolean latexDocFlag = true;

   protected boolean autoRange;
   protected double[] manualRange;

   protected boolean grid = false;
   protected double ystepGrid;

   final protected double BOR = 0.1;

   /**
    * Returns the `JFreeChart` object associated with this chart.
    *  @return the associated JFreeChart object.
    */
   public JFreeChart getJFreeChart() {
      return chart;
   }

   /**
    * Returns the chart’s range axis (@f$y@f$-axis) object.
    *  @return chart’s range axis (@f$y@f$-axis) object.
    */
   public Axis getYAxis() {
      return YAxis;
   }

  /**
    * Displays the chart on the screen using Swing. 
    *
    *  @param width     frame width
    *  @param height    frame height
    *
    *  @return frame containing the chart
    */
   public abstract JFrame view (int width, int height);

   /**
    * Gets the current chart title.
    *  @return Chart title.
    */
   public String getTitle() {
      return chart.getTitle().getText();
   }

   /**
    * Sets a title to this chart. This title will appear on the chart
    * displayed by method @ref #view.
    *
    *  @param title        chart title
    */
   public void setTitle (String title) {
      chart.setTitle(title);
   }

   /**
    * Sets chart @f$y@f$ range to automatic values.
    */
   public void setAutoRange () {
      autoRange = true;

      double BorneMin = Math.abs((dataset.getRangeBounds())[0]);
      double BorneMax = Math.abs((dataset.getRangeBounds())[1]);

      double max = Math.max(BorneMin,BorneMax) * BOR;
      YAxis.getAxis().setLowerBound(BorneMin - max);
      YAxis.getAxis().setUpperBound(BorneMax + max);
      YAxis.setLabelsAuto();
   }

   /**
    * Sets new @f$y@f$-axis bounds, using the format: `range` = [<tt>ymin,
    * ymax</tt>].
    *  @param range        new axis ranges.
    */
   private void setManualRange (double[] range)  {
      if(range.length != 2)
         throw new IllegalArgumentException (
             "range must have the format: [ymin, ymax]");
      autoRange = false;
      YAxis.getAxis().setLowerBound(Math.min(range[0],range[1]));
      YAxis.getAxis().setUpperBound(Math.max(range[0],range[1]));
   }

   /**
    * Puts a grid on the background. It is important to note that the grid
    * is always shifted in such a way that it contains the axes. Thus, the
    * grid does not always have an intersection at the corner points; this
    * occurs only if the corner points are multiples of the steps: `xstep`
    * and `ystep` sets the step in each direction.
    *  @param xstep        sets the step in the x-direction.
    *  @param ystep        sets the step in the y-direction.
    */
   public void enableGrid (double xstep, double ystep) {
      this.grid = true;
      this.ystepGrid = ystep;
   }

   /**
    * Disables the background grid.
    */
   public void disableGrid() {
      this.grid = false;
   }

   /**
    * @name Latex-specific methods
    * @{
    */

   /**
    * Transforms the chart into LaTeX form and returns it as a `String`.
    */
   public abstract String toLatex (double width, double height);

   /**
    * Same as in  @ref XYChart.
    */
   public void setLatexDocFlag (boolean flag) {
      latexDocFlag = flag;
   }


   protected double computeYScale (double position) {
      double[] bounds = new double[2];
      bounds[0] = YAxis.getAxis().getLowerBound();
      bounds[1] = YAxis.getAxis().getUpperBound();

      if (position < bounds[0])
         bounds[0] = position;
      if (position > bounds[1])
         bounds[1] = position;
      bounds[0] -= position;
      bounds[1] -= position;
      return computeScale (bounds);
   }

   protected double computeScale (double[] bounds) {
      int tenPowerRatio = 0;
      // echelle < 1 si les valeurs sont grandes
      while (bounds[1] > 1000 || bounds[0] < -1000) {
         bounds[1] /= 10;
         bounds[0] /= 10;
         tenPowerRatio++;
      }
      // echelle > 1 si les valeurs sont petites
      while (bounds[1] < 100 && bounds[0] > -100) {
         bounds[1] *= 10;
         bounds[0] *= 10;
         tenPowerRatio--;
      }
      return 1/Math.pow(10, tenPowerRatio);
   }
}

/**
 * @}
 */