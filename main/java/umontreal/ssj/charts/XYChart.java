/*
 * Class:        XYChart
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
import java.io.*;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.*;
import org.jfree.chart.plot.XYPlot;
import org.jfree.ui.TextAnchor;
import javax.swing.JFrame;

/**
 * This class provides tools to create charts from data in a simple way. Its
 * main feature is to produce TikZ/PGF (see WWW link
 * [http://sourceforge.net/projects/pgf/](http://sourceforge.net/projects/pgf/))
 * compatible source code which can be included in LaTeX documents, but it
 * can also produce charts in other formats. One can easily create a new
 * chart, and customize its appearance using methods of this class, with the
 * encapsulated  @ref umontreal.ssj.charts.SSJXYSeriesCollection object
 * representing the data, and the two  @ref umontreal.ssj.charts.Axis objects
 * representing the axes. All these classes depend on the `JFreeChart` API
 * (see WWW link
 * [http://www.jfree.org/jfreechart/](http://www.jfree.org/jfreechart/))
 * which provides tools to build charts with Java, to draw them, and export
 * them to files. However, only basic features are used here.
 *
 * Moreover, `XYChart` provides methods to plot data using a MATLAB friendly
 * syntax. None of these methods provides new features; they just propose a
 * different syntax to create charts. Therefore some features are unavailable
 * when using these methods only.
 *
 * <div class="SSJ-bigskip"></div>
 */
public abstract class XYChart {
   protected Axis XAxis;
   protected Axis YAxis;

   protected SSJXYSeriesCollection dataset;
   protected JFreeChart chart;
   protected boolean latexDocFlag = true;

   protected boolean autoRange;
   protected double[] manualRange;

   protected boolean grid = false;
   protected double xstepGrid;
   protected double ystepGrid;

   // this flag is set true when plotting probabilities. In that case,
   // y is always >= 0.
   protected boolean probFlag = false;

   protected double chartMargin = 0.02;   // margin around the chart

   /**
    * Returns the `JFreeChart` object associated with this chart.
    *  @return the associated JFreeChart object.
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
    * Displays chart on the screen using Swing.
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
    * displayed by method  #view.
    *  @param title        chart title.
    */
   public void setTitle (String title) {
      chart.setTitle(title);
   }

   /**
    * Must be set `true` when plotting probabilities, `false` otherwise.
    *  @param flag         `true` for plotting probabilities
    */
   public void setprobFlag (boolean flag) {
      probFlag = flag;
   }

   /**
    * The @f$x@f$ and the @f$y@f$ ranges of the chart are set
    * automatically.
    */
   public void setAutoRange() {
      setAutoRange (false, false, true, true);
   }

   /**
    * The @f$x@f$ and the @f$y@f$ ranges of the chart are set
    * automatically. If `right` is `true`, the vertical axis will be on
    * the left of the points, otherwise on the right. If `top` is `true`,
    * the horizontal axis will be under the points, otherwise above the
    * points.
    *  @param right        true if the x-values on the right of axis.
    *  @param top          true if the y-values on the top of axis.
    */
   public void setAutoRange (boolean right, boolean top) {
      setAutoRange (false, false, right, top);
   }

   private double[] adjustRangeBounds (double bmin, double bmax) {
      // resets chart lower and upper bounds to round values.
      // Returns corrected [lowerBound, upperBound]

      double del = (bmax - bmin)/20.0;    // Choose 20 intervals to round
      int a = (int) Math.floor(0.5 + Math.log10(del));
      double d = Math.pow(10.0, (double) a);     // power of 10
      double lower = d*Math.ceil((bmin - del)/d);
      if (lower > bmin)
         lower -= d;
      if (0 == Math.abs(bmin))
         lower = 0;
      double upper = d*Math.floor((bmax + del)/d);
      if (upper < bmax)
         upper += d;
      double [] range = new double[2];
      range[0] = lower;
      range[1] = upper;
      return range;
   }

   protected void setAutoRange (boolean xZero, boolean yZero, boolean right, boolean top) {
      // see description of setAxesZero
      autoRange = true;
      double BorneMin = (dataset.getDomainBounds())[0];
      double BorneMax = (dataset.getDomainBounds())[1];
      double del;
      if (BorneMax - BorneMin < 1)
         del = (BorneMax - BorneMin) * chartMargin;
      else
         del = chartMargin;
      if (BorneMin < 0.0) BorneMin *= 1.0 + del;
      else BorneMin *= 1.0 - del;
      if (BorneMax < 0.0) BorneMax *= 1.0 - del;
      else BorneMax *= 1.0 + del;
      double [] newRange = new double[2];
      newRange = adjustRangeBounds (BorneMin, BorneMax);
      if (probFlag && (BorneMin == 0.0))
         newRange[0] = 0.0;
      XAxis.getAxis().setLowerBound(newRange[0]);
      XAxis.getAxis().setUpperBound(newRange[1]);

      BorneMin = (dataset.getRangeBounds())[0];
      BorneMax = (dataset.getRangeBounds())[1];
      if (BorneMax - BorneMin < 1)
         del = (BorneMax - BorneMin) * chartMargin;
      else
         del = chartMargin;
      if (BorneMin < 0.0) BorneMin *= 1.0 + del;
      else BorneMin *= 1.0 - del;
      if (BorneMax < 0.0) BorneMax *= 1.0 - del;
      else BorneMax *= 1.0 + del;
      newRange = adjustRangeBounds (BorneMin, BorneMax);
      if (probFlag && (newRange[0] <= 0.0))   // probabilities are always >= 0
         newRange[0] = 0.0;
      YAxis.getAxis().setLowerBound(newRange[0]);
      YAxis.getAxis().setUpperBound(newRange[1]);

      if (xZero)
         XAxis.setTwinAxisPosition(0);
      else {
         if (right)
            XAxis.setTwinAxisPosition(XAxis.getAxis().getLowerBound());
         else
            XAxis.setTwinAxisPosition(XAxis.getAxis().getUpperBound());
      }

      if (yZero)
         YAxis.setTwinAxisPosition(0);
      else {
         if (top)
            YAxis.setTwinAxisPosition(YAxis.getAxis().getLowerBound());
         else
            YAxis.setTwinAxisPosition(YAxis.getAxis().getUpperBound());
      }
   }

/**
 * The @f$x@f$ and the @f$y@f$ ranges of the chart are set automatically. If
 * `xZero` is `true`, the vertical axis will pass through the point @f$(0,
 * y)@f$. If `yZero` is `true`, the horizontal axis will pass through the
 * point @f$(x, 0)@f$.
 *  @param xZero        true if vertical axis passes through point 0
 *  @param yZero        true if horizontal axis passes through point 0
 */
public void setAutoRange00 (boolean xZero, boolean yZero) {
      setAutoRange (xZero, yZero, true, true);
   }

   /**
    * Sets the @f$x@f$ and @f$y@f$ ranges of the chart using the format:
    * `range = [xmin, xmax, ymin, ymax]`.
    *  @param range        new axis ranges.
    */
   public void setManualRange (double[] range) {
      setManualRange (range, false, false, true, true);
   }

   /**
    * Sets the @f$x@f$ and @f$y@f$ ranges of the chart using the format:
    * `range = [xmin, xmax, ymin, ymax]`. If `right` is `true`, the
    * vertical axis will be on the left of the points, otherwise on the
    * right. If `top` is `true`, the horizontal axis will be under the
    * points, otherwise above the points.
    *  @param range        new axis ranges.
    *  @param right        true if the x-values on the right.
    *  @param top          true if the y-values on the top.
    */
   public void setManualRange (double[] range, boolean right, boolean top) {
      setManualRange (range, false, false, right, top);
   }


   private void setManualRange (double[] range, boolean xZero, boolean yZero,
                                boolean right, boolean top) {
      if (range.length != 4)
         throw new IllegalArgumentException (
             "range must have the format: [xmin, xmax, ymin, ymax]");
      autoRange = false;
      XAxis.getAxis().setLowerBound(Math.min(range[0],range[1]));
      XAxis.getAxis().setUpperBound(Math.max(range[0],range[1]));
      YAxis.getAxis().setLowerBound(Math.min(range[2],range[3]));
      YAxis.getAxis().setUpperBound(Math.max(range[2],range[3]));

      if (xZero)
         XAxis.setTwinAxisPosition(0);
      else {
         if (right)
            XAxis.setTwinAxisPosition(XAxis.getAxis().getLowerBound());
         else
            XAxis.setTwinAxisPosition(XAxis.getAxis().getUpperBound());
      }

      if (yZero)
         YAxis.setTwinAxisPosition(0);
      else {
         if (top)
            YAxis.setTwinAxisPosition(YAxis.getAxis().getLowerBound());
         else
            YAxis.setTwinAxisPosition(YAxis.getAxis().getUpperBound());
      }
   }

/**
 * Sets the @f$x@f$ and @f$y@f$ ranges of the chart using the format: `range
 * = [xmin, xmax, ymin, ymax]`. If `xZero` is `true`, the vertical axis will
 * pass through the point @f$(0, y)@f$. If `yZero` is `true`, the horizontal
 * axis will pass through the point @f$(x, 0)@f$.
  * @param range        new axis ranges.
 *  @param xZero        true if vertical axis passes through point 0
 *  @param yZero        true if horizontal axis passes through point 0
 */
public void setManualRange00 (double[] range, boolean xZero, boolean yZero) {
      setManualRange (range, xZero, yZero, true, true);
   }

   /**
    * Returns the chart margin, which is the fraction by which the chart
    * is enlarged on its borders. The default value is @f$0.02@f$.
    */
   public double getChartMargin() {
      return chartMargin;
   }

   /**
    * Sets the chart margin to `margin`. It is the fraction by which the
    * chart is enlarged on its borders. Restriction: @f$\mathtt{margin}
    * \ge0@f$.
    *  @param margin       margin percentage amount.
    */
   public void setChartMargin (double margin) {
      if (margin < 0.0)
         throw new IllegalArgumentException ("margin < 0");
      chartMargin = margin;
   }

   /**
    * Synchronizes @f$x@f$-axis ticks to the @f$s@f$-th series
    * @f$x@f$-values.
    *  @param s            series.
    */
   public abstract void setTicksSynchro (int s);

   /**
    * Draws a vertical line on the chart at @f$x@f$-coordinate `x`. `name`
    * is written near the line at @f$y@f$ position `yfrac` (a fraction of
    * the @f$y@f$-size of the chart, 0 is the bottom, 1 is the top); if
    * `right` is `true`, `name` is written on the right of the line, else
    * on the left.
    *  @param x            @f$x@f$-coordinate of the line
    *  @param name         description of the line
    *  @param yfrac        @f$y@f$-position of name
    *  @param right        @f$x@f$-position of name
    */
   public void drawVerticalLine (double x, String name, double yfrac,
                                 boolean right) {
      double ybottom = YAxis.getAxis().getLowerBound();
      final Object o = this;
      if (this instanceof HistogramChart)
         ybottom = 0;
      double ytop = YAxis.getAxis().getUpperBound();
      XYLineAnnotation line = new XYLineAnnotation(x, ybottom, x, ytop);
      XYTextAnnotation text = new XYTextAnnotation(name, x, ytop*yfrac);
      if (!right)
         text.setTextAnchor(TextAnchor.HALF_ASCENT_RIGHT);
      else
         text.setTextAnchor(TextAnchor.HALF_ASCENT_LEFT);
      XYPlot plot = getJFreeChart().getXYPlot();
      plot.addAnnotation(line);
      plot.addAnnotation(text);
   }

   /**
    * @name Latex-specific methods
    * @{
    */

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
      this.xstepGrid = xstep;
      this.ystepGrid = ystep;
   }

   /**
    * Disables the background grid.
    */
   public void disableGrid() {
      this.grid = false;
   }

   /**
    * Exports the chart to a LaTeX source code using PGF/TikZ. This method
    * constructs and returns a string that can be written to a LaTeX
    * document to render the plot. `width` and `height` represents the
    * width and the height of the produced chart. These dimensions do not
    * take into account the axes and labels extra space. The `width` and
    * the `height` of the chart are measured in centimeters.
    *  @param width        Chart’s width in centimeters.
    *  @param height       Chart’s height in centimeters.
    *  @return LaTeX source code.
    */
   public abstract String toLatex (double width, double height);

   /**
    * Transforms the chart to LaTeX form and writes it in file `fileName`.
    * The chart’s width and height (in centimeters) are `width` and
    * `height`.
    */
   public void toLatexFile (String fileName, double width, double height) {
      String output = toLatex(width, height);
      Writer file = null;
      try {
         file = new FileWriter(fileName);
         file.write(output);
         file.close();
     } catch (IOException e) {
         System.err.println ("   toLatexFile:  cannot write to  " + fileName);
         e.printStackTrace();
         try {
            if (file != null)
               file.close();
         } catch (IOException ioe) {}
      }
  }

   /**
    * Flag to remove the <code>\\documentclass</code> (and other) commands
    * in the created LaTeX files. If `flag` is `true`, then when charts
    * are translated into LaTeX form, it will be as a self-contained file
    * that can be directly compiled with LaTeX. However, in this form, the
    * file cannot be included in another LaTeX file without causing
    * compilation errors because of the multiple instructions
    * <code>\\documentclass</code> and <code>\\begin{document}</code>. By
    * setting `flag` to `false`, these instructions will be removed from
    * the LaTeX chart files, which can then be included in a master LaTeX
    * file. By default, the flag is `true`.
    */
   public void setLatexDocFlag (boolean flag) {
      latexDocFlag = flag;
   }


   protected void setTick0Flags() {
      // Set flag true if first or last label is on perpendicular axis.
      // The label will be moved a little to the right (x-label), or above
      // (y-label) to prevent it from being on the perpendicular axis.
      // But it is unnecessary when graph begins or ends where label is;
      // in this case, flag is false.
      // We cannot put this method in Axis because it depends on the
      // other axis.
     double minAxis = Math.min (XAxis.getAxis().getRange().getLowerBound(),
                                                   XAxis.getTwinAxisPosition());
     double maxAxis = Math.max (XAxis.getAxis().getRange().getUpperBound(),
                                                     XAxis.getTwinAxisPosition());
     if (XAxis.getTwinAxisPosition() == minAxis ||
         XAxis.getTwinAxisPosition() == maxAxis)
        YAxis.setTick0Flag(false);
     else
        YAxis.setTick0Flag(true);

     minAxis = Math.min (YAxis.getAxis().getRange().getLowerBound(),
                                                   YAxis.getTwinAxisPosition());
     maxAxis = Math.max (YAxis.getAxis().getRange().getUpperBound(),
                                                     YAxis.getTwinAxisPosition());
     if (YAxis.getTwinAxisPosition() == minAxis ||
         YAxis.getTwinAxisPosition() == maxAxis)
        XAxis.setTick0Flag(false);
     else
        XAxis.setTick0Flag(true);
   }


   protected double computeXScale (double position) {
      double[] bounds = new double[2];
      bounds[0] = XAxis.getAxis().getLowerBound();
      bounds[1] = XAxis.getAxis().getUpperBound();

      if (position < bounds[0])
         bounds[0] = position;
      if (position > bounds[1])
         bounds[1] = position;
      bounds[0] -= position;
      bounds[1] -= position;
      return computeScale (bounds);
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