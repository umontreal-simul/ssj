/*
 * Class:        ContinuousDistChart
 * Description:  
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        May 2008
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
   import umontreal.ssj.probdist.ContinuousDistribution;

import org.jfree.chart.JFreeChart;
import javax.swing.JFrame;

/**
 * This class provides tools to plot the density and the cumulative
 * probability of a continuous probability distribution.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class ContinuousDistChart {
   protected ContinuousDistribution dist;
   protected double a,b;
   protected int m;
   protected XYLineChart cdfChart;
   protected XYLineChart densityChart;

   private void init() {
      double[][] cdf = new double[2][m+1];
      double[][] density = new double[2][m+1];
      double h = (b - a) / m;
      double x;
      int coex = 0;

      try {
         for (int i = 0; i <= m; i++) {
            x = a + i*h;
            cdf[0][i] = x;
            cdf[1][i] = dist.cdf (x);
         }
         cdfChart = new XYLineChart("cdf: " + dist.toString(), "", "", cdf);
      } catch (UnsupportedOperationException e) {
         coex++;
         System.err.println (e);
//         e.printStackTrace();
      }

      try {
         for (int i = 0; i <= m; i++) {
            x = a + i*h;
            density[0][i] = x;
            density[1][i] = dist.density (x);
         }
         densityChart = new XYLineChart("density: " + dist.toString(),
                                        "", "", density);
      } catch (UnsupportedOperationException e) {
         System.err.println (e);
         if (coex == 1)
            throw e;
      }
      cdfChart.setprobFlag (true);
      densityChart.setprobFlag (true);
   }

   /**
    * Constructor for a new `ContinuousDistChart` instance. It will plot
    * the continuous distribution `dist` over the interval @f$[a,b]@f$,
    * using @f$m+1@f$ equidistant sample points.
    *  @param dist         continuous distribution to plot
    *  @param a            lower bound of interval
    *  @param b            upper bound of interval
    *  @param m            number of steps
    */
   public ContinuousDistChart (ContinuousDistribution dist, double a,
                               double b, int m) {
      this.dist = dist;
      this.a = a;
      this.b = b;
      this.m = m;
      init();
   }

   /**
    * Displays a chart of the cumulative distribution function (cdf) on
    * the screen using Swing. This method creates an application
    * containing a chart panel displaying the chart. The created frame is
    * positioned on-screen, and displayed before it is returned. The
    * `width` and the `height` of the chart are measured in pixels.
    *  @param width        frame width in pixels
    *  @param height       frame height in pixels
    *  @return frame containing the chart
    */
   public JFrame viewCdf (int width, int height) {
      return cdfChart.view(width, height);
   }

   /**
    * Similar to  #viewCdf, but for the probability density instead of the
    * cdf.
    *  @param width        frame width in pixels
    *  @param height       frame height in pixels
    *  @return frame containing the chart
    */
   public JFrame viewDensity (int width, int height) {
      return densityChart.view(width, height);
   }

   /**
    * Exports a chart of the cdf to a LaTeX source code using PGF/TikZ.
    * This method constructs and returns a string that can be written to a
    * LaTeX document to render the plot. `width` and `height` represents
    * the width and the height of the produced chart. These dimensions do
    * not take into account the axes and labels extra space. The `width`
    * and the `height` of the chart are measured in centimeters.
    *  @param width        Chart’s width in centimeters
    *  @param height       Chart’s height in centimeters
    *  @return LaTeX source code
    */
   public String toLatexCdf (int width, int height) {
      return cdfChart.toLatex(width, height);
   }

   /**
    * Similar to  #toLatexCdf, but for the probability density instead of
    * the cdf.
    *  @param width        Chart’s width in centimeters
    *  @param height       Chart’s height in centimeters
    *  @return LaTeX source code
    */
   public String toLatexDensity (int width, int height) {
      return densityChart.toLatex(width, height);
   }

   /**
    * Sets the parameters @f$a@f$, @f$b@f$ and @f$m@f$ for this object.
    *  @param a            lower bound of interval
    *  @param b            upper bound of interval
    *  @param m            number of points in the plot minus one
    */
   public void setParam (double a, double b, int m) {
      this.a = a;
      this.b = b;
      this.m = m;
      init();
   }

}