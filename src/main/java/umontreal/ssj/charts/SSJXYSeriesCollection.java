/*
 * Class:        SSJXYSeriesCollection
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

import   org.jfree.data.xy.XYDataset;
import   org.jfree.chart.renderer.xy.XYItemRenderer;

import   java.util.Locale;
import   java.util.Formatter;
import   java.awt.Color;

/**
 * Stores data used in a `XYChart`. This class provides tools to manage data
 * sets and rendering options, and modify plot color, plot style, and marks
 * on points for each series.
 *
 * <div class="SSJ-bigskip"></div>
 */
public abstract class SSJXYSeriesCollection {
   protected XYItemRenderer renderer;
   protected XYDataset seriesCollection;

   /**
    * @name Data control methods
    * @{
    */

   /**
    * Returns the @f$x@f$-value at the specified index in the specified
    * series.
    *  @param series       required series value.
    *  @param index        value’s index.
    *  @return @f$x@f$-value at the specified index in the specified
    * series.
    */
   public double getX (int series, int index) {
      return seriesCollection.getXValue(series, index);
   }

   /**
    * Returns the @f$y@f$-value at the specified index in the specified
    * series.
    *  @param series       required series value.
    *  @param index        value’s index.
    *  @return @f$y@f$-value at the specified index in the specified
    * series.
    */
   public double getY (int series, int index) {
      return seriesCollection.getYValue(series, index);
   }

   /**
    * Returns the `XYDataset` object associated with the current object.
    *  @return `XYDataset` object associated with the current variable.
    */
   public XYDataset getSeriesCollection() {
      return seriesCollection;
   }

   /**
    * Returns domain (@f$x@f$-coordinates) min and max values.
    *  @return domain min and max values.
    */
   public double[] getDomainBounds() {
      double max=-1.0e307, min=1.0e307;

      if(seriesCollection.getSeriesCount() != 0 && seriesCollection.getItemCount(0) != 0)
         max = min = seriesCollection.getXValue(0, 0);

      for(int i = 0; i < seriesCollection.getSeriesCount(); i++) {
         for( int j = 0; j < seriesCollection.getItemCount(i); j++) {
            max = Math.max(max, seriesCollection.getXValue(i, j));
            min = Math.min(min, seriesCollection.getXValue(i, j));
         }
      }

      double[] retour = {min, max};
      return retour;
   }

   /**
    * Returns range (@f$y@f$-coordinates) min and max values.
    *  @return range min and max values.
    */
   public double[] getRangeBounds() {
      double max=-1.7e307, min=1.7e307;

      if(seriesCollection.getSeriesCount() != 0 && seriesCollection.getItemCount(0) != 0)
         max = min = seriesCollection.getYValue(0, 0);

      for(int i = 0; i < seriesCollection.getSeriesCount(); i++) {
         for( int j = 0; j < seriesCollection.getItemCount(i); j++) {
            max = Math.max(max, seriesCollection.getYValue(i, j));
            min = Math.min(min, seriesCollection.getYValue(i, j));
         }
      }

      double[] retour = {min, max};
      return retour;
   }

   /**
    * Returns in a `String` all data contained in the current object.
    *  @return All data contained in the current object as a  String.
    */
   public String toString() {
      Formatter formatter = new Formatter(Locale.US);
      for(int i = 0; i < seriesCollection.getSeriesCount(); i++) {
         formatter.format(" Series " + i + " : %n");
         for(int j = 0; j < seriesCollection.getItemCount(i); j++)
            formatter.format("%15e,%15e%n", getX(i, j), getY(i, j));
      }
      return formatter.toString();
   }

   /**
    * @}
    */

   /**
    * @name Rendering methods
    * @{
    */

   /**
    * Returns the `XYItemRenderer` object associated with the current
    * object.
    *  @return `XYItemRenderer` object associated with the current
    * variable.
    */
   public XYItemRenderer getRenderer() {
      return renderer;
   }

   /**
    * Sets the `XYItemRenderer` object associated with the current
    * variable. This object determines the chart JFreeChart look, produced
    * by method `view` in class  @ref umontreal.ssj.charts.XYChart.
    *  @param renderer      new `XYItemRenderer` object.
    */
   public void setRenderer (XYItemRenderer renderer) {
      this.renderer = renderer;
   }

   /**
    * Gets the current plotting color of the selected series.
    *  @return current plotting color.
    */
   public Color getColor (int series) {
      return (Color)renderer.getSeriesPaint(series);
   }

   /**
    * Sets a new plotting color to the series @f$series@f$.
    *  @param series       series index.
    *  @param color        plotting color.
    */
   public void setColor (int series, Color color) {
      renderer.setSeriesPaint(series, color);
   }

   /**
    * Formats and returns a string containing a LaTeX-compatible source
    * code which represents this data series collection. The original
    * datasets are shifted and scaled with the `XShift`, `YShift`,
    * `XScale` and `YScale` parameters. `xmin`, `xmax`, `ymin` and `ymax`
    * represent the chart bounds.
    *  @param XScale       Domain original data scale.
    *  @param YScale       Range original data scale.
    *  @param XShift       Domain original data shift value.
    *  @param YShift       Range original data shift value.
    *  @param xmin         Domain min bound.
    *  @param xmax         Domain nax bound.
    *  @param ymin         Range min bound.
    *  @param ymax         Range nax bound.
    *  @return TikZ code.
    */
   public abstract String toLatex (double XScale, double YScale,
                                   double XShift, double YShift,
                                   double xmin, double xmax,
                                   double ymin, double ymax);


   /**
    * Converts a java Color object into a friendly and readable LaTeX/xcolor string.
    *
    * @param   color    in color.
    * @return           friendly color with string format as possible, null otherwise.
    */
   protected static String detectXColorClassic(Color color) {
      String retour = null;

      int red = color.getRed();
      int green = color.getGreen();
      int blue = color.getBlue();

      // On utilise pas la method Color.equals(Color ) car on ne veut pas tester le parametre de transparence : Alpha
      if (   red == Color.GREEN.getRed()
          && blue == Color.GREEN.getBlue()
          && green == Color.GREEN.getGreen())
         return "green";
      else if (   red == Color.RED.getRed()
               && blue == Color.RED.getBlue()
               && green == Color.RED.getGreen())
         return "red";
      else if (   red == Color.WHITE.getRed()
               && blue == Color.WHITE.getBlue()
               && green == Color.WHITE.getGreen())
         return "white";
      else if (   red == Color.GRAY.getRed()
               && blue == Color.GRAY.getBlue()
               && green == Color.GRAY.getGreen())
          return "gray";
      else if (   red == Color.BLACK.getRed()
               && blue == Color.BLACK.getBlue()
               && green == Color.BLACK.getGreen())
          return "black";
      else if (   red == Color.YELLOW.getRed()
               && blue == Color.YELLOW.getBlue()
               && green == Color.YELLOW.getGreen())
          return "yellow";
      else if (   red == Color.MAGENTA.getRed()
               && blue == Color.MAGENTA.getBlue()
               && green == Color.MAGENTA.getGreen())
          return "magenta";
      else if (   red == Color.CYAN.getRed()
               && blue == Color.CYAN.getBlue()
               && green == Color.CYAN.getGreen())
          return "cyan";
      else if (   red == Color.BLUE.getRed()
               && blue == Color.BLUE.getBlue()
               && green == Color.BLUE.getGreen())
          return "blue";
      else if (   red == Color.DARK_GRAY.getRed()
               && blue == Color.DARK_GRAY.getBlue()
               && green == Color.DARK_GRAY.getGreen())
          return "darkgray";
      else if (   red == Color.LIGHT_GRAY.getRed()
               && blue == Color.LIGHT_GRAY.getBlue()
               && green == Color.LIGHT_GRAY.getGreen())
          return "lightgray";
      else if (   red == Color.ORANGE.getRed()
               && blue == Color.ORANGE.getBlue()
               && green == Color.ORANGE.getGreen())
          return "orange";
      else if (   red == Color.PINK.getRed()
               && blue == Color.PINK.getBlue()
               && green == Color.PINK.getGreen())
          return "pink";


      if (red == 192 && blue == 128 && green == 64)
         return "brown";
     else if (red == 128 && blue == 128 && green == 0)
         return "olive";
      else if (red == 128 && blue == 0 && green == 128)
         return "violet";
      else if (red == 192 && blue == 0 && green ==64)
         return "purple";
      else return null;
   }

   /**
    * Gives the default color associated with a series
    *
    * @param   index Index of the series in the XYDataset object.
    * @return        default color object.
    */
   protected static Color getDefaultColor(int index) {
      if(index%6 == 0)
         return Color.RED;
      else if(index%6 == 1)
         return Color.BLUE;
      else if(index%6 == 2)
         return Color.GREEN;
      else if(index%6 == 3)
         return Color.YELLOW;
      else if(index%6 == 4)
         return Color.MAGENTA;
      else
         return Color.CYAN;
   }

   // Returns maximum value in table t
   protected static double max (double[] t) {
      double aux = t[0];
      for (int i=1 ; i < t.length ; i++)
         if (t[i] > aux)
            aux = t[i];
      return aux ;
   }
   // Returns minimum value in table t
   protected static double min (double[] t) {
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