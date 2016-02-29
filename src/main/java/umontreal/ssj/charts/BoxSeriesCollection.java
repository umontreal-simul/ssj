/*
 * Class:        BoxSeriesCollection
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

import   org.jfree.chart.renderer.category.BoxAndWhiskerRenderer;
import   org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;

import   java.util.List;
import   java.util.ArrayList;
import   java.util.Locale;
import   java.util.Formatter;

/**
 * This class stores data used in a  @ref umontreal.ssj.charts.CategoryChart.
 * It also provides complementary tools to draw box-and-whisker plots; for
 * example, one may add or remove plots series and modify plot style. This
 * class is linked with the JFreeChart `DefaultBoxAndWhiskerCategoryDataset`
 * class to store data plots, and linked with the JFreeChart
 * `BoxAndWhiskerRenderer` to render the plots.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class BoxSeriesCollection extends SSJCategorySeriesCollection {
   final double BARWIDTH = 0.1;

   /**
    * Creates a new `BoxSeriesCollection` instance with an empty dataset.
    */
   public BoxSeriesCollection () {
      renderer = new BoxAndWhiskerRenderer();
      seriesCollection = new DefaultBoxAndWhiskerCategoryDataset ();
      ((BoxAndWhiskerRenderer)renderer).setMaximumBarWidth(BARWIDTH);
   }

   /**
    * Creates a new `BoxSeriesCollection` instance with default parameters
    * and input series `data`. Only *the first* `numPoints` of `data` will
    * taken into account.
    *  @param data         point sets.
    *  @param numPoints    Number of points
    */
   public BoxSeriesCollection (double[] data, int numPoints) {
      renderer = new BoxAndWhiskerRenderer();

      ((BoxAndWhiskerRenderer)renderer).setMaximumBarWidth(BARWIDTH);
      seriesCollection = new DefaultBoxAndWhiskerCategoryDataset ();

      DefaultBoxAndWhiskerCategoryDataset tempSeriesCollection =
            (DefaultBoxAndWhiskerCategoryDataset)seriesCollection;

      final List<Double> list = new ArrayList<Double>();
      for (int i = 0; i < numPoints; i ++)
         list.add(data[i]);

      tempSeriesCollection.add(list, 0, 0);
   }

   /**
    * Creates a new `BoxSeriesCollection` instance with default parameters
    * and given data series. The input parameter represents series of
    * point sets.
    *  @param data         series of point sets.
    */
   public BoxSeriesCollection (double[]... data) {
      renderer = new BoxAndWhiskerRenderer();
      seriesCollection = new DefaultBoxAndWhiskerCategoryDataset ();

      DefaultBoxAndWhiskerCategoryDataset tempSeriesCollection =
            (DefaultBoxAndWhiskerCategoryDataset)seriesCollection;

      for (int i = 0; i < data.length; i ++) {
         if (data[i].length == 0)
            throw new IllegalArgumentException("Unable to render the plot. data["
                                                   + i +"] contains no row");
         final List<Double> list = new ArrayList<Double>();
         for (int j = 0; j < data[i].length-1; j ++)
            list.add(data[i][j]);
         tempSeriesCollection.add(list, 0, "Serie " + i);
         list.clear();
      }
      ((BoxAndWhiskerRenderer)renderer).setMaximumBarWidth(BARWIDTH);
   }

   /**
    * Creates a new `BoxSeriesCollection` instance with default parameters
    * and given data series. The input parameter represents a
    * `DefaultBoxAndWhiskerCategoryDataset`.
    *  @param data         series of point sets.
    */
   public BoxSeriesCollection (DefaultBoxAndWhiskerCategoryDataset data) {
      renderer = new BoxAndWhiskerRenderer();
      ((BoxAndWhiskerRenderer)renderer).setFillBox(false);
      seriesCollection = data;
      ((BoxAndWhiskerRenderer)renderer).setMaximumBarWidth(BARWIDTH);
   }

   /**
    * @name Data control methods
    * @{
    */

   /**
    * Adds a data series into the series collection. Vector `data`
    * represents a point set.
    *  @param data         point sets.
    *  @return Integer that represent the new point set’s position in the
    * JFreeChart `DefaultBoxAndWhiskerXYDataset` object.
    */
   public int add (double[] data) {
      return add(data, data.length);
   }

   /**
    * Adds a data series into the series collection. Vector `data`
    * represents a point set. Only *the first* `numPoints` of `data` will
    * be added to the new series.
    *  @param data         Point set
    *  @param numPoints    Number of points to add
    *  @return Integer that represent the new point set’s position in the
    * JFreeChart `DefaultBoxAndWhiskerXYDataset` object.
    */
   public int add (double[] data, int numPoints) {
      DefaultBoxAndWhiskerCategoryDataset tempSeriesCollection =
            (DefaultBoxAndWhiskerCategoryDataset)seriesCollection;

      final List<Double> list = new ArrayList<Double>();
      for (int i = 0; i < numPoints; i ++)
         list.add(data[i]);

      int count = tempSeriesCollection.getColumnCount();
      tempSeriesCollection.add(list, 0, "Serie " + count);
      return count;
   }

   /**
    * Gets the current name of the selected series.
    *  @param series       series index.
    *  @return current name of the series.
    */
   public String getName (int series) {
      return (String)((DefaultBoxAndWhiskerCategoryDataset)seriesCollection).getColumnKey(series);
   }

   /**
    * Returns the range (@f$y@f$-coordinates) min and max values.
    *  @return range min and max values.
    */
   public double[] getRangeBounds() {
      double max=0, min=0;
      DefaultBoxAndWhiskerCategoryDataset tempSeriesCollection =
            (DefaultBoxAndWhiskerCategoryDataset)seriesCollection;

      if(tempSeriesCollection.getColumnCount() != 0 && tempSeriesCollection.getRowCount() != 0) {
         max = tempSeriesCollection.getItem(0, 0).getMaxOutlier().doubleValue() ;
         min = tempSeriesCollection.getItem(0, 0).getMinOutlier().doubleValue() ;
      }

      for(int i = 0; i < tempSeriesCollection.getRowCount(); i++) {
         for( int j = 0; j < tempSeriesCollection.getColumnCount(); j++) {
            max = Math.max(max, tempSeriesCollection.getItem(i, j).getMaxOutlier().doubleValue() );
            min = Math.min(min, tempSeriesCollection.getItem(i, j).getMinOutlier().doubleValue() );
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
      for(int i = 0; i < seriesCollection.getRowCount(); i++) {
         formatter.format(" Series " + i + " : %n");
         for(int j = 0; j < seriesCollection.getColumnCount(); j++)
            formatter.format(",%15e%n", seriesCollection.getValue(i, j));
      }
      return formatter.toString();
   }

   /**
    * NOT IMPLEMENTED: To do.
    *  @param ymin
    *  @param ymax
    *  @return LaTeX source code
    */
   public String toLatex (double YScale, double YShift, 
                          double ymin, double ymax) {
      throw new UnsupportedOperationException(" NOT implemented yet");
/*         
      // Calcule les bornes reelles du graphique, en prenant en compte la position des axes
      
      ymin = Math.min(YShift, ymin);
      ymax = Math.max(YShift, ymax);

      DefaultBoxAndWhiskerCategoryDataset tempSeriesCollection = (DefaultBoxAndWhiskerCategoryDataset)seriesCollection;
      Formatter formatter = new Formatter(Locale.US);
//       double var;
//       double margin = ((BoxAndWhiskerRenderer)renderer).getMargin();

//       for (int i = tempSeriesCollection.getColumnCount() - 1; i >= 0; i--) {
//          List temp = tempSeriesCollection.getBins(i);
//          ListIterator iter = temp.listIterator();
//          
//          Color color = (Color)renderer.getSeriesPaint(i);
//          String colorString = detectXColorClassic(color);
//          if (colorString == null) {
//             colorString = "color"+i;
//             formatter.format( "\\definecolor{%s}{rgb}{%.2f, %.2f, %.2f}%n",
//                               colorString, color.getRed()/255.0, color.getGreen()/255.0, color.getBlue()/255.0);
//          }
//          
//          HistogramBin currentBin=null;
//          while(iter.hasNext()) {
//             double currentMargin;
//             currentBin = (HistogramBin)iter.next();
//             currentMargin = ((margin*(currentBin.getEndBoundary()-currentBin.getStartBoundary()))-XShift)*XScale;
//             if ((currentBin.getStartBoundary() >= xmin && currentBin.getStartBoundary() <= xmax) 
//                && (currentBin.getCount() >= ymin && currentBin.getCount() <= ymax) )
//             {
//                var = Math.min( currentBin.getEndBoundary(), xmax);
//                if (filled[i]) {
//                   formatter.format("\\filldraw [line width=%.2fpt, opacity=%.2f, color=%s] ([xshift=%.4f] %.4f, %.4f) rectangle ([xshift=-%.4f] %.4f, %.4f); %%%n",
//                         lineWidth[i], (color.getAlpha()/255.0), colorString,
//                         currentMargin, (currentBin.getStartBoundary()-XShift)*XScale, 0.0,
//                         currentMargin, (var-XShift)*XScale, (currentBin.getCount()-YShift)*YScale);
//               }
//               else {
//                   formatter.format("\\draw [line width=%.2fpt, color=%s] ([xshift=%.4f] %.4f, %.4f) rectangle ([xshift=-%.4f] %.4f, %.4f); %%%n",
//                         lineWidth[i], colorString,
//                         currentMargin, (currentBin.getStartBoundary()-XShift)*XScale, 0.0,
//                         currentMargin, (var-XShift)*XScale, (currentBin.getCount()-YShift)*YScale);
//               }
//             }
//             else if (   (currentBin.getStartBoundary() >= xmin && currentBin.getStartBoundary() <= xmax) 
//                         && (currentBin.getCount() >= ymin && currentBin.getCount() > ymax) )
//             { // Cas ou notre rectangle ne peut pas etre affiche en entier (trop haut)
//                var = Math.min( currentBin.getEndBoundary(), xmax);
//                if (filled[i]) {
//                   formatter.format("\\filldraw [line width=%.2fpt,  opacity=%.2f, color=%s] ([xshift=%.4f] %.4f, %.4f) rectangle ([xshift=-%.4f] %.4f, %.4f); %%%n",
//                         lineWidth[i], (color.getAlpha()/255.0), colorString,
//                         currentMargin, (currentBin.getStartBoundary()-XShift)*XScale, 0.0,
//                         currentMargin, (var-XShift)*XScale, (ymax-YShift)*YScale);
//               formatter.format("\\draw [line width=%.2fpt, color=%s, style=dotted] ([xshift=%.4f] %.4f, %.4f) rectangle ([yshift=3mm, xshift=-%.4f] %.4f, %.4f); %%%n",
//                         lineWidth[i], colorString,
//                         currentMargin, (currentBin.getStartBoundary()-XShift)*XScale, (ymax-YShift)*YScale,
//                         currentMargin, (var-XShift)*XScale, (ymax-YShift)*YScale);
//                }
//                else {
//                   formatter.format("\\draw [line width=%.2fpt, color=%s] ([xshift=%.4f] %.4f, %.4f) rectangle ([xshift=-%.4f] %.4f, %.4f); %%%n",
//                         lineWidth[i], colorString,
//                         currentMargin, (currentBin.getStartBoundary()-XShift)*XScale, 0.0,
//                         currentMargin, (var-XShift)*XScale, (ymax-YShift)*YScale);
//                   
//               formatter.format("\\draw [line width=%.2fpt, color=%s, style=dotted] ([xshift=%.4f] %.4f, %.4f) rectangle ([yshift=3mm, xshift=-%.4f] %.4f, %.4f); %%%n",
//                         lineWidth[i], colorString,
//                         currentMargin, (currentBin.getStartBoundary()-XShift)*XScale, (ymax-YShift)*YScale,
//                         currentMargin, (var-XShift)*XScale, (ymax-YShift)*YScale);
//                }
//             }
//          }
//       }
      return formatter.toString();
*/
}

}

/**
 * @}
 */