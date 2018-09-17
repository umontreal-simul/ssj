/*
 * Class:        EmpiricalSeriesCollection
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

import   org.jfree.data.xy.XYSeriesCollection;
import   org.jfree.data.xy.XYSeries;
import   cern.colt.list.DoubleArrayList;

import   java.util.Locale;
import   java.util.Formatter;
import   java.util.List;
import   java.util.ListIterator;
import   java.awt.Color;

/**
 * Stores data used in a `EmpiricalChart`. `EmpiricalSeriesCollection`
 * provides complementary tools to draw empirical distribution charts, like
 * add plots series. This class is linked with `XYSeriesCollection` class
 * from JFreeChart to store data plots, and linked with JFreeChart
 * `EmpiricalRenderer` to render the plot. `EmpiricalRenderer` has been
 * developed at the Université de Montréal to extend the JFreeChart API, and
 * is used to render charts with an empirical chart style in a JFreeChart
 * chart.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class EmpiricalSeriesCollection extends SSJXYSeriesCollection {

   // for the zero point: fix its x a little smaller than the first point;
   // later when we know where x-axis begins, reset its x at the
   // beginning of x-axis; its y is always 0.
   private final double EMPIR_EPS = 0.015625;

   private String[] marksType;   //marks on points (+, x, *...)
   private String[] dashPattern; //line dashing  (solid, dotted, densely dotted, loosely dotted,
                                 //               dashed, densely dashed, loosely dashed)

   private double setZeroPoint (double x1) {
      // set temporary 0-point with x EPS smaller than x1
      return x1 - EMPIR_EPS*Math.abs(x1);
   }

   /**
    * Creates a new `EmpiricalSeriesCollection` instance with empty
    * dataset.
    */
   public EmpiricalSeriesCollection() {
      renderer = new EmpiricalRenderer();
      seriesCollection = new XYSeriesCollection();
   }

   /**
    * Creates a new `EmpiricalSeriesCollection` instance with default
    * parameters and given data series. Each input parameter represents an
    * observation set. The values of this observation set *must be sorted*
    * in increasing order.
    *  @param data         series of point sets.
    */
   public EmpiricalSeriesCollection (double[]... data) {
      seriesCollection = new XYSeriesCollection();
      renderer = new EmpiricalRenderer();
      XYSeriesCollection tempSeriesCollection = (XYSeriesCollection)seriesCollection;
      for (int j = 0; j < data.length; j++) {
         XYSeries serie = new XYSeries(" ");
         serie.add(setZeroPoint(data[j][0]), 0); // correct x-value of 0-point will be set later
         for (int k = 0; k < data[j].length; k++)
            serie.add(data[j][k], (double)(k + 1)/data[j].length);
         tempSeriesCollection.addSeries(serie);
      }

      /*set default colors*/
      for (int i = 0; i < tempSeriesCollection.getSeriesCount(); i++) {
         renderer.setSeriesPaint(i, getDefaultColor(i));
      }

      /* set default plot style*/
      marksType = new String[tempSeriesCollection.getSeriesCount()];
      dashPattern = new String[tempSeriesCollection.getSeriesCount()];
      for (int i = 0; i < tempSeriesCollection.getSeriesCount(); i++) {
         marksType[i] = "*";
         dashPattern[i] = "solid";
      }
   }

   /**
    * Creates a new `EmpiricalSeriesCollection` instance with default
    * parameters and a given series `data`. The values of `data` *must be
    * sorted* in increasing order. However, only *the first* `numPoints`
    * of `data` will be considered for the series.
    *  @param data         series of point sets
    *  @param numPoints    number of points to plot
    */
   public EmpiricalSeriesCollection (double[] data, int numPoints) {
      seriesCollection = new XYSeriesCollection();
      renderer = new EmpiricalRenderer();
      XYSeriesCollection tempSeriesCollection = (XYSeriesCollection)seriesCollection;
      XYSeries serie = new XYSeries(" ");
      serie.add(setZeroPoint(data[0]), 0); // correct x-value of zero point will be set later
      for (int k = 0; k < numPoints; k++)
         serie.add(data[k], (double)(k + 1) / numPoints);
      tempSeriesCollection.addSeries(serie);

      // set default colors
      renderer.setSeriesPaint(0, getDefaultColor(0));

      // set default plot style
      marksType = new String[1];
      dashPattern = new String[1];
      marksType[0] = "*";
      dashPattern[0] = "solid";
   }

   /**
    * Creates a new `EmpiricalSeriesCollection` instance with default
    * parameters and given data. The input parameter represents a
    * collection of data observation sets. Each <tt>DoubleArrayList</tt>
    * input parameter represents an observation set. The values of this
    * observation set *must be sorted in increasing order*. Each
    * DoubleArrayList variable corresponds to an observation set.
    *  @param data         series of point sets.
    */
   public EmpiricalSeriesCollection (DoubleArrayList... data) {
      seriesCollection = new XYSeriesCollection();
      renderer = new EmpiricalRenderer();
      XYSeriesCollection tempSeriesCollection = (XYSeriesCollection)seriesCollection;

      for (int j = 0; j < data.length; j++) {
         XYSeries serie = new XYSeries(" ");
         serie.add(setZeroPoint(data[j].get(0)), 0); // correct x-value of zero point will be set later
         for (int k = 0; k < data[j].size(); k++)
            serie.add(data[j].get(k), (double)(k + 1)/data[j].size());
         tempSeriesCollection.addSeries(serie);
      }

      /*set default colors*/
      for (int i = 0; i < tempSeriesCollection.getSeriesCount(); i++) {
         renderer.setSeriesPaint(i, getDefaultColor(i));
      }

      /* set default plot style*/
      marksType = new String[tempSeriesCollection.getSeriesCount()];
      dashPattern = new String[tempSeriesCollection.getSeriesCount()];
      for (int i = 0; i < tempSeriesCollection.getSeriesCount(); i++) {
         marksType[i] = "*";
         dashPattern[i] = "solid";
      }
   }

   /**
    * Creates a new `EmpiricalSeriesCollection` instance with default
    * parameters and given data. The input parameter represents a
    * collection of data observation sets. Each `TallyStore` input
    * parameter represents an observation set.
    *  @param tallies      series of point sets.
    */
   public EmpiricalSeriesCollection (TallyStore... tallies) {
      seriesCollection = new XYSeriesCollection();
      renderer = new EmpiricalRenderer();
      XYSeriesCollection tempSeriesCollection = (XYSeriesCollection)seriesCollection;

      for (int j = 0; j < tallies.length; j++) {
         TallyStore temp = tallies[j];
         temp.quickSort();
         double[] array = temp.getArray();
         XYSeries serie = new XYSeries(" ");
         serie.add(setZeroPoint(array[0]), 0); // correct x-value of zero point will be set later
         for (int k = 0; k < tallies[j].numberObs(); k++)
            serie.add(array[k], (double)(k + 1)/tallies[j].numberObs());
         tempSeriesCollection.addSeries(serie);
      }

      /*set default colors*/
      for (int i = 0; i < tempSeriesCollection.getSeriesCount(); i++) {
         renderer.setSeriesPaint(i, getDefaultColor(i));
      }

      /* set default plot style*/
      marksType = new String[tempSeriesCollection.getSeriesCount()];
      dashPattern = new String[tempSeriesCollection.getSeriesCount()];
      for (int i = 0; i < tempSeriesCollection.getSeriesCount(); i++) {
         marksType[i] = "*";
         dashPattern[i] = "solid";
      }
   }

   /**
    * Creates a new `EmpiricalSeriesCollection` instance with default
    * parameters and given data series. The input parameter represents a
    * set of plotting data. Each series of the given collection
    * corresponds to a plot on the chart.
    *  @param data         series of point sets.
    */
   public EmpiricalSeriesCollection (XYSeriesCollection data) {
      renderer = new EmpiricalRenderer();
      seriesCollection = data;

      /*set default colors*/
      for (int i = 0; i < data.getSeriesCount(); i++) {
         renderer.setSeriesPaint(i, getDefaultColor(i));
         XYSeries ser = data.getSeries(i);
         ser.add(setZeroPoint(ser.getX(0).doubleValue()), 0); // add zero point; will set its correct x-value later
     }

      /* set default plot style*/
      marksType = new String[data.getSeriesCount()];
      dashPattern = new String[data.getSeriesCount()];
      for (int i = 0; i < data.getSeriesCount(); i++) {
         marksType[i] = "*";
         dashPattern[i] = "solid";
      }
   }

   /**
    * @name Data control methods
    * @{
    */

   /**
    * Adds a data series into the series collection.
    *  @param observationSet new series values.
    *  @return Integer that represent the new point set’s position in the
    * `XYSeriesCollection` object.
    */
   public int add (double[] observationSet) {
      return add(observationSet, observationSet.length);
   }

   /**
    * Adds a data series into the series collection. Only *the first*
    * `numPoints` of `observationSet` will be added to the new series.
    *  @param observationSet new series values.
    *  @param numPoints    number of points to add.
    *  @return Integer that represent the new point set’s position in the
    * `XYSeriesCollection` object.
    */
   public int add (double[] observationSet, int numPoints) {
      XYSeriesCollection tempSeriesCollection = (XYSeriesCollection)seriesCollection;

      XYSeries serie = new XYSeries(" ");
      serie.add(setZeroPoint(observationSet[0]), 0); // coordinates of first point will be reset later
      for (int k = 0; k < numPoints; k++)
         serie.add(observationSet[k], (double)(k + 1)/numPoints);
      tempSeriesCollection.addSeries(serie);

      // color
      int j = seriesCollection.getSeriesCount()-1;
      renderer.setSeriesPaint(j, getDefaultColor(j));

      String[] newMarksType = new String[seriesCollection.getSeriesCount()];
      String[] newDashPattern = new String[seriesCollection.getSeriesCount()];
      for (j = 0; j < seriesCollection.getSeriesCount()-1; j++) {
         newMarksType[j] = marksType[j];
         newDashPattern[j] = dashPattern[j];
      }

      newMarksType[j] = "*";
      newDashPattern[j] = "solid";
      marksType = newMarksType;
      dashPattern = newDashPattern;

      return seriesCollection.getSeriesCount()-1;
   }

   /**
    * Adds a data series into the series collection.
    *  @param observationSet new series values.
    *  @return Integer that represent the new point set’s position in the
    * `XYSeriesCollection` object.
    */
   public int add (DoubleArrayList observationSet) {
      return add(observationSet.elements(), observationSet.size());
   }

   /**
    * Adds a data series into the series collection.
    *  @param tally         @ref umontreal.ssj.stat.TallyStore to add
    *                       values.
    *  @return Integer that represent the new point set’s position in the
    * `XYSeriesCollection` object.
    */
   public int add (TallyStore tally) {
      tally.quickSort();
      return add(tally.getArray(), tally.numberObs());
   }

   /**
    * @}
    */

   /**
    * @name Rendering methods
    * @{
    */

   /**
    * Returns the mark type associated with the <tt>series</tt>-th data
    * series.
    *  @param series       series index.
    *  @return mark type.
    */
   public String getMarksType (int series) {
      return marksType[series];
   }

   /**
    * Adds marks on points to a data series. It is possible to use all the
    * marks provided by the TikZ package, some of which are `*`, `+` and
    * `x`. A blank character, used by default, will disable marks. The
    * PGF/TikZ documentation provides more information about placing marks
    * on plots.
    *  @param series       series index.
    *  @param marksType    mark type.
    */
   public void setMarksType (int series, String marksType) {
      this.marksType[series] = marksType;
   }

   /**
    * Returns the dash pattern associated with the <tt>series</tt>-th data
    * series.
    *  @param series       series index.
    *  @return dash style.
    */
   public String getDashPattern (int series) {
      return dashPattern[series];
   }

   /**
    * Selects dash pattern for a data series. It is possible to use all
    * the dash options provided by the TikZ package: solid, dotted,
    * densely dotted, loosely dotted, dashed, densely dashed and loosely
    * dashed.
    *  @param series       series index.
    *  @param dashPattern  dash style.
    */
   public void setDashPattern (int series, String dashPattern) {
      this.dashPattern[series] = dashPattern;
   }
   public String toLatex (double XScale, double YScale,
                          double XShift, double YShift,
                          double xmin, double xmax,
                          double ymin, double ymax) {

      // Calcule les bornes reelles du graphique, en prenant en compte la position des axes
      xmin = Math.min(XShift, xmin);
      xmax = Math.max(XShift, xmax);
      ymin = Math.min(YShift, ymin);
      ymax = Math.max(YShift, ymax);

      XYSeriesCollection tempSeriesCollection = (XYSeriesCollection)seriesCollection;
      Formatter formatter = new Formatter(Locale.US);
      double var;
      for (int i = tempSeriesCollection.getSeriesCount()-1; i >= 0; i--) {

         Color color = (Color)renderer.getSeriesPaint(i);
         String colorString = detectXColorClassic(color);
         if (colorString == null) {
            colorString = "color"+i;
            formatter.format( "\\definecolor{%s}{rgb}{%.2f, %.2f, %.2f}%n",
                              colorString, color.getRed()/255.0, color.getGreen()/255.0, color.getBlue()/255.0);
         }

         double currentX, currentY, nextX, nextY;
         for (int j = 0; j < tempSeriesCollection.getItemCount(i); j++) {
            currentX = tempSeriesCollection.getXValue(i, j);
            if (j == tempSeriesCollection.getItemCount(i)-1)
               nextX = xmax;
            else
               nextX = Math.min(tempSeriesCollection.getXValue(i, j+1), xmax);

            if ((currentX >= xmin && currentX <= xmax) )
            {
               formatter.format("\\draw [color=%s] plot[mark=%s] (%.4f, %.4f) --plot[style=%s] (%.4f, %.4f); %%%n",
                        colorString, marksType[i],
                        (currentX-XShift)*XScale, (tempSeriesCollection.getYValue(i, j)-YShift)*YScale,
                        dashPattern[i],
                        (nextX-XShift)*XScale, (tempSeriesCollection.getYValue(i, j)-YShift)*YScale);
            }
         }
      }
      return formatter.toString();
   }
}

/**
 * @}
 */