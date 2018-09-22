
package umontreal.ssj.charts;

import java.util.*;
import org.jfree.data.general.DatasetChangeEvent;
import org.jfree.data.statistics.HistogramBin;
import org.jfree.data.statistics.HistogramType;
import org.jfree.data.xy.AbstractIntervalXYDataset;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.util.ObjectUtilities;
import org.jfree.util.PublicCloneable;

/*
 * Class:        CustomHistogramDataset
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


/**
 * A dataset that can be used for creating histograms.
 * This class is inspired from JFreeChart HistogramDataset class,
 * and provides tools to customize histogram bins.
 *
 */
public class CustomHistogramDataset extends AbstractIntervalXYDataset
         implements IntervalXYDataset, Cloneable, PublicCloneable {

   /** A list of maps. */
   private List list;

   /** The histogram type. */
   public HistogramType type;

   /**
    * Creates a new (empty) dataset with a default type of
    * `HistogramType.FREQUENCY`.
    */
   public CustomHistogramDataset()
   {
      list = new ArrayList();
      type = HistogramType.FREQUENCY;
   }

   /**
    * Returns the histogram type.
    */
   public HistogramType getType()
   {
      return type;
   }

   /**
    * Sets the histogram type and sends a @ref DatasetChangeEvent to all
    * registered listeners.
    *
    * @param type  the type (<code>null</code> not permitted).
    */
   public void setType(HistogramType type)
   {
      if (type == null) {
         throw new IllegalArgumentException("Null 'type' argument");
      } else {
         this.type = type;
         notifyListeners(new DatasetChangeEvent(this, this));
         return ;
      }
   }

   /**
    * Adds a series to the dataset, using the specified number of bins.
    *
    * @param key  the series key (<code>null</code> not permitted).
    * @param values the values (<code>null</code> not permitted).
    * @param bins  the number of bins (must be at least 1).
    */
   public void addSeries(Comparable key, double values[], int bins)
   {
      double minimum = getMinimum(values);
      double maximum = getMaximum(values);
      addSeries(key, values, bins, minimum, maximum);
   }

   /**
    * Adds a series to the dataset, using the specified number of bins.
    *
    * @param key  the series key (<code>null</code> not permitted).
    * @param values the values (<code>null</code> not permitted).
    * @param numPoints  only the first numPoints values are used.
    * @param bins  the number of bins (must be at least 1).
    */
   public void addSeries(Comparable key, double values[], int numPoints, int bins)
   {
      double minimum = getMinimum(values);
      double maximum = getMaximum(values);
      addSeries(key, values, numPoints, bins, minimum, maximum);
   }

   /**
    * Adds a series to the dataset. Any data value less than minimum will be
    * assigned to the first bin, and any data value greater than maximum will
    * be assigned to the last bin.  Values falling on the boundary of
    * adjacent bins will be assigned to the higher indexed bin.
    *
    * @param key  the series key (<code>null</code> not permitted).
    * @param values  the raw observations.
    * @param bins  the number of bins (must be at least 1).
    * @param minimum  the lower bound of the bin range.
    * @param maximum  the upper bound of the bin range.
    */
   public void addSeries(Comparable key, double values[], int bins,
                         double minimum, double maximum)
   {
      addSeries(key, values, values.length, bins, minimum, maximum);
   }

   /**
    * Adds a series to the dataset. Any data value less than minimum will be
    * assigned to the first bin, and any data value greater than maximum will
    * be assigned to the last bin.  Values falling on the boundary of
    * adjacent bins will be assigned to the higher indexed bin.
    * Only the first numPoints values are used.
    *
    * @param key  the series key (<code>null</code> not permitted).
    * @param values  the raw observations.
    * @param numPoints  only the first numPoints values are used.
    * @param bins  the number of bins (must be at least 1).
    * @param minimum  the lower bound of the bin range.
    * @param maximum  the upper bound of the bin range.
    */
   public void addSeries(Comparable key, double values[], int numPoints,
                         int bins, double minimum, double maximum)
   {
      if (key == null)
         throw new IllegalArgumentException("Null 'key' argument.");
      if (values == null)
         throw new IllegalArgumentException("Null 'values' argument.");
      if (bins < 1)
         throw new IllegalArgumentException("The 'bins' value must be at least 1.");
      double binWidth = (maximum - minimum) / (double)bins;
      double lower = minimum;
      final double EPS = 1.0e-15;
      List binList = new ArrayList(bins);
      for (int i = 0; i < bins; i++) {
         HistogramBin bin;
         if (i == bins - 1) {
            bin = new HistogramBin(lower, maximum*(1.0 + EPS));
         } else {
            double upper = minimum + (double)(i + 1) * binWidth;
            bin = new HistogramBin(lower, upper);
            lower = upper;
         }
         binList.add(bin);
      }

      ArrayList valuesList = new ArrayList(numPoints);
      for (int i = 0; i < numPoints; i++)
         valuesList.add(new Double(values[i]));

      synchronizeValuesAndBins(binList, valuesList);
      Map map = new HashMap();
      map.put("key", key);
      map.put("values", valuesList);
      map.put("bins", binList);
      map.put("numPoints", new Integer(numPoints));
      map.put("bin width", new Double(binWidth));
      list.add(map);
   }

   /**
    * Adds a series to the dataset. Values falling on the boundary of
    * adjacent bins will be assigned to the higher indexed bin.
    *
    * @param key  the series key (<code>null</code> not permitted).
    * @param values  the raw observations.
    * @param bins  new bins (size must be at least 1).
    */
   public void addSeries(Comparable key, double values[], HistogramBin bins[])
   {
      addSeries(key, values, values.length, bins);
   }

   /**
    * Adds a series to the dataset. Values falling on the boundary of
    * adjacent bins will be assigned to the higher indexed bin.
    *
    * @param key  the series key (<code>null</code> not permitted).
    * @param values  the raw observations.
    * @param numPoints  only the first numPoints values are used.
    * @param bins  new bins (size must be at least 1).
    */
   public void addSeries(Comparable key, double values[], int numPoints,
                         HistogramBin bins[])
   {
      if (key == null)
         throw new IllegalArgumentException("Null 'key' argument.");
      if (values == null)
         throw new IllegalArgumentException("Null 'values' argument.");
      if (bins == null || bins.length < 2)
         throw new IllegalArgumentException
             ("The 'bins' table must contain at least 1 org.jfree.data.statistics.HistogramBin.");
      List binList = new ArrayList(bins.length);
      for (int i = 0; i < bins.length; i++)
         binList.add(bins[i]);

      ArrayList valuesList = new ArrayList(numPoints);
      for (int i = 0; i < numPoints; i++)
         valuesList.add(new Double(values[i]));

      synchronizeValuesAndBins(binList, valuesList);
      Map map = new HashMap();
      map.put("key", key);
      map.put("values", valuesList);
      map.put("bins", binList);
      map.put("numPoints", new Integer(numPoints));
      map.put("bin width", new Double( -1D));
      list.add(map);
   }

   /**
    * Returns the minimum value in an array of values.
    *
    * @param values  the values (<code>null</code> not permitted and
    *                zero-length array not permitted).
    *
    * @return The minimum value.
    */
   private double getMinimum(double values[])
   {
      if (values == null || values.length < 1)
         throw new IllegalArgumentException
            ("Null or zero length 'values' argument.");
      double min = 1.7E+308D;
      for (int i = 0; i < values.length; i++)
         if (values[i] < min)
            min = values[i];

      return min;
   }

   /**
    * Returns the maximum value in an array of values.
    *
    * @param values  the values (<code>null</code> not permitted and
    *                zero-length array not permitted).
    *
    * @return The maximum value.
    */
   private double getMaximum(double values[])
   {
      if (values == null || values.length < 1)
         throw new IllegalArgumentException
            ("Null or zero length 'values' argument.");
      double max = -1.7E+308D;
      for (int i = 0; i < values.length; i++)
         if (values[i] > max)
            max = values[i];

      return max;
   }

   /**
    * Returns the bins for a series.
    *
    * @param series  the series index (in the range <code>0</code> to
    *     <code>getSeriesCount() - 1</code>).
    *
    * @return A list of bins.
    *
    * @throws IndexOutOfBoundsException if <code>series</code> is outside the
    *     specified range.
    */
   public List getBins(int series)
   {
      Map map = (Map)list.get(series);
      return (List)map.get("bins");
   }

   /**
    * Sets the bins for a series.
    *
    * @param series  the series index (in the range <code>0</code> to
    *     <code>getSeriesCount() - 1</code>).
    * @param bins  the number of bins (must be at least 1).
    *
    * @throws IndexOutOfBoundsException if <code>series</code> is outside the
    *     specified range.
    */
   public void setBins(int series, int bins)
   {
      double minimum = getMinimum(getValues(series));
      double maximum = getMaximum(getValues(series));
      setBins(series, bins, minimum, maximum);
   }

   /**
    * Sets the bins for a series.
    *
    * @param series  the series index (in the range <code>0</code> to
    *     <code>getSeriesCount() - 1</code>).
    * @param bins  the number of bins (must be at least 1).
    * @param minimum  the lower bound of the bin range.
    * @param maximum  the upper bound of the bin range.
    *
    * @throws IndexOutOfBoundsException if <code>series</code> is outside the
    *     specified range.
    */
   public void setBins(int series, int bins, double minimum, double maximum)
   {
      Map map = (Map)list.get(series);
      List currentValues = (List)map.get("values");
      double binWidth = (maximum - minimum) / (double)bins;
      double lower = minimum;
      List binList = new ArrayList(bins);
      double EPS = 1.0e-15;
      for (int i = 0; i < bins; i++) {
         HistogramBin bin;
         if (i == bins - 1) {
            bin = new HistogramBin(lower, maximum*(1.0 + EPS));
         } else {
            double upper = minimum + (double)(i + 1) * binWidth;
            bin = new HistogramBin(lower, upper);
            lower = upper;
         }
         binList.add(bin);
      }

      synchronizeValuesAndBins(binList, currentValues);
      map.put("values", currentValues);
      map.put("bins", binList);
   }

   /**
    * Sets the bins for a series.
    *
    * @param series  the series index (in the range <code>0</code> to
    *     <code>getSeriesCount() - 1</code>).
    * @param bins  the number of bins (must be at least 1).
    *
    * @throws IndexOutOfBoundsException if <code>series</code> is outside the
    *     specified range.
    */
   public void setBins(int series, HistogramBin bins[])
   {
      Map map = (Map)list.get(series);
      List currentValues = (List)map.get("values");
      ArrayList binList = new ArrayList(bins.length);
      for (int i = 0; i < bins.length; i++)
         binList.add(bins[i]);

      synchronizeValuesAndBins(binList, currentValues);
      map.put("values", currentValues);
      map.put("bins", binList);
   }

   /**
    * Returns the values for a series.
    *
    * @param series  the series index (in the range <code>0</code> to
    *     <code>getSeriesCount() - 1</code>).
    *
    * @return A list of values.
    *
    * @throws IndexOutOfBoundsException if <code>series</code> is outside the
    *     specified range.
    */
   public List getValuesList(int series)
   {
      Map map = (Map)list.get(series);
      return (List)map.get("values");
   }

   /**
    * Returns the values for a series.
    *
    * @param series  the series index (in the range <code>0</code> to
    *     <code>getSeriesCount() - 1</code>).
    *
    * @return A table of values.
    *
    * @throws IndexOutOfBoundsException if <code>series</code> is outside the
    *     specified range.
    */
   public double[] getValues(int series)
   {
      List valuesList = (List)((Map)list.get(series)).get("values");
      ListIterator iter = valuesList.listIterator();
      double retour[] = new double[valuesList.size()];
      for (int i = 0; iter.hasNext(); i++)
         retour[i] = ((Double)iter.next()).doubleValue();

      return retour;
   }

   /**
    * Sets the values for a series.
    *
    * @param series  the series index (in the range <code>0</code> to
    *     <code>getSeriesCount() - 1</code>).
    * @param valuesList  List of new values.
    *
    * @throws IndexOutOfBoundsException if <code>series</code> is outside the
    *     specified range.
    */
   public void setValues(int series, List valuesList)
   {
      Map map = (Map)list.get(series);
      List currentBins = (List)map.get("bins");
      synchronizeValuesAndBins(currentBins, valuesList);
      map.put("values", valuesList);
      map.put("bins", currentBins);
   }

   /**
    * Sets the values for a series.
    *
    * @param series  the series index (in the range <code>0</code> to
    *     <code>getSeriesCount() - 1</code>).
    * @param values  Table of new values.
    *
    * @throws IndexOutOfBoundsException if <code>series</code> is outside the
    *     specified range.
    */
   public void setValues(int series, double values[])
   {
      ArrayList valuesList = new ArrayList(values.length);
      for (int i = 0; i < values.length; i++)
         valuesList.add(new Double(values[i]));

      setValues(series, ((List) (valuesList)));
   }

   /**
    * Synchronize values to bins. Compute bins values.
    *
    * @param bins   List of bins.
    * @param values List of values.
    */
   private void synchronizeValuesAndBins(List bins, List values)
   {
      ListIterator iterBins = bins.listIterator(0);
      ListIterator iterValues = values.listIterator();
      HistogramBin bin;
      for (; iterBins.hasNext(); iterBins.set(
            new HistogramBin(bin.getStartBoundary(), bin.getEndBoundary())))
         bin = (HistogramBin)iterBins.next();

      iterBins = bins.listIterator(0);
      while (iterValues.hasNext()) {
         double currentValue = ((Double)iterValues.next()).doubleValue();
         boolean continu = true;
         iterBins = bins.listIterator(0);
         while (continu && iterBins.hasNext()) {
            HistogramBin tempBin = (HistogramBin)iterBins.next();
            if (currentValue >= tempBin.getStartBoundary() &&
                currentValue <  tempBin.getEndBoundary()) {
               tempBin.incrementCount();
               continu = false;
            }
         }
      }
   }

   /**
    * Returns the total number of observations for a series.
    *
    * @param series  the series index.
    *
    * @return The total.
    */
   public int getTotal(int series)
   {
      Map map = (Map)list.get(series);
      return ((Integer)map.get("numPoints")).intValue();
   }

   /**
    * Returns the bin width for a series.
    *
    * @param series  the series index (zero based).
    *
    * @return The bin width.
    */
   public double getBinWidth(int series)
   {
      Map map = (Map)list.get(series);
      return ((Double)map.get("bin width")).doubleValue();
   }

   /**
    * Returns the number of series in the dataset.
    *
    * @return The series count.
    */
   public int getSeriesCount()
   {
      return list.size();
   }

   /**
    * Returns the key for a series.
    *
    * @param series  the series index (in the range <code>0</code> to
    *     <code>getSeriesCount() - 1</code>).
    *
    * @return The series key.
    *
    * @throws IndexOutOfBoundsException if <code>series</code> is outside the
    *     specified range.
    */
   public Comparable getSeriesKey(int series)
   {
      Map map = (Map)list.get(series);
      return (Comparable)map.get("key");
   }

   /**
    * Returns the number of data items for a series.
    *
    * @param series  the series index (in the range <code>0</code> to
    *     <code>getSeriesCount() - 1</code>).
    *
    * @return The item count.
    *
    * @throws IndexOutOfBoundsException if <code>series</code> is outside the
    *     specified range.
    */
   public int getItemCount(int series)
   {
      return getBins(series).size();
   }

   /**
    * Returns the X value for a bin.  This value won't be used for plotting
    * histograms, since the renderer will ignore it.  But other renderers can
    * use it (for example, you could use the dataset to create a line
    * chart).
    *
    * @param series  the series index (in the range <code>0</code> to
    *     <code>getSeriesCount() - 1</code>).
    * @param item  the item index (zero based).
    *
    * @return The start value.
    *
    * @throws IndexOutOfBoundsException if <code>series</code> is outside the
    *     specified range.
    */
   public Number getX(int series, int item)
   {
      List bins = getBins(series);
      HistogramBin bin = (HistogramBin)bins.get(item);
      double x = (bin.getStartBoundary() + bin.getEndBoundary()) / 2D;
      return new Double(x);
   }

   /**
    * Returns the y-value for a bin (calculated to take into account the
    * histogram type).
    *
    * @param series  the series index (in the range <code>0</code> to
    *     <code>getSeriesCount() - 1</code>).
    * @param item  the item index (zero based).
    *
    * @return The y-value.
    *
    * @throws IndexOutOfBoundsException if <code>series</code> is outside the
    *     specified range.
    */
   public Number getY(int series, int item)
   {
      List bins = getBins(series);
      HistogramBin bin = (HistogramBin)bins.get(item);
      double total = getTotal(series);
      double binWidth = getBinWidth(series);
      if (type == HistogramType.FREQUENCY)
         return new Double(bin.getCount());
      if (type == HistogramType.RELATIVE_FREQUENCY)
         return new Double((double)bin.getCount() / total);
      if (type == HistogramType.SCALE_AREA_TO_1)
         return new Double((double)bin.getCount() / (binWidth * total));
      else
         throw new IllegalStateException();
   }

   /**
    * Returns the start value for a bin.
    *
    * @param series  the series index (in the range <code>0</code> to
    *     <code>getSeriesCount() - 1</code>).
    * @param item  the item index (zero based).
    *
    * @return The start value.
    *
    * @throws IndexOutOfBoundsException if <code>series</code> is outside the
    *     specified range.
    */
   public Number getStartX(int series, int item)
   {
      List bins = getBins(series);
      HistogramBin bin = (HistogramBin)bins.get(item);
      return new Double(bin.getStartBoundary());
   }

   /**
    * Returns the end value for a bin.
    *
    * @param series  the series index (in the range <code>0</code> to
    *     <code>getSeriesCount() - 1</code>).
    * @param item  the item index (zero based).
    *
    * @return The end value.
    *
    * @throws IndexOutOfBoundsException if <code>series</code> is outside the
    *     specified range.
    */
   public Number getEndX(int series, int item)
   {
      List bins = getBins(series);
      HistogramBin bin = (HistogramBin)bins.get(item);
      return new Double(bin.getEndBoundary());
   }

   /**
    * Returns the start y-value for a bin (which is the same as the y-value).
    * This method exists only to support the general form of the
    * {IntervalXYDataset} interface.
    *
    * @param series  the series index (in the range <code>0</code> to
    *     <code>getSeriesCount() - 1</code>).
    * @param item  the item index (zero based).
    *
    * @return The y-value.
    *
    * @throws IndexOutOfBoundsException if <code>series</code> is outside the
    *     specified range.
    */
   public Number getStartY(int series, int item)
   {
      return getY(series, item);
   }

   /**
    * Returns the end y-value for a bin (which is the same as the y-value).
    * This method exists only to support the general form of the
    * IntervalXYDataset} interface.
    *
    * @param series  the series index (in the range <code>0</code> to
    *     <code>getSeriesCount() - 1</code>).
    * @param item  the item index (zero based).
    *
    * @return The Y value.
    *
    * @throws IndexOutOfBoundsException if <code>series</code> is outside the
    *     specified range.
    */
   public Number getEndY(int series, int item)
   {
      return getY(series, item);
   }

   /**
    * Tests this dataset for equality with an arbitrary object.
    *
    * @param obj  the object to test against (<code>null</code> permitted).
    *
    * @return A boolean.
    */
   public boolean equals(Object obj)
   {
      if (obj == this)
         return true;
      if (!(obj instanceof CustomHistogramDataset))
         return false;
      CustomHistogramDataset that = (CustomHistogramDataset)obj;
      if (!ObjectUtilities.equal(type, that.type))
         return false;
      return ObjectUtilities.equal(list, that.list);
   }

   /**
    * Returns a clone of the dataset.
    *
    * @return A clone of the dataset.
    *
    * @throws CloneNotSupportedException if the object cannot be cloned.
    */
   public Object clone()
   throws CloneNotSupportedException
   {
      return super.clone();
   }
}
