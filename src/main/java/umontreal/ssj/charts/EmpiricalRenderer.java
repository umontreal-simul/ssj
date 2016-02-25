package umontreal.ssj.charts;

import java.awt.Graphics2D;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.Shape;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.entity.EntityCollection;
import org.jfree.chart.entity.XYItemEntity;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.chart.urls.XYURLGenerator;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRendererState;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.plot.PlotRenderingInfo;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.CrosshairState;
import org.jfree.data.xy.XYDataset;
import org.jfree.util.PublicCloneable;
import org.jfree.util.ShapeUtilities;

/*
 * Class:        EmpiricalRenderer
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
 * A renderer that draws horizontal lines between points and/or draws shapes
 * at each data point to provide an empirical style chart.  This renderer is
 * designed for use with the {XYPlot} class.
 */
public class EmpiricalRenderer extends XYLineAndShapeRenderer
         implements XYItemRenderer, Cloneable, PublicCloneable
{

   /**
    * Creates a new renderer.
    */
   public EmpiricalRenderer()
   {
      this(null, null);
   }


   /**
    * Creates a new renderer with selected tool tip and url generators.
    *
    * @param  toolTipGenerator   Tool tip generator.
    * @param  urlGenerator       Url generator.
    */
   public EmpiricalRenderer(XYToolTipGenerator toolTipGenerator, XYURLGenerator urlGenerator)
   {
      setBaseToolTipGenerator(toolTipGenerator);
      setURLGenerator(urlGenerator);
      setShapesFilled(true);
      setShapesVisible(true);
   }

   /**
    * Draws the visual representation of a single data item.
    *
    * @param g2           the graphics device.
    * @param state        the renderer state.
    * @param dataArea     the area within which the data is being drawn.
    * @param info         collects information about the drawing.
    * @param plot         the plot (can be used to obtain standard color
    *                     information etc).
    * @param domainAxis   the domain axis.
    * @param rangeAxis    the range axis.
    * @param dataset      the dataset.
    * @param series       the series index (zero-based).
    * @param item         the item index (zero-based).
    * @param crosshairState  crosshair information for the plot 
    *                        (<code>null</code> permitted).
    * @param pass         the pass index.
    */
   public void drawItem(Graphics2D g2, XYItemRendererState state, Rectangle2D dataArea,
                        PlotRenderingInfo info, XYPlot plot, ValueAxis domainAxis, ValueAxis rangeAxis,
                        XYDataset dataset, int series, int item, CrosshairState crosshairState, int pass)
   {

      if (!getItemVisible(series, item))
         return ;
      PlotOrientation orientation = plot.getOrientation();
      java.awt.Paint seriesPaint = getItemPaint(series, item);
      java.awt.Stroke seriesStroke = getItemStroke(series, item);
      g2.setPaint(seriesPaint);
      g2.setStroke(seriesStroke);
      double x0 = dataset.getXValue(series, item);
      double y0 = dataset.getYValue(series, item);
      if (java.lang.Double.isNaN(y0))
         return ;
      org.jfree.ui.RectangleEdge xAxisLocation = plot.getDomainAxisEdge();
      org.jfree.ui.RectangleEdge yAxisLocation = plot.getRangeAxisEdge();
      double transX0 = domainAxis.valueToJava2D(x0, dataArea, xAxisLocation);
      double transY0 = rangeAxis.valueToJava2D(y0, dataArea, yAxisLocation);

      double x1 = 0, y1 = 0;
      if (item < dataset.getItemCount(series) - 1) {
         x1 = dataset.getXValue(series, item + 1);
         y1 = dataset.getYValue(series, item + 1);
      } else {
         x1 = dataArea.getMaxX();
         y1 = dataArea.getMaxY();
      }

      boolean useFillPaint = getUseFillPaint();
      ;
      boolean drawOutlines = getDrawOutlines();
      if (!java.lang.Double.isNaN(y0)) {
         double transX1;
         double transY1;
         if (item < dataset.getItemCount(series) - 1) {
            transX1 = domainAxis.valueToJava2D(x1, dataArea, xAxisLocation);
            transY1 = rangeAxis.valueToJava2D(y1, dataArea, yAxisLocation);
         } else {
            transX1 = x1;
            transY1 = y1;
         }
         Line2D line = state.workingLine;
         if (orientation == PlotOrientation.HORIZONTAL) {
            line.setLine(transY0, transX0, transY0, transX1);
            g2.draw(line);
         } else if (orientation == PlotOrientation.VERTICAL) {
            line.setLine(transX0, transY0, transX1, transY0);
            g2.draw(line);
         }
      }
      if (getItemShapeVisible(series, item)) {
         Shape shape = getItemShape(series, item);
         if (orientation == PlotOrientation.HORIZONTAL)
            shape = ShapeUtilities.createTranslatedShape(shape, transY0, transX0);
         else if (orientation == PlotOrientation.VERTICAL)
            shape = ShapeUtilities.createTranslatedShape(shape, transX0, transY0);
         if (shape.intersects(dataArea)) {
            if (getItemShapeFilled(series, item)) {
               if (useFillPaint)
                  g2.setPaint(getItemFillPaint(series, item));
               else
                  g2.setPaint(getItemPaint(series, item));
               g2.fill(shape);
            }
            if (drawOutlines) {
               if (getUseOutlinePaint())
                  g2.setPaint(getItemOutlinePaint(series, item));
               else
                  g2.setPaint(getItemPaint(series, item));
               g2.setStroke(getItemOutlineStroke(series, item));
               g2.draw(shape);
            }
         }
      }
      if (isItemLabelVisible(series, item)) {
         double xx = transX0;
         double yy = transY0;
         if (orientation == PlotOrientation.HORIZONTAL) {
            xx = transY0;
            yy = transX0;
         }
         drawItemLabel(g2, orientation, dataset, series, item, xx, yy, y0 < 0.0D);
      }
      int domainAxisIndex = plot.getDomainAxisIndex(domainAxis);
      int rangeAxisIndex = plot.getRangeAxisIndex(rangeAxis);
      updateCrosshairValues(crosshairState, x0, y0, domainAxisIndex, rangeAxisIndex, transX0, transY0, orientation);
      if (state.getInfo() != null) {
         EntityCollection entities = state.getEntityCollection();
         if (entities != null) {
            int r = getDefaultEntityRadius();
            java.awt.Shape shape = orientation != PlotOrientation.VERTICAL ? ((java.awt.Shape) (new java.awt.geom.Rectangle2D.Double(transY0 - (double)r, transX0 - (double)r, 2 * r, 2 * r))) : ((java.awt.Shape) (new java.awt.geom.Rectangle2D.Double(transX0 - (double)r, transY0 - (double)r, 2 * r, 2 * r)));
            if (shape != null) {
               String tip = null;
               XYToolTipGenerator generator = getToolTipGenerator(series, item);
               if (generator != null)
                  tip = generator.generateToolTip(dataset, series, item);
               String url = null;
               if (getURLGenerator() != null)
                  url = getURLGenerator().generateURL(dataset, series, item);
               XYItemEntity entity = new XYItemEntity(shape, dataset, series, item, tip, url);
               entities.add(entity);
            }
         }
      }
   }

   /**
    * Returns a clone of the renderer.
    *
    * @return A clone.
    *
    * @throws CloneNotSupportedException if the clone cannot be created.
    */
   public Object clone() throws CloneNotSupportedException
   {
      return super.clone();
   }
}
