/*
 * Class:        Axis
 * Description:  An axis of a chart and its properties
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

import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import java.util.Locale;
import java.util.Formatter;
import java.lang.Math;

/**
 * Represents an axis of a chart encapsulated by an instance of `XYChart`.
 * `Axis` uses the JFreeChart class `NumberAxis` to store some axis
 * properties. This class represents the @f$x@f$-axis or the @f$y@f$-axis of
 * a `XYChart` and, consequently, is drawn when calling the  #toLatex method.
 * It provides tools to customize the axis in modifying labels and
 * description.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class Axis {

   public static final boolean ORIENTATION_VERTICAL = false;
   public static final boolean ORIENTATION_HORIZONTAL = true;

   private NumberAxis axis;        // Uses JFreeChart API
   private boolean orientation;    // Vertical or horizontal
   private double twinAxisPosition = 0;

   private boolean tick0Flag; // if true, label number at origin is shifted
   // right or above a little to prevent the label from being on the
   // perpendicular axis.

   private boolean  labelsFlag;    // true : enable manual labels
   private String[] labelsName;    // Sets manual labels
   private double[] labelsValue;   // Where putting the labels on the axes

   private double fixStep (double minA, double maxA, double step) {
      // reset step for labels on axes if it is too large or too small
      final double H = 10;
      final double inter = maxA - minA;
      while (step > inter)   // step is too large, make it smaller
         step /= H;
      while (step < inter / H)   // step is too small, make it larger
         step *= H;
      if (step > inter/2)
         step /= 2;
      return step;
   }

   /**
    * Create a new `Axis` instance from an existing `NumberAxis` instance
    * with vertical (@f$y@f$-axis) or horizontal (@f$x@f$-axis)
    * orientation.
    *  @param inAxis       NumberAxis instance associated to the new
    *                      variable.
    *  @param orientation  axis direction, horizontal or vertical
    */
   public Axis (NumberAxis inAxis, boolean orientation) {
      this.axis = inAxis;
      this.orientation = orientation;
      this.labelsFlag = false;
      this.labelsValue = null;
      this.labelsName = null;
      inAxis.setAutoRangeIncludesZero(false);
      inAxis.setLowerMargin(0);
      inAxis.setUpperMargin(0);
      axis.setTickUnit(new NumberTickUnit(computeAutoTickValue()));
   }

   /**
    * Returns the `NumberAxis` instance (from JFreeChart) linked with the
    * current variable.
    *  @return the `NumberAxis` instance (from JFreeChart) linked with the
    * current variable.
    */
   protected NumberAxis getAxis() {
      return axis;
   }

   /**
    * Returns the axis description.
    *  @return the axis description.
    */
   public String getLabel() {
      return axis.getLabel();
   }

   /**
    * Sets the axis description. This description will be displayed on the
    * chart, near the axis.
    *  @param label        axis label.
    */
   public void setLabel (String label) {
      axis.setLabel(label);
   }

   /**
    * Sets a periodic label display. Labels will be shown every `tick`
    * unit. This tick unit replaces the default unit.
    *  @param tick         tick unit.
    */
   public void setLabels (double tick) {
      axis.setTickUnit(new NumberTickUnit(tick));
      labelsFlag = false;
   }

   /**
    * Calculates and sets an automatic tick unit.
    */
   public void setLabelsAuto() {
      axis.setTickUnit(new NumberTickUnit(computeAutoTickValue()));
      labelsFlag = false;
   }

   /**
    * @name LaTex-specific methods
    * @{
    */

   /**
    * Not used anymore.
    */
   @Deprecated
   public void enableCustomLabels() {
      labelsFlag = true;
   }

   /**
    * Not used anymore.
    */
   @Deprecated
   public void disableCustomLabels() {
      labelsFlag = false;
   }

   /**
    * Sets the position of each label on this axis. This method requires
    * an array containing an increasing sequence of numbers representing
    * the positions at which labels will appear on the axis. It is
    * designed to export the axis to a LaTeX source code; it has no effect
    * on the chart appearance displayed with `XYChart.view()`.
    *  @param position     new label positions.
    */
   public void setLabels (double[] position) {
      this.labelsName = null;
      this.labelsValue = position.clone();
      double max = max(position);
      if (axis.getUpperBound() < max)
         axis.setUpperBound(max);
      double min = min(position);
      if (axis.getLowerBound() > min)
         axis.setLowerBound(min);
      labelsFlag = true;
   }

   /**
    * Assigns custom labels to user-defined positions on the axis. This
    * method requires an array of positions as well as an array of labels.
    * The label `label[i]` will be used at position `position[i]`. It is
    * designed to export the axis to a LaTeX source code, and to use
    * LaTeX/TikZ commands to write prettier characters; it has no effect
    * on the chart appearance displayed with `XYChart.view()`.
    *  @param position     label series position on the axis.
    *  @param label        label series name on the axis.
    */
   public void setLabels (double[] position, String[] label) {
      if (label.length != position.length)
         throw new IllegalArgumentException
            ("A label is required for each given position");
      this.labelsName = label.clone();
      this.labelsValue = position.clone();
      double max = max(position);
      if (axis.getUpperBound() < max)
         axis.setUpperBound(max);
      double min = min(position);
      if (axis.getLowerBound() > min)
         axis.setLowerBound(min);
      labelsFlag = true;
   }

   /**
    * Returns the drawing position parameter (default equals 0).
    *  @return drawing position parameter.
    */
   public double getTwinAxisPosition() {
      return this.twinAxisPosition;
   }

   /**
    * Defines where the opposite axis must be drawn on the current axis,
    * where it should appear, and on which label.
    *  @param position     new drawing position.
    */
   public void setTwinAxisPosition (double position) {
      this.twinAxisPosition = position;
   }

   /**
    * Formats and returns a string containing a LaTeX-compatible source
    * code which represents this axis and its parameters.
    *  @return LaTeX source code in a String.
    *
    *  @param scale        current axis wished scale.
    */
   public String toLatex (double scale) {
      Formatter formatter = new Formatter(Locale.US);
      double maxAxis = Math.max(axis.getRange().getUpperBound(), twinAxisPosition);  //valeur du label le plus grand
      double minAxis = Math.min(axis.getRange().getLowerBound(), twinAxisPosition);  //valeur du label le plus petit

      String precisionAffichageLabel;  //determine combien de decimales seront affichees pour les labels
      double pas = axis.getTickUnit().getSize();   //step d'affichage des labels
      pas = fixStep (minAxis, maxAxis, pas);

      int puissDix;                                //echelle des valeurs sur l'axe
      if (Math.log10(pas) < 0)
         puissDix = (int)Math.log10(pas) - 1;
      else
         puissDix = (int)Math.log10(pas);

      //Placement des fleches, on pourrait facilement les rendre personnalisables...
      String arrowLeftType = "latex";
      String arrowRightType = "latex";
      String arrowLeftMargin = "3mm";
      String arrowRightMargin = "3mm";
      if (twinAxisPosition == minAxis) {
         arrowLeftType = "";
         arrowLeftMargin = "0mm";
      } else if (twinAxisPosition == maxAxis) {
         arrowRightType = "";
         arrowRightMargin = "0mm";
      }

      // Label pour l'axe, avec eventuellement l'echelle
      String label = axis.getLabel();
      if (label == null)
         label = " ";

      if (puissDix < -2 || puissDix > 2)
         label = label + " $(10^{" + puissDix + "})$";
      else
         puissDix = 0;

      if (orientation) { //on est sur l'axe des abscisses

         //affichage de l'axe
         formatter.format("\\draw [%s-%s] ([xshift=-%s] %s,0) -- ([xshift=%s] %s,0) node[right] {%s};%n", arrowLeftType, arrowRightType, arrowLeftMargin, (minAxis - twinAxisPosition)*scale, arrowRightMargin, (maxAxis - twinAxisPosition)*scale, label);

         if (labelsFlag) {   //labels manuels
            String name;
            double value;
            double labelTemp;
            for (int i = 0; i < labelsValue.length; i++) {
               value = labelsValue[i];
               if (labelsName == null) {
                  labelTemp = (value * 100.0 / Math.pow(10, puissDix));
                  //                   if(labelTemp == (int)(labelTemp))
                  name = Integer.toString((int)(labelTemp / 100.0));
                  //                   else
                  //                      name = Double.toString(Math.round(labelTemp)/100.0);
               } else
                  name = labelsName[i];
               if (value == twinAxisPosition && tick0Flag)
                  formatter.format("\\draw (%s,00) -- +(0mm,1mm) -- +(0mm,-1mm) node[below right] {%s};%n", (value - twinAxisPosition)*scale, name);
               else
                  formatter.format("\\draw (%s,0) -- +(0mm,1mm) -- +(0mm,-1mm) node[below] {%s};%n", (value - twinAxisPosition)*scale, name);
            }
         } else {   //labels automatiques
            double labelTemp;
            double k = twinAxisPosition;
            labelTemp = (Math.round(k * 100 / Math.pow(10, puissDix))) / 100.0;
            if (labelTemp == (int)(labelTemp))
               precisionAffichageLabel = "0";
            else if (labelTemp*10 == (int)(labelTemp*10))
               precisionAffichageLabel = "1";
            else
               precisionAffichageLabel = "2";
            if (tick0Flag)
               formatter.format("\\draw (%s,0) -- +(0mm,1mm) -- +(0mm,-1mm) node[below right] {%." + precisionAffichageLabel + "f};%n", (k - twinAxisPosition)*scale, labelTemp);
            else
               formatter.format("\\draw (%s,0) -- +(0mm,1mm) -- +(0mm,-1mm) node[below] {%." + precisionAffichageLabel + "f};%n", (k - twinAxisPosition)*scale, labelTemp);
            k += pas;
            while (k <= maxAxis) { //cote positif
               labelTemp = (Math.round(k * 100 / Math.pow(10, puissDix))) / 100.0;
               if (labelTemp == (int)(labelTemp))
                  precisionAffichageLabel = "0";
               else if (labelTemp*10 == (int)(labelTemp*10))
                  precisionAffichageLabel = "1";
               else
                  precisionAffichageLabel = "2";
               formatter.format("\\draw (%s,0) -- +(0mm,1mm) -- +(0mm,-1mm) node[below] {%." + precisionAffichageLabel + "f};%n", (k - twinAxisPosition)*scale, labelTemp);
               k += pas;
            }

            k = twinAxisPosition - pas;
            while (k >= minAxis) { //cote negatif
               labelTemp = (Math.round(k * 100 / Math.pow(10, puissDix))) / 100.0;
               if (labelTemp == (int)(labelTemp))
                  precisionAffichageLabel = "0";
               else if (labelTemp*10 == (int)(labelTemp*10))
                  precisionAffichageLabel = "1";
               else
                  precisionAffichageLabel = "2";
               formatter.format("\\draw (%s,0) -- +(0mm,1mm) -- +(0mm,-1mm) node[below] {%." + precisionAffichageLabel + "f};%n", (k - twinAxisPosition)*scale, labelTemp);
               k -= pas;
            }
         }

      } else { //On est sur l'axe des ordonnees

         //affichage de l'axe
         formatter.format("\\draw [%s-%s] ([yshift=-%s] 0,%s) -- ([yshift=%s] 0, %s) node[above] {%s};%n", arrowLeftType, arrowRightType, arrowLeftMargin, (minAxis - twinAxisPosition)*scale, arrowRightMargin, (maxAxis - twinAxisPosition)*scale, label);

         if (labelsFlag) {
            //labels manuels
            String name;
            double value;
            double labelTemp;
            for (int i = 0; i < labelsValue.length; i++) {
               value = labelsValue[i];
               if (labelsName == null) {
                  labelTemp = (value * scale * 100 / Math.pow(10, puissDix));
                  if (labelTemp == (int)(labelTemp))
                     name = Integer.toString((int)(labelTemp / 100.0));
                  else
                     name = Double.toString(Math.round(labelTemp) / 100.0);
               } else
                  name = labelsName[i];
               if (value == twinAxisPosition && tick0Flag)
                  formatter.format("\\draw (%s,%s) -- +(1mm,0mm) -- +(-1mm,0mm) node[above left] {%s};%n", 0, (value - twinAxisPosition)*scale, name);
               else
                  formatter.format("\\draw (%s,%s) -- +(1mm,0mm) -- +(-1mm,0mm) node[left] {%s};%n", 0, (value - twinAxisPosition)*scale, name);
            }
         } else {
            //Les labels automatiques
            double k = twinAxisPosition;
            double labelTemp;
            labelTemp = (Math.round(k * 100 / Math.pow(10, puissDix))) / 100.0;
            if (labelTemp == (int)(labelTemp))
               precisionAffichageLabel = "0";
            else if (labelTemp*10 == (int)(labelTemp*10))
               precisionAffichageLabel = "1";
            else
               precisionAffichageLabel = "2";
            if (tick0Flag)
               formatter.format("\\draw (0,%s) -- +(1mm,0mm) -- +(-1mm,0mm) node[above left] {%." + precisionAffichageLabel + "f};%n", (k - twinAxisPosition)*scale, labelTemp);
            else
               formatter.format("\\draw (0,%s) -- +(1mm,0mm) -- +(-1mm,0mm) node[left] {%." + precisionAffichageLabel + "f};%n", (k - twinAxisPosition)*scale, labelTemp);
            k += pas;
            while (k <= maxAxis) { //cote positif de l'axe
               labelTemp = (Math.round(k * 100 / Math.pow(10, puissDix))) / 100.0;
               if (labelTemp == (int)(labelTemp))
                  precisionAffichageLabel = "0";
               else if (labelTemp*10 == (int)(labelTemp*10))
                  precisionAffichageLabel = "1";
               else
                  precisionAffichageLabel = "2";
               formatter.format("\\draw (0,%s) -- +(1mm,0mm) -- +(-1mm,0mm) node[left] {%." + precisionAffichageLabel + "f};%n", (k - twinAxisPosition)*scale, labelTemp);
               k += pas;
            }
            k = twinAxisPosition - pas;
            while (k >= minAxis) { //cote negatif de l'axe
               labelTemp = (Math.round(k * 100 / Math.pow(10, puissDix))) / 100.0;
               if (labelTemp == (int)(labelTemp))
                  precisionAffichageLabel = "0";
               else if (labelTemp*10 == (int)(labelTemp*10))
                  precisionAffichageLabel = "1";
               else
                  precisionAffichageLabel = "2";
               formatter.format("\\draw (0,%s) -- +(1mm,0mm) -- +(-1mm,0mm) node[left] {%." + precisionAffichageLabel + "f};%n", (k - twinAxisPosition)*scale, labelTemp);
               k -= pas;
            }
         }
      }
      return formatter.toString();
   }

   private double computeAutoTickValue() {
      // Calcule le pas d'affichage par defaut des labels
      double pas;
      double[] bounds = new double[2];
      bounds[0] = Math.min(axis.getRange().getLowerBound(), twinAxisPosition);
      bounds[1] = Math.max(axis.getRange().getUpperBound(), twinAxisPosition);
      if ( bounds[1] - bounds[0] >= 1) {
         int nbMaxDeFois = (int)(bounds[1] - bounds[0]);
         int puissDix = (int)Math.log10(nbMaxDeFois);
         nbMaxDeFois = (int)(nbMaxDeFois / Math.pow(10, puissDix)) + 1;
         if (nbMaxDeFois > 5)
            pas = Math.pow(10, puissDix);
         else if (nbMaxDeFois > 3)
            pas = 0.5 * Math.pow(10, puissDix);
         else if (nbMaxDeFois > 1)
            pas = 0.25 * Math.pow(10, puissDix);
         else //nbMaxDeFois==1
            pas = 0.1 * Math.pow(10, puissDix);
      } else {
         double nbMaxDeFois = bounds[1] - bounds[0];
         int puissDix = 0;
         while (nbMaxDeFois < 1) {
            nbMaxDeFois *= 10;
            puissDix++;
         }
         if (nbMaxDeFois > 5)
            pas = 1 / Math.pow(10, puissDix);
         else if (nbMaxDeFois > 3)
            pas = 0.5 / Math.pow(10, puissDix);
         else if (nbMaxDeFois > 1)
            pas = 0.3 / Math.pow(10, puissDix);
         else //nbMaxDeFois==1
            pas = 0.1 / Math.pow(10, puissDix);
      }
      return pas;
   }

   private static double max (double[] t) {
      if (t == null)
         throw new IllegalArgumentException ("max:   null argument.");
      double aux = t[0];
      for (int i = 1; i < t.length; i++)
         if (t[i] > aux)
            aux = t[i];
      return aux;
   }

   private static double min (double[] t) {
      if (t == null)
         throw new IllegalArgumentException ("min:   null argument.");
      double aux = t[0];
      for (int i = 1; i < t.length; i++)
         if (t[i] < aux)
            aux = t[i];
      return aux;
   }


   void setTick0Flag(boolean flag) {
        tick0Flag = flag;
   }
}

/**
 * @}
 */