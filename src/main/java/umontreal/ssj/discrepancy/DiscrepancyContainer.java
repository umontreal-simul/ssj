/*
 * Class:        DiscrepancyContainer
 * Description:
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Adam L'Archevêque Gaudet
 * @since        January 2009

 * SSJ is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, either version 3 of the License, or
 * any later version.

 * SSJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * A copy of the GNU General Public License is available at
   <a href="http://www.gnu.org/licenses">GPL licence site</a>.
 */
package umontreal.ssj.discrepancy;
import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.util.Num;
import umontreal.ssj.charts.XYLineChart;
import umontreal.ssj.functionfit.LeastSquares;
import java.io.*;

/**
 * This class is used to compute, store and display discrepancies. The method
 * #add computes the discrepancy of a given array of points for the selected
 * discrepancies and associates them with a given index. One can assign a
 * parameter value to this index (using  #setParam ) so that a graph or a
 * data file is created showing the discrepancies as functions of the
 * parameter. One can also scale the discrepancies of an index with a given
 * scale factor or take the logarithm of the discrepancies.
 *
 * The discrepancies are computed, assuming that the theoretical distribution
 * of the points is over the unit hypercube @f$[0,1]^s@f$; thus all the
 * coordinates of the points must be in @f$[0,1]@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class DiscrepancyContainer {
   private double[][] disc;
   private Discrepancy[] discrepancies;
   private double[] regression;
   private int n,nbDisc;
   private String title,xLabel, yLabel;

   /**
    * Creates a  @ref DiscrepancyContainer for the given discrepancies.
    *  @param discrepancies array containing the discrepancies to be used
    */
   public DiscrepancyContainer (Discrepancy[] discrepancies) {
      this.discrepancies = discrepancies;
      nbDisc = discrepancies.length;
      regression = new double[nbDisc + 1];
   }

   /**
    * Initialize the container with enough space for `n` values of the
    * parameter and sets the values to 0. Also sets the labels for the
    * parameter and the discrepancy that will be used for graphs and data
    * file.
    *  @param n            number of values the parameter will take
    *  @param title        title
    *  @param xLabel       label name for the parameter
    *  @param yLabel       label name for the discrepancy
    */
   public void init (int n, String title, String xLabel, String yLabel) {
      this.xLabel = xLabel;
      this.yLabel = yLabel;
      this.n = n;
      disc = new double[nbDisc+1][n];
      reset();
   }

   /**
    * Calls  {@link #init(int,String,String,String) init(n,"","Parameter",
    * "Discrepancy")}.
    *  @param n            Number of values the parameter will take.
    */
   public void init (int n) {
      init(n,"", "Parameter", "Discrepancy");
   }

   /**
    * Resets the values of the discrepancies at index `i` to 0.
    *  @param i            index to reset
    */
   public void reset (int i) {
      for (int j=0;j<nbDisc;++j)
         disc[j+1][i] = 0.0;
   }

   /**
    * Calls  {@link #reset(int) reset(i)} for all indices `i`.
    */
   public void reset() {
      for (int i=0; i<n; ++i)
         reset(i);
   }

   /**
    * Computes the discrepancies of the first `n` values contained in
    * `points` and sets the values at index `i`.
    *  @param i            index where to add the discrepancies
    *  @param points       values for which to compute the discrepancies
    *  @param n            number of points to use
    */
   public void compute (int i, double[] points, int n) {
      for (int j = 0; j < nbDisc; ++j)
         disc[j+1][i] = discrepancies[j].compute(points, n);
   }

   /**
    * Computes the discrepancies of the first `n` values contained in
    * `points` using the first `s` coordinates and sets the values at
    * index `i`.
    *  @param i            index where to add the discrepancies
    *  @param points       values for which to compute the discrepancies
    *  @param n            number of points to use
    *  @param s            number of coordinates to use for each point
    */
   public void compute (int i, double[][] points, int n, int s) {
      for (int j = 0; j < nbDisc; ++j)
         disc[j+1][i] = discrepancies[j].compute(points, n, s);
   }

   /**
    * Computes the discrepancies of the first `n` values contained in
    * `points`, and adds the values at index `i`. **NOTE: This method does
    * not replace the values, it adds to them. Can be used with  #scale to
    * calculate an average discrepancy over more than one point set.
    *  @param i            index where to add the discrepancies
    *  @param points       values for which to compute the discrepancies
    *  @param n            number of points to use
    */
   public void add (int i, double[] points, int n) {
      for (int j = 0; j < nbDisc; ++j)
         disc[j+1][i] += discrepancies[j].compute(points, n);
   }

   /**
    * Computes the square of the discrepancies of the first `n` values
    * contained in `points`, and adds the values at index `i`. **NOTE:
    * This method does not replace the values, it adds to them. Can be
    * used with  #scale to calculate an average square discrepancy over
    * more than one point set.
    *  @param i            index where to add the discrepancies
    *  @param points       values for which to compute the discrepancies
    *  @param n            number of points to use
    */
   public void addSquare (int i, double[] points, int n) {
      double discrepancy;
      for (int j = 0; j < nbDisc; ++j){
         discrepancy = discrepancies[j].compute(points, n);
         disc[j+1][i] += discrepancy*discrepancy;
      }
   }

   /**
    * Computes the discrepancies of the first `n` values contained in
    * `points` using the first `s` coordinates, and adds the values at
    * index `i`. **NOTE: This method does not replace the values, it adds
    * to them. Can be used with  #scale to calculate an average
    * discrepancy over more than one point set.
    *  @param i            index where to add the discrepancies
    *  @param points       values for which to compute the discrepancies
    *  @param n            number of points to use
    *  @param s            number of coordinates to use for each point
    */
   public void add (int i, double[][] points, int n, int s) {
      for (int j = 0; j < nbDisc; ++j)
         disc[j+1][i] += discrepancies[j].compute(points, n, s);
   }

   /**
    * Computes the square discrepancies of the first `n` values contained
    * in `points` using the first `s` coordinates, and adds the values at
    * index `i`. **NOTE: This method does not replace the values, it adds
    * to them. Can be used with  #scale to calculate an average square
    * discrepancy over more than one point set.
    *  @param i            index where to add the discrepancies
    *  @param points       values for which to compute the discrepancies
    *  @param n            number of points to use
    *  @param s            number of coordinates to use for each point
    */
   public void addSquare (int i, double[][] points, int n, int s) {
      double discrepancy;
      for (int j = 0; j < nbDisc; ++j){
         discrepancy = discrepancies[j].compute(points, n, s);
         disc[j+1][i] += discrepancy*discrepancy;
      }
   }

   /**
    * Multiplies all the discrepancies at index `i` by `scale`;
    *  @param i            index where to scale
    *  @param scale        scale factor
    */
   public void scale (int i, double scale) {
      for (int j = 0; j < nbDisc; ++j)
         disc[j+1][i] *= scale;
   }

   /**
    * Calls  {@link #scale(int,double) scale(i,scale)} for all indices
    * `i`.
    *  @param scale        scale factor
    */
   public void scale (double scale) {
      for (int i = 0; i < n; ++i)
         scale(i, scale);
   }

   /**
    * Takes the logarithm in base 2 of the discrepancy values at index
    * `i`.
    *  @param i            index where to take the logarithm
    */
   public void log2 (int i) {
      for (int j = 0; j < nbDisc; ++j)
         disc[j+1][i] = Num.log2(disc[j+1][i]);
   }

   /**
    * Squares the discrepancy values at index `i`.
    *  @param i            index where to take the logarithm
    */
   public void square (int i) {
      for (int j = 0; j < nbDisc; ++j)
         disc[j+1][i] *= disc[j+1][i];
   }

   /**
    * Sets the parameter value at index `i` to `parmValue`.
    *  @param i            index where to set the parameter value
    *  @param paramValue   value to set the parameter
    */
   public void setParam (int i, double paramValue) {
      disc[0][i] = paramValue;
   }

   // Removes infinite values caused by loss of precision in
   // discrepancy computations.
   private void removeInfinite(){
      for (int i=0; i<n; ++i){
         for (int j = 0; j < nbDisc; ++j)
            if (disc[j+1][i] == Double.NEGATIVE_INFINITY)
               disc[j+1][i] = 0.0;
      }
   }

/**
 * Computes the linear regression slope for the discrepancies as function of
 * the parameter. Used by  #regressionToString.
 */
protected void calcRegressionSlope() {
      for (int j = 0; j < nbDisc; ++j)
         regression[j] = LeastSquares.calcCoefficients(disc[0],disc[j+1],1)[1];
   }

   /**
    * Formats and returns a `String` containing the linear regression
    * slopes for the discrepancies as function of the parameter.
    */
   public String regressionToString() {
      calcRegressionSlope();
      StringBuffer sb = new StringBuffer
      ("***************************************************************" +
           PrintfFormat.NEWLINE);
      sb.append("Linear regression slope"+PrintfFormat.NEWLINE);
      for (int j = 0; j < nbDisc; ++j)
         sb.append (PrintfFormat.s(15, discrepancies[j].toString()) +
                    PrintfFormat.g(15, 6, regression[j]) +
                    PrintfFormat.NEWLINE);
      return sb.toString();
   }

   /**
    * Creates a file named `filename.tex` containing LaTeX code that can
    * be compiled by `pdfLaTeX` to a graph of the discrepancies as
    * function of the parameter.
    *  @param filename     name of the LaTeX file to be created (without
    *                      the `.tex` extension)
    */
   public void toTexFile (String filename) {
      removeInfinite();
 //     XYLineChart chart = new XYLineChart(title, xLabel, yLabel); // ambigu
      XYLineChart chart = new XYLineChart();
      chart.setTitle(title);
      chart.getXAxis().setLabel(xLabel);
      chart.getYAxis().setLabel(yLabel);
      for (int j = 0; j < nbDisc; ++j)
         chart.add(disc[0], disc[j+1], discrepancies[j].toString(),  "sharp plot");

      try {
         Writer file = new FileWriter (filename+".tex");
         file.write(chart.toLatex(9,5));
         file.close();
      } catch (IOException e) {
         e.printStackTrace();
      }
   }

   /**
    * Creates a file named `filename.dat` and writes in it the given
    * `header`, a table showing the discrepancies for the different values
    * of the parameter and the linear regression slopes.
    *  @param filename     name of the output file (without any extension)
    *  @param header       header of the file (can be empty)
    */
   public void toDatFile (String filename, String header) {
      try {
         Writer file = new FileWriter (filename+".dat");
         file.write(header);
         file.write(toString());
         file.write(regressionToString());
         file.close();
      } catch (IOException e) {
         e.printStackTrace();
      }
   }

   /**
    * Returns a `String` containing a table showing the discrepancies for
    * the different values of the parameter.
    */
   public String toString() {
      StringBuffer sb = new StringBuffer(PrintfFormat.NEWLINE);
      sb.append (PrintfFormat.s(35,xLabel));
      for(int i=0; i < n ; ++i)
         sb.append (PrintfFormat.g(15,6,  disc[0][i]));
      sb.append(PrintfFormat.NEWLINE);

      for (int j = 0; j < nbDisc; ++j) {
         sb.append (PrintfFormat.s(15,yLabel));
         sb.append (PrintfFormat.s(20, discrepancies[j].toString()));
         for(int i=0; i < n; ++i)
            sb.append (PrintfFormat.g(15, 6, disc[j+1][i]));
         sb.append(PrintfFormat.NEWLINE);
      }
      return sb.toString();
   }

}