/*
 * Class:        PlotFormat
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

import   java.util.Formatter;
import   java.util.Locale;
import   java.util.StringTokenizer;
import   java.util.regex.Pattern;

import   java.io.Reader;
import   java.io.StringReader;
import   java.io.IOException;

import   org.jfree.data.xy.XYSeriesCollection;
import   org.jfree.data.category.CategoryDataset;
import   org.jfree.data.io.CSV;

/**
 * Provide tools to import and export data set tables to and from Gnuplot,
 * MATLAB and Mathematica compatible formats, or customized format.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class PlotFormat {
   private PlotFormat() {}

   /**
    * @name Data import functions
    * @{
    */

   /**
    * Parses `data` according to the standard GNUPlot format and stores
    * the extracted values in the returned table. All unrecognized
    * characters and empty lines will be cut. Deleting comments inside the
    * input `String` is recommended.
    *  @param data         String to parse.
    *  @return Table that represent data contained in `data`.
    *
    *  @exception IllegalArgumentException if each line of the input
    * `String` doesn’t share the same number of values
    */
   public static double[][] fromGNUPlot (String data) {
      int maxLineNumber = 16;
      int columnNumber = 0;
      int lineNumber = 0;
      double[][] retour = new double[1][1]; //initialisation temporaire

      Pattern p = Pattern.compile("\\n\\n*");
      String[] string = p.split(data); // A chaque String correspond une ligne

      for (int j = 0; j < string.length; j++) {
         if (string[j].charAt(0) != '#' && string[j].length() != 0) {
            //On ne traite pas les lignes de commentaires ni les lignes

            p = Pattern.compile("\\s\\s*");
            String[] result = p.split(string[j]);
            double[] line = new double[result.length];

            int k = 0; // Compte le nombre de valeurs sur la ligne courante
            for (int i = 0; i < result.length; i++) {
               try {
                  line[k] = Double.valueOf(result[i]);   // les chaines de caracteres qui ne sont pas des nombres sont ignorees
                  k++;
               } catch (NumberFormatException nfe) {}
            }

            if (k != 0) { // il y a des valeurs sur la ligne
               //La premiere ligne determine combien de colonnes aura notre tableau de retour
               if (lineNumber == 0) {
                  columnNumber = k;
                  retour = new double[columnNumber][maxLineNumber]; // Allocation du tableau de retour
               } else if ( columnNumber != k ) // Toutes les lignes de la String n'ont pas le meme nombre de valeurs
                  throw new IllegalArgumentException("Each line must have the same number of values");

               if (lineNumber >= maxLineNumber) { // Cas ou notre tableau initial n'est plus assez grand
                  maxLineNumber = maxLineNumber * 2;
                  double[][] temp = new double[columnNumber][maxLineNumber];
                  for (int i = 0; i < columnNumber; i++)
                     for (int l = 0; l < retour[0].length; l++)
                        temp[i][l] = retour[i][l];
                  retour = temp;
               }

               for (int i = 0; i < columnNumber; i++)
                  retour[i][lineNumber] = line[i];
               lineNumber++;
            }
         }
      }
      // Cree un tableau aux dimensions exactes
      double[][] temp;
      if (columnNumber == 0)
         temp = null;
      else
         temp = new double[columnNumber][lineNumber];
      for (int i = 0; i < columnNumber; i++)
         for (int l = 0; l < lineNumber; l++)
            temp[i][l] = retour[i][l];
      retour = null;
      return temp;
   }

   /**
    * Parses `data` according to the standard CSV format and stores the
    * extracted values in the returned table.
    *  @param data         String to parse.
    *  @return Table that represent data contained in `data`.
    *
    *  @exception IllegalArgumentException if each line of the input
    * `String` doesn’t share the same number of values
    */
   public static double[][] fromCSV (String data) {
      int maxLineNumber = 16;
      int columnNumber = 0;
      int lineNumber = 0;
      double[][] retour = new double[1][1]; //initialisation temporaire

      Pattern p = Pattern.compile("\\n\\n*");
      String[] string = p.split(data); // A chaque String correspond une ligne

      for (int j = 0; j < string.length; j++) {
         if (string[j].length() != 0) {
            //On ne traite pas les lignes vides

            p = Pattern.compile(",,*");
            String[] result = p.split(string[j]);
            double[] line = new double[result.length];

            int k = 0; // Compte le nombre de valeurs sur la ligne courante
            for (int i = 0; i < result.length; i++) {
               try {
                  line[k] = Double.valueOf(result[i]);   // les chaines de caracteres qui ne sont pas des nombres sont ignorees
                  k++;
               } catch (NumberFormatException nfe) {}
            }

            if (k != 0) { // il y a des valeurs sur la ligne
               //La premiere ligne determine combien de colonnes aura notre tableau de retour
               if (lineNumber == 0) {
                  columnNumber = k;
                  retour = new double[columnNumber][maxLineNumber]; // Allocation du tableau de retour
               } else if ( columnNumber != k ) // Toutes les lignes de la String n'ont pas le meme nombre de valeurs
                  throw new IllegalArgumentException("Each line must have the same number of values");

               if (lineNumber >= maxLineNumber) { // Cas ou notre tableau initial n'est plus assez grand
                  maxLineNumber = maxLineNumber * 2;
                  double[][] temp = new double[columnNumber][maxLineNumber];
                  for (int i = 0; i < columnNumber; i++)
                     for (int l = 0; l < retour[0].length; l++)
                        temp[i][l] = retour[i][l];
                  retour = temp;
               }

               for (int i = 0; i < columnNumber; i++)
                  retour[i][lineNumber] = line[i];
               lineNumber++;
            }
         }
      }
      // Cree un tableau aux dimensions exactes
      double[][] temp;
      if (columnNumber == 0)
         temp = null;
      else
         temp = new double[columnNumber][lineNumber];
      for (int i = 0; i < columnNumber; i++)
         for (int l = 0; l < lineNumber; l++)
            temp[i][l] = retour[i][l];
      retour = null;
      return temp;
   }

   /**
    * Parses `data` according to a user defined format and stores the
    * extracted values in the returned table. `betweenValues` sets
    * characters between values on the same line and `endLine` sets
    * characters which separates each line of the input data set. Usually,
    * this parameter contains the classic end of line character
    * <tt>%n</tt>. This method uses regular expressions, so it is possible
    * to use regular characters as defined in the standard java API,
    * specially in the class  Pattern (package java.util.regex).
    *  @param betweenValues `String` which separates values on the same
    *                       line.
    *  @param endLine      `String` which separates lines.
    *  @param data         String to parse.
    *  @return Table that represent data contained in `data`.
    *
    *  @exception IllegalArgumentException if each line of the input
    * `String` doesn’t share the same number of values
    */
   public static double[][] fromCustomizedFormat (String betweenValues,
                                           String endLine, String data) {
      int maxLineNumber = 16;
      int columnNumber = 0;
      int lineNumber = 0;
      double[][] retour = new double[1][1]; //initialisation temporaire

      Pattern p = Pattern.compile("(" + endLine + ")(" + endLine + ")*");
      String[] string = p.split(data); // A chaque String correspond une ligne

      for (int j = 0; j < string.length; j++) {
         if (string[j].length() != 0) {
            //On ne traite pas les lignes vides

            p = Pattern.compile("(" + betweenValues + ")(" + betweenValues + ")*");
            String[] result = p.split(string[j]);
            double[] line = new double[result.length];

            int k = 0; // Compte le nombre de valeurs sur la ligne courante
            for (int i = 0; i < result.length; i++) {
               try {
                  line[k] = Double.valueOf(result[i]);   // les chaines de caracteres qui ne sont pas des nombres sont ignorees
                  k++;
               } catch (NumberFormatException nfe) {}
            }

            if (k != 0) { // il y a des valeurs sur la ligne
               //La premiere ligne determine combien de colonnes aura notre tableau de retour
               if (lineNumber == 0) {
                  columnNumber = k;
                  retour = new double[columnNumber][maxLineNumber]; // Allocation du tableau de retour
               } else if ( columnNumber != k ) // Toutes les lignes de la String n'ont pas le meme nombre de valeurs
                  throw new IllegalArgumentException("Each line must have the same number of values");

               if (lineNumber >= maxLineNumber) { // Cas ou notre tableau initial n'est plus assez grand
                  maxLineNumber = maxLineNumber * 2;
                  double[][] temp = new double[columnNumber][maxLineNumber];
                  for (int i = 0; i < columnNumber; i++)
                     for (int l = 0; l < retour[0].length; l++)
                        temp[i][l] = retour[i][l];
                  retour = temp;
               }

               for (int i = 0; i < columnNumber; i++)
                  retour[i][lineNumber] = line[i];
               lineNumber++;
            }
         }
      }
      // Cree un tableau aux dimensions exactes
      double[][] temp;
      if (columnNumber == 0)
         temp = null;
      else
         temp = new double[columnNumber][lineNumber];
      for (int i = 0; i < columnNumber; i++)
         for (int l = 0; l < lineNumber; l++)
            temp[i][l] = retour[i][l];
      retour = null;
      return temp;
   }

   /**
    * @}
    */

   /**
    * @name Data export functions
    * @{
    */

   /**
    * Stores data tables `data` into a `String`, with format
    * understandable by GNUPlot.
    *  @param data         data tables.
    *  @return String that represent data tables in GNUPlot format.
    */
   public static String toGNUPlot (double[]... data) {
      checkData(data);
      Formatter formatter = new Formatter(Locale.US);

      for (int i = 0; i < data[0].length; i++) {
         for (int j = 0; j < data.length; j++) {
            formatter.format("%20f", data[j][i]);
         }
         formatter.format("%n");
      }
      formatter.format("%n%n");
      return formatter.toString();
   }

   /**
    * Stores series collection `data` into a `String`, with format
    * understandable by GNUPlot.
    *  @param data         data tables.
    *  @return String that represent data tables in GNUPlot format.
    */
   public static String toGNUPlot (XYSeriesCollection data) {
      return toGNUPlot(toTable(data));
   }

   /**
    * Stores data tables `data` into a `String` with format CSV
    * (comma-separated value tabular data). The output string could be
    * imported from a file to a matrix into Mathematica with Mathematica’s
    * function `Import["fileName", "CSV"]`, or into MATLAB with MATLAB’s
    * function `csvread(’fileName’)`.
    *  @param data         data tables.
    *  @return String that represent data tables in CSV format.
    */
   public static String toCSV (double[]...data) {
      checkData(data);
      Formatter formatter = new Formatter(Locale.US);

      for (int i = 0; i < data[0].length; i++) {
         for (int j = 0; j < data.length - 1; j++) {
            formatter.format("%-20f, ", data[j][i]);
         }
         formatter.format("%-20f%n", data[data.length-1][i]);      //le dernier
      }
      formatter.format("%n");

      return formatter.toString();
   }

   /**
    * Stores series collection `data` into a `String` with format CSV
    * (comma-separated value tabular data).
    *  @param data         data tables.
    *  @return String that represent data tables in CSV format.
    */
   public static String toCSV (XYSeriesCollection data) {
      return toCSV(toTable(data));
   }

   /**
    * Stores data tables `data` into a `String` with customized format.
    * `heading` sets the head of the returned `String`, `footer` sets its
    * end, `betweenValues` sets characters between values on the same line
    * and finally `endLine` sets characters which separates each line of
    * the input data set. Normally, this parameter contains the classic
    * end of line character ‘%n’.
    *  @param heading      head of the returned `String`.
    *  @param footer       end of the returned `String`.
    *  @param betweenValues `String` which separates values on the same
    *                       line.
    *  @param endLine      `String` which separates lines.
    *  @param precision    number of decimal
    *  @param data         data tables.
    *  @return String that represent data tables in customized format.
    */
   public static String toCustomizedFormat (String heading, String footer,
                                   String betweenValues, String endLine,
                                   int precision, double[]...data) {
      checkData(data);
      Formatter formatter = new Formatter(Locale.US);
      String myString = "%20."+ precision +"f";

      formatter.format("%s", heading);
      for (int i = 0; i < data[0].length; i++) {
         for (int j = 0; j < data.length - 1; j++) {
            formatter.format(myString+"%s", data[j][i], betweenValues);
         }
         formatter.format(myString+"%s", data[data.length-1][i], endLine);      //le dernier de la ligne
      }
      formatter.format("%s", footer);
      return formatter.toString();
   }

   /**
    * Stores data tables `data` into a `String` with customized format
    * from an `XYSeriesCollection` variable.
    *  @param heading      head of the returned `String`.
    *  @param footer       end of the returned `String`.
    *  @param betweenValues `String` which separates values on the same
    *                       line.
    *  @param endLine      `String` which separates lines.
    *  @param precision    number of decimal
    *  @param data         data tables.
    *
    *  @return String that represent data tables in customized format.
    */
   public static String toCustomizedFormat (String heading, String footer,
                                   String betweenValues, String endLine,
                                   int precision, XYSeriesCollection data) {
      return toCustomizedFormat (heading, footer, betweenValues, endLine,
                       precision, toTable(data));
   }

   /**
    * Converts a <TT>XYSeriesCollection</TT> object into a <TT>double[][]</TT>.
    */
   private static double[][] toTable (XYSeriesCollection data) {
      double[][] transform = new double[data.getSeriesCount()*2][];

      for(int i = 0; i<data.getSeriesCount(); i++) {
         double[][] temp = data.getSeries(i).toArray();
         transform[2*i] = temp[0];
         transform[2*i+1] = temp[1];
      }
      return transform;
   }


   /**
    * Check the data table <TT>data</TT>.
    *
    * @exception  IllegalArgumentException   If the input tables doesn't have the same length.
    */
   private static void checkData(double[]...data) {
      for(int i = 0; i < data.length-1; i++) {
         if(data[i].length != data[i+1].length)
            throw new IllegalArgumentException("Data tables " + i + " and " + (i+1) + " must share the same length");
      }
   }

/* //Inutile si on exporte en CSV
   private static String mathematicaFormatPoint (double x, double y) {
      // Writes the pair (x, y) in returned string, in a format understood
      // by Mathematica
      StringBuffer sb = new StringBuffer();
      String S;

      sb.append ("   { ");
      if ((x != 0.0) && (x < 0.1 || x > 1.0)) {
         S = PrintfFormat.E (16, 7, x);
         int exppos = S.indexOf ('E');
         if (exppos != -1)
            S = S.substring (0, exppos) + "*10^(" +
                             S.substring (exppos+1) + ")";
      }
      else
         S = PrintfFormat.g (16, 8, x);

      sb.append (S + ",     ");

      if (y != 0.0 && (y < 0.1 || y > 1.0)) {
         S = PrintfFormat.E (16, 7, y);
         int exppos = S.indexOf ('E');
         if (exppos != -1)
            S = S.substring (0, exppos) + "*10^(" +
                             S.substring (exppos+1) + ")";
      }
      else
        S = PrintfFormat.g (16, 8, y);

      sb.append (S + " }");
      return sb.toString();
   }*/

}

/**
 * @}
 */