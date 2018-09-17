/*
 * Class:        TextDataReader
 * Description:  Provides static methods to read data from text files
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
package umontreal.ssj.util;

import java.io.LineNumberReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.net.URL;
import java.util.logging.Logger;

/**
 * Provides static methods to read data from text files.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class TextDataReader {
   private static Logger log = Logger.getLogger ("umontreal.ssj.util");

/**
 * Reads an array of double-precision values from the reader `input`. For
 * each line of text obtained from the given reader, this method trims
 * whitespaces, and parses the remaining text as a double-precision value.
 * This method ignores every character other than the digits, the plus and
 * minus signs, the period (<tt>.</tt>), and the letters `e` and `E`.
 * Moreover, lines starting with a pound sign (<tt>#</tt>) are considered as
 * comments and thus skipped. The method returns an array containing all the
 * parsed values.
 *  @param input        the reader to obtain data from.
 *  @return the obtained array of double-precision values.
 *
 *  @exception IOException if an I/O error occurs.
 */
public static double[] readDoubleData (Reader input) throws IOException {
      LineNumberReader inb = new LineNumberReader (input);
      double[] data = new double[5];
      int n = 0;
      String li;
      while ((li = inb.readLine()) != null) {
        li = li.trim();
        if (li.startsWith ("#"))
           continue;

         // look for the first non-digit character on the read line
         int index = 0;
         while (index < li.length() &&
            (li.charAt (index) == '+' || li.charAt (index) == '-' ||
             li.charAt (index) == 'e' || li.charAt (index) == 'E' ||
             li.charAt (index) == '.' || Character.isDigit (li.charAt (index))))
           ++index; 

         // truncate the line
         li = li.substring (0, index);
         if (!li.equals ("")) {
            try {
               data[n++] = Double.parseDouble (li);
               if (n >= data.length) {
                  double[] newData = new double[2*n];
                  System.arraycopy (data, 0, newData, 0, data.length);
                  data = newData;
               }
            }
            catch (NumberFormatException nfe) {
               log.warning ("Invalid line " + inb.getLineNumber() + ": " + li);
            }
         }
      }
      if (data.length != n) {
         double[] data2 = new double[n];
         System.arraycopy (data, 0, data2, 0, n);
         return data2;
      }
      return data;
   }

   /**
    * Connects to the URL referred to by the URL object `url`, and calls
    * #readDoubleData(Reader) to obtain an array of double-precision
    * values from the resource.
    *  @param url          the URL object representing the resource to
    *                      read.
    *  @return the obtained array of double-precision values.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static double[] readDoubleData (URL url) throws IOException {
      Reader reader = new InputStreamReader (url.openStream());
      try {
         return readDoubleData (reader);
      }
      finally {
         reader.close();
      }
   }

   /**
    * Opens the file referred to by the file object `file`, and calls
    * #readDoubleData(Reader) to obtain an array of double-precision
    * values from the file.
    *  @param file         the file object representing the file to read.
    *  @return the obtained array of double-precision values.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static double[] readDoubleData (File file) throws IOException {
      FileReader reader = new FileReader (file);
      try {
         return readDoubleData (reader);
      }
      finally {
         reader.close();
      }
   }

   /**
    * Opens the file with name `file`, and calls  #readDoubleData(Reader)
    * to obtain an array of double-precision values from the file.
    *  @param file         the name of the file to read.
    *  @return the obtained array of double-precision values.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static double[] readDoubleData (String file) throws IOException {
      FileReader reader = new FileReader (file);
      try {
         return readDoubleData (reader);
      }
      finally {
         reader.close();
      }
   }

   /**
    * This is equivalent to  #readDoubleData(Reader), for reading
    * integers.
    *  @param input        the reader to obtain data from.
    *  @return the obtained array of integers.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static int[] readIntData (Reader input) throws IOException {
      LineNumberReader inb = new LineNumberReader (input);
      int[] data = new int[5];
      int n = 0;
      String li;
      while ((li = inb.readLine()) != null) {
        li = li.trim();
        if (li.startsWith ("#"))
           continue;

         // look for the first non-digit character on the read line
         int index = 0;
         while (index < li.length() &&
            (li.charAt (index) == '+' || li.charAt (index) == '-' ||
             Character.isDigit (li.charAt (index))))
           ++index; 

         // truncate the line
         li = li.substring (0, index);
         if (!li.equals ("")) {
            try {
               data[n++] = Integer.parseInt (li);
               if (n >= data.length) {
                  int[] newData = new int[2*n];
                  System.arraycopy (data, 0, newData, 0, data.length);
                  data = newData;
               }
            }
            catch (NumberFormatException nfe) {
               log.warning ("Invalid line " + inb.getLineNumber() + ": " + li);
            }
         }
      }
      if (data.length != n) {
         int[] data2 = new int[n];
         System.arraycopy (data, 0, data2, 0, n);
         return data2;
      }
      return data;
   }

   /**
    * Connects to the URL referred to by the URL object `url`, and calls
    * #readIntData(Reader) to obtain an array of integers from the
    * resource.
    *  @param url          the URL object representing the resource to
    *                      read.
    *  @return the obtained array of integers.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static int[] readIntData (URL url) throws IOException {
      Reader reader = new InputStreamReader (url.openStream());
      try {
         return readIntData (reader);
      }
      finally {
         reader.close();
      }
   }

   /**
    * This is equivalent to  #readDoubleData(File), for reading integers.
    *  @param file         the file object represented to file to read.
    *  @return the array of integers.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static int[] readIntData (File file) throws IOException {
      FileReader reader = new FileReader (file);
      try {
         return readIntData (reader);
      }
      finally {
         reader.close();
      }
   }

   /**
    * This is equivalent to  #readDoubleData(String), for reading
    * integers.
    *  @param file         the name of the file to read.
    *  @return the array of integers.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static int[] readIntData (String file) throws IOException {
      FileReader reader = new FileReader (file);
      try {
         return readIntData (reader);
      }
      finally {
         reader.close();
      }
   }

   /**
    * Reads an array of strings from the reader `input`. For each line of
    * text obtained from the given reader, this method trims leading and
    * trailing whitespaces, and stores the remaining string. Lines
    * starting with a pound sign (<tt>#</tt>) are considered as comments
    * and thus skipped. The method returns an array containing all the
    * read strings.
    *  @param input        the reader to obtain data from.
    *  @return the obtained array of strings.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static String[] readStringData (Reader input) throws IOException {
      LineNumberReader inb = new LineNumberReader (input);
      String[] data = new String[5];
      int n = 0;
      String li;
      while ((li = inb.readLine()) != null) {
        li = li.trim();
        if (li.startsWith ("#"))
           continue;

        data[n++] = li;
        if (n >= data.length) {
           String[] newData = new String[2*n];
           System.arraycopy (data, 0, newData, 0, data.length);
           data = newData;
        }
      }
      if (data.length != n) {
         String[] data2 = new String[n];
         System.arraycopy (data, 0, data2, 0, n);
         return data2;
      }
      return data;
   }

   /**
    * Connects to the URL referred to by the URL object `url`, and calls
    * #readStringData(Reader) to obtain an array of integers from the
    * resource.
    *  @param url          the URL object representing the resource to
    *                      read.
    *  @return the obtained array of strings.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static String[] readStringData (URL url) throws IOException {
      Reader reader = new InputStreamReader (url.openStream());
      try {
         return readStringData (reader);
      }
      finally {
         reader.close();
      }
   }

   /**
    * This is equivalent to  #readDoubleData(File), for reading strings.
    *  @param file         the file object represented to file to read.
    *  @return the array of strings.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static String[] readStringData (File file) throws IOException {
      FileReader reader = new FileReader (file);
      try {
         return readStringData (reader);
      }
      finally {
         reader.close();
      }
   }

   /**
    * This is equivalent to  #readDoubleData(String), for reading strings.
    *  @param file         the name of the file to read.
    *  @return the array of strings.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static String[] readStringData (String file) throws IOException {
      FileReader reader = new FileReader (file);
      try {
         return readStringData (reader);
      }
      finally {
         reader.close();
      }
   }

   /**
    * Uses the reader `input` to obtain a 2-dimensional array of
    * double-precision values. For each line of text obtained from the
    * given reader, this method trims whitespaces, and parses the
    * remaining text as an array of double-precision values. Every
    * character other than the digits, the plus (<tt>+</tt>) and minus
    * (<tt>-</tt>) signs, the period (<tt>.</tt>), and the letters `e` and
    * `E` are ignored and can be used to separate numbers on a line.
    * Moreover, lines starting with a pound sign (<tt>#</tt>) are
    * considered as comments and thus skipped. The lines containing only a
    * semicolon sign (<tt>;</tt>) are considered as empty lines. The
    * method returns a 2D array containing all the parsed values. The
    * returned array is not always rectangular.
    *  @param input        the reader to obtain data from.
    *  @return the 2D array of double-precison values.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static double[][] readDoubleData2D (Reader input)
                                              throws IOException {
      LineNumberReader inb = new LineNumberReader (input);
      double[][] data = new double[5][];
      int n = 0;
      String li;
      String number;

      while ((li = inb.readLine()) != null) {
         li = li.trim();
         if (li.startsWith ("#"))
            continue;

         if (li.equals(";")) {
            data[n++] = new double[0];
         }
         else {

            int index = 0;
            int begin = 0;
            boolean end = false;

            double[] row = new double[5];
            int k = 0;

            while (index < li.length() && (! end))
            {
               while (index < li.length() &&
                  (li.charAt (index) == '+' || li.charAt (index) == '-' ||
                   li.charAt (index) == 'e' || li.charAt (index) == 'E' ||
                   li.charAt (index) == '.' || Character.isDigit (li.charAt (index))))
                  ++index;

               if (index >= li.length() || (Character.isWhitespace (li.charAt (index))))
               {
                  number = li.substring (begin, index);
                  begin = ++index;

                  if (! number.equals("")) {
                     try {
                        row[k++] = Double.parseDouble (number);
                        if (k >= row.length) {
                           double[] newRow = new double[2*k];
                           System.arraycopy (row, 0, newRow, 0, row.length);
                           row = newRow;
                        }
                     }
                     catch (NumberFormatException nfe) {
                        log.warning ("Invalid column " + k + " at line " + inb.getLineNumber() + ": " + number);
                     }
                  }
               }
               else {
                  end = true;
               }
            }

            if (k > 0) {
               data[n] = new double[k];
               System.arraycopy (row, 0, data[n], 0, k);
               n++;
            }
            else {
               log.warning ("Invalid line " + inb.getLineNumber() + ": " + li);
            }
         }

         if (n == data.length) {
            double[][] newData = new double[2*n][];
            System.arraycopy (data, 0, newData, 0, n);
            data = newData;
         }
      }

      double[][] data2 = new double[n][];
      System.arraycopy (data, 0, data2, 0, n);
      return data2;
   }

   /**
    * Connects to the URL referred to by the URL object `url`, and calls
    * #readDoubleData2D(Reader) to obtain a matrix of double-precision
    * values from the resource.
    *  @param url          the URL object representing the resource to
    *                      read.
    *  @return the obtained matrix of double-precision values.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static double[][] readDoubleData2D (URL url) throws IOException {
      Reader reader = new InputStreamReader (url.openStream());
      try {
         return readDoubleData2D (reader);
      }
      finally {
         reader.close();
      }
   }

   /**
    * Opens the file referred to by the file object `file`, and calls
    * #readDoubleData2D(Reader) to obtain a matrix of double-precision
    * values from the file.
    *  @param file         the file object representing the file to read.
    *  @return the obtained matrix of double-precision values.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static double[][] readDoubleData2D (File file) throws IOException {
      FileReader reader = new FileReader (file);
      try {
         return readDoubleData2D (reader);
      }
      finally {
         reader.close();
      }
   }

   /**
    * Opens the file with name `file`, and calls
    * #readDoubleData2D(Reader) to obtain a matrix of double-precision
    * values from the file.
    *  @param file         the name of the file to read.
    *  @return the obtained matrix of double-precision values.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static double[][] readDoubleData2D (String file)
                                              throws IOException {
      FileReader reader = new FileReader (file);
      try {
         return readDoubleData2D (reader);
      }
      finally {
         reader.close();
      }
   }

   /**
    * This is equivalent to  #readDoubleData2D(Reader), for reading
    * integers.
    *  @param input        the reader to obtain data from.
    *  @return the obtained 2D array of integers.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static int[][] readIntData2D (Reader input) throws IOException {
      LineNumberReader inb = new LineNumberReader (input);
      int[][] data = new int[5][];
      int n = 0;
      String li;
      String number;

      while ((li = inb.readLine()) != null) {
         li = li.trim();
         if (li.startsWith ("#"))
            continue;

         if (li.equals(";")) {
            data[n++] = new int[0];
         }
         else {

            int index = 0;
            int begin = 0;
            boolean end = false;

            int[] row = new int[5];
            int k = 0;

            while (index < li.length() && (! end))
            {
               while (index < li.length() &&
                  (li.charAt (index) == '+' || li.charAt (index) == '-' ||
                   Character.isDigit (li.charAt (index))))
                  ++index;

               if (index >= li.length() || (Character.isWhitespace (li.charAt (index))))
               {
                  number = li.substring (begin, index);
                  begin = ++index;

                  if (! number.equals("")) {
                     try {
                        row[k++] = Integer.parseInt (number);
                        if (k >= row.length) {
                           int[] newRow = new int[2*k];
                           System.arraycopy (row, 0, newRow, 0, row.length);
                           row = newRow;
                        }
                     }
                     catch (NumberFormatException nfe) {
                        log.warning ("Invalid column " + k + " at line " + inb.getLineNumber() + ": " + number);
                     }
                  }
               }
               else {
                  end = true;
               }
            }

            if (k > 0) {
               data[n] = new int[k];
               System.arraycopy (row, 0, data[n], 0, k);
               n++;
            }
            else {
               log.warning ("Invalid line " + inb.getLineNumber() + ": " + li);
            }
         }

         if (n == data.length) {
            int[][] newData = new int[2*n][];
            System.arraycopy (data, 0, newData, 0, n);
            data = newData;
         }
      }

      int[][] data2 = new int[n][];
      System.arraycopy (data, 0, data2, 0, n);
      return data2;
   }

   /**
    * Connects to the URL referred to by the URL object `url`, and calls
    * #readDoubleData(Reader) to obtain a matrix of integers from the
    * resource.
    *  @param url          the URL object representing the resource to
    *                      read.
    *  @return the obtained matrix of integers.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static int[][] readIntData2D (URL url) throws IOException {
      Reader reader = new InputStreamReader (url.openStream());
      try {
         return readIntData2D (reader);
      }
      finally {
         reader.close();
      }
   }

   /**
    * This is equivalent to  #readDoubleData2D(File), for reading
    * integers.
    *  @param file         the file object represented to file to read.
    *  @return the obtained matrix of integer values.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static int[][] readIntData2D (File file) throws IOException {
      FileReader reader = new FileReader (file);
      try {
         return readIntData2D (reader);
      }
      finally {
         reader.close();
      }
   }

   /**
    * This is equivalent to  #readDoubleData2D(String), for reading
    * integers.
    *  @param file         the name of the file to read.
    *  @return the obtained matrix of integer values.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static int[][] readIntData2D (String file) throws IOException {
      FileReader reader = new FileReader (file);
      try {
         return readIntData2D (reader);
      }
      finally {
         reader.close();
      }
   }

   /**
    * Reads comma-separated values (CSV) from reader `input`, and returns
    * a 2D array of strings corresponding to the read data. Lines are
    * delimited using line separators <tt>\\r</tt>, <tt>\\n</tt>, and
    * <tt>\\r\\n</tt>. Each line contains one or more values, separated by
    * the column delimiter `colDelim`. If a string of characters is
    * surrounded with the string delimiter `stringDelim`, any line
    * separator and column separator appear in the string. The string
    * delimiter can be inserted in such a string by putting it twice.
    * Usually, the column delimiter is the comma, and the string delimiter
    * is the quotation mark. The following example uses these default
    * delimiters. <tt>
    * <pre>
    *          "One","Two","Three"
    *           1,2,3
    *          "String with "" delimiter",n,m
    * </pre>
    * </tt> This produces a matrix of strings with dimensions
    * @f$3\times3@f$. The first row contains the strings `One`, `Two`,
    * and `Three` while the second row contains the strings `1`, `2`, and
    * `3`. The first column of the last row contains the string `String
    * with " delimiter`.
    *  @param input        the reader to obtain data from.
    *  @param colDelim     the column delimiter.
    *  @param stringDelim  the string delimiter.
    *  @return the obtained 2D array of strings.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static String[][] readCSVData (Reader input, char colDelim,
                                         char stringDelim)
                                         throws IOException {
      // Using a buffered reader is important here for performance
      // LineNumberReader is a subclass of BufferedReader
      LineNumberReader inb = new LineNumberReader (input);
      StringBuffer sb = new StringBuffer();
      boolean stringMode = false;
      String[][] data = new String[5][];
      int numRows = 0;
      int numColumns = 0;
      boolean newRow = false, newColumn = false;
      int ich = -1;
      char ch = ' ';
      boolean readDone = false;
      while (!readDone) {
         if (ich == -2)
            // A character is pending
            ich = 0;
         else {
            ich = inb.read();
            if (ich == -1)
               // End of stream: process the last column and row, and exit
               newRow = newColumn = readDone = true;
            else
               ch = (char)ich;
         }
         if (ich != -1) {
            if (stringMode) {
               if (ch == stringDelim) {
                  // Check if there is a second string delimiter
                  int ichNext = inb.read();
                  if (ichNext >= 0) {
                     char chNext = (char)ichNext;
                     if (chNext == stringDelim)
                        // Append the quoted string delimiter
                        sb.append (stringDelim);
                     else {
                        // Indicate the end of the string, and a new pending character
                        stringMode = false;
                        ich = -2;
                        ch = chNext;
                     }
                  }
               }
               else
                  sb.append (ch);
            }
            else {
              if (ch == '\n' || ch == '\r') {
                 int ichNext = inb.read();
                 if (ichNext >= 0) {
                    char chNext = (char)ichNext;
                    if (ch == '\r' && chNext == '\n') {
                       ichNext = inb.read();
                       if (ichNext >= 0) {
                          chNext = (char)ichNext;
                          ich = -2;
                          ch = chNext;
                          newRow = true;
                       }
                    }
                    else {
                       ich = -2;
                       ch = chNext;
                       newRow = true;
                    }
                 }
              }
              else if (ch == colDelim)
                 newColumn = true;
              else if (ch == stringDelim)
                 stringMode = true;
              else
                 sb.append (ch);
            }
         }
         if (newColumn || newRow) {
            if (numColumns == 0) {
               ++numRows;
               numColumns = 1;
            }
            else
               ++numColumns;
            if (data.length < numRows) {
               String[][] newData = new String[2*data.length][];
               System.arraycopy (data, 0, newData, 0, data.length);
               data = newData;
            }
            if (data[numRows - 1] == null)
               data[numRows - 1] = new String[5];
            else if (data[numRows - 1].length < numColumns) {
               String[] newData = new String[2*data[numRows - 1].length];
               System.arraycopy (data[numRows - 1], 0, newData, 0, data[numRows - 1].length);
               data[numRows - 1] = newData;
            }
            data[numRows - 1][numColumns - 1] = sb.toString();
            sb.delete (0, sb.length());
            newColumn = false;
         }
         if (newRow) {
            if (data[numRows - 1].length != numColumns) {
               String[] data2 = new String[numColumns];
               System.arraycopy (data[numRows - 1], 0, data2, 0, numColumns);
               data[numRows - 1] = data2;
            }
            numColumns = 0;
            newRow = false;
         }
      }

      if (stringMode)
         throw new IllegalArgumentException ("Too many string delimiters " + stringDelim);
      if (data.length != numRows) {
         String[][] data2 = new String[numRows][];
         System.arraycopy (data, 0, data2, 0, numRows);
         return data2;
      }
      return data;
   }

   /**
    * Connects to the URL referred to by the URL object `url`, and calls
    * #readCSVData(Reader,char,char) to obtain a matrix of strings from
    * the resource.
    *  @param url          the URL object representing the resource to
    *                      read.
    *  @param colDelim     the column delimiter.
    *  @param stringDelim  the string delimiter.
    *  @return the obtained matrix of strings.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static String[][] readCSVData (URL url, char colDelim,
                                         char stringDelim)
                                         throws IOException {
      Reader reader = new InputStreamReader (url.openStream());
      try {
         return readCSVData (reader, colDelim, stringDelim);
      }
      finally {
         reader.close();
      }
   }

   /**
    * This is equivalent to  #readDoubleData2D(File), for reading strings.
    *  @param file         the file object represented to file to read.
    *  @param colDelim     the column delimiter.
    *  @param stringDelim  the string delimiter.
    *  @return the obtained matrix of string values.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static String[][] readCSVData (File file, char colDelim,
                                         char stringDelim)
                                         throws IOException {
      FileReader reader = new FileReader (file);
      try {
         return readCSVData (reader, colDelim, stringDelim);
      }
      finally {
         reader.close();
      }
   }

   /**
    * This is equivalent to  #readDoubleData2D(String), for reading
    * strings.
    *  @param file         the name of the file to read.
    *  @param colDelim     the column delimiter.
    *  @param stringDelim  the string delimiter.
    *  @return the obtained matrix of string values.
    *
    *  @exception IOException if an I/O error occurs.
    */
   public static String[][] readCSVData (String file, char colDelim,
                                         char stringDelim)
                                         throws IOException {
      FileReader reader = new FileReader (file);
      try {
         return readCSVData (reader, colDelim, stringDelim);
      }
      finally {
         reader.close();
      }
   }

}