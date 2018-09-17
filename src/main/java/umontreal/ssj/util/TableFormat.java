/*
 * Class:        TableFormat
 * Description:  Provides methods to format arrays into String's
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

/**
 * This class provides methods to format arrays and matrices into  String s
 * in different styles. This could be useful for printing arrays and
 * subarrays, or for putting them in files for further treatment by other
 * softwares such as *Mathematica*, *Matlab*, etc.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class TableFormat {
   private TableFormat() {}

   /**
    * @name Formating styles
    * @{
    */

   /**
    * Plain text matrix printing style
    */
   public static final int PLAIN = 0;

   /**
    * Mathematica matrix printing style
    */
   public static final int MATHEMATICA = 1;

   /**
    * Matlab matrix printing style
    */
   public static final int MATLAB = 2;

   /**
    * @}
    */

   /**
    * @name Functions to convert arrays to `String`
    * @{
    */

   /**
    * Formats a  String containing the elements `n1` to `n2` (inclusive)
    * of table `V`, `k` elements per line, `p` positions per element. If
    * `k` = 1, the array index will also appear on the left of each
    * element, i.e., each line `i` will have the form `i V[i]`.
    *  @param V            array to be formated
    *  @param n1           index of the first element being formated
    *  @param n2           index of the last element being formated
    *  @param k            number of elements per line
    *  @param p            number of positions per element
    *  @return formated string repreenting the elements
    */
   public static String format (int V[], int n1, int n2, int k, int p) {
      int i;
      StringBuffer sb = new StringBuffer();
      if (k > 1) {
         sb.append ("Elements  " + n1 + "  to  " + n2 +
                     PrintfFormat.NEWLINE + PrintfFormat.NEWLINE);
         for (i = n1; i <= n2; i++) {
            sb.append (PrintfFormat.d (p, V[i]));
            if (((i + 1 - n1) % k) == 0)
               sb.append (PrintfFormat.NEWLINE);
         }
         sb.append (PrintfFormat.NEWLINE);
      }
      else {
         sb.append (PrintfFormat.NEWLINE + " Index        Element" +
                    PrintfFormat.NEWLINE);
         for (i = n1; i <= n2; i++)
            sb.append (PrintfFormat.d (6, i) + "   " +
                      PrintfFormat.d (12, V[i]) + PrintfFormat.NEWLINE);
      }
      sb.append (PrintfFormat.NEWLINE);
      return sb.toString();
   }

   /**
    * Similar to the previous method, but for an array of
    * <tt>double</tt>’s. Gives at least `p1` positions per element, `p2`
    * digits after the decimal point, and at least `p3` significant
    * digits.
    *  @param V            array to be formated
    *  @param n1           index of the first element being formated
    *  @param n2           index of the last element being formated
    *  @param k            number of elements per line
    *  @param p1           number of positions per element
    *  @param p2           number of digits after the decimal point
    *  @param p3           number of significant digits
    *  @return formated string repreenting the elements
    */
   public static String format (double V[], int n1, int n2,
                                int k, int p1, int p2, int p3) {
      int i;
      StringBuffer sb = new StringBuffer();
      if (k > 1) {
         sb.append ("Elements  " + n1 + "  to  " + n2  +
                     PrintfFormat.NEWLINE + PrintfFormat.NEWLINE);
         for (i = n1; i <= n2; i++) {
            sb.append (PrintfFormat.format (p1, p2, p3, V[i]));
            if (((i + 1 - n1) % k) == 0)
               sb.append (PrintfFormat.NEWLINE);
         }
         sb.append (PrintfFormat.NEWLINE);

      } else {
         sb.append (PrintfFormat.NEWLINE + " Index            Element" +
                    PrintfFormat.NEWLINE);
         for (i = n1; i <= n2; i++)
            sb.append (PrintfFormat.d (6, i) + "   " +
                       PrintfFormat.format (p1, p2, p3, V[i]) +
                       PrintfFormat.NEWLINE);
      }
      sb.append (PrintfFormat.NEWLINE);
      return sb.toString();
   }

   private static int Style = PLAIN;

   private static char OuvrantMat = ' ';     // Matrix delimitors
   private static char FermantMat = ' ';

   private static char OuvrantVec = ' ';     // Vector delimitors
   private static char FermantVec = ' ';

   private static char SepareVec = ' ';      // Element separators
   private static char SepareElem = ' ';

   private static void fixeDelim (int style) {
      /* Fixe les delimiteurs pour imprimer une matrice selon un format
         approprie */
      Style = style;
      switch (style) {
      case MATHEMATICA:
         OuvrantMat = '{';
         FermantMat = '}';
         OuvrantVec = '{';
         FermantVec = '}';
         SepareVec = ',';
         SepareElem = ',';
         break;
      case MATLAB:
         OuvrantMat = '[';
         FermantMat = ']';
         OuvrantVec = ' ';
         FermantVec = ' ';
         SepareVec = ' ';
         SepareElem = ' ';
         break;
      default:
         OuvrantMat = ' ';
         FermantMat = ' ';
         OuvrantVec = ' ';
         FermantVec = ' ';
         SepareVec = ' ';
         SepareElem = ' ';
         break;
      }
   }


   @Deprecated
   public static String format (int[][] Mat, int i1, int i2,
                                int j1, int j2, int w, int p,
                                int style, String Name) {
      return format (Mat, i1, i2, j1, j2, w, style, Name);
   }

/**
 * Formats the submatrix with lines `i1` @f$\le i \le@f$ `i2` and columns
 * `j1` @f$\le j \le@f$ `j2` of the matrix `Mat`, using the formatting style
 * `style`. The elements are formated in `w` positions each, with a precision
 * of `p` digits. `Name` provides an identifier for the submatrix. To be
 * treated by `Matlab`, the returned string must be copied to a file with
 * extension `.m`. If the file is named `poil.m`, for example, it can be
 * accessed by calling `poil` in `Matlab`. For `Mathematica`, if the file is
 * named `poil`, it will be read using `<< poil;`.
 *  @param Mat          matrix to be formated
 *  @param i1           index of the first row being formated
 *  @param i2           index of the last row being formated
 *  @param j1           index of the first column being formated
 *  @param j2           index of the last column being formated
 *  @param w            number of positions for each element
 *  @param p            number of digits after the decimal point of the
 *                      elements
 *  @param style        formating style of the submatrix, being one of
 *                      #MATHEMATICA,  #MATLAB, or  #PLAIN
 *  @param Name         descriptive name of the submatrix
 *  @return formated string representing the submatrix
 */
public static String format (double[][] Mat, int i1, int i2,
                                int j1, int j2, int w, int p,
                                int style, String Name) {
      int k;
      int j;
      int i;
      double x;
      String S;

      fixeDelim (style);
      StringBuffer sb = new StringBuffer();
      if (Name.length() > 0)
         sb.append (Name + " = ");

      double prec = Math.pow (10.0, (double)p);
      sb.append (OuvrantMat + PrintfFormat.NEWLINE);
      for (i = i1; i <= i2; i++) {
         sb.append (OuvrantVec);
         for (j = j1; j <= j2; j++) {
            sb.append (' ');
            switch (style) {
            case MATHEMATICA:
               x = Mat[i][j];
               if (((x != 0.0) && (Math.abs (x) < 0.1)) ||
                   (Math.abs (x) > prec)) {
                  S = PrintfFormat.G (0, p, x);
                  int exppos = S.indexOf ('E');
                  if (exppos != -1)
                     S = S.substring (0, exppos) + "*10^(" +
                          S.substring (exppos+1) + ")";
               }
               else
                  S = PrintfFormat.f (0, p, x);
               S = PrintfFormat.s (w, S);
               break;
            default:
               // MATLAB, Default */
               sb.append (PrintfFormat.G (w, p, Mat[i][j]));
               break;
            }
            if (j < j2)
               sb.append (SepareElem);
         }
         sb.append (FermantVec);
         if (i < i2)
            sb.append (SepareVec + PrintfFormat.NEWLINE);
      }
      sb.append (FermantMat + PrintfFormat.NEWLINE);
      return sb.toString();
   }

   /**
    * Similar to the previous method, but for a matrix of <tt>int</tt>’s.
    *  @param Mat          matrix to be formated
    *  @param i1           index of the first row being formated
    *  @param i2           index of the last row being formated
    *  @param j1           index of the first column being formated
    *  @param j2           index of the last column being formated
    *  @param w            number of positions for each element
    *  @param style        formating style of the submatrix, being one of
    *                      #MATHEMATICA,  #MATLAB, or  #PLAIN
    *  @param Name         descriptive name of the submatrix
    *  @return formated string representing the submatrix
    */
   public static String format (int[][] Mat, int i1, int i2, int j1, int j2,
                                int w, int style, String Name) {
      int i;
      int j;

      fixeDelim (style);
      StringBuffer sb = new StringBuffer();
      if (Name.length() > 0)
         sb.append (Name + " = ");

      sb.append (OuvrantMat + PrintfFormat.NEWLINE);
      for (i = i1; i <= i2; i++) {
         sb.append (OuvrantVec);
         for (j = j1; j <= j2; j++) {
            sb.append (PrintfFormat.d (w, Mat[i][j]));
            if (j < j2)
               sb.append (SepareElem);
         }
         sb.append (FermantVec);
         if (i < i2)
            sb.append (SepareVec + PrintfFormat.NEWLINE);
      }
      sb.append (FermantMat + PrintfFormat.NEWLINE + PrintfFormat.NEWLINE);
      return sb.toString();
   }
 }

/**
 * @}
 */