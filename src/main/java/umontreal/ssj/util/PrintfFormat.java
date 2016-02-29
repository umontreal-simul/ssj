/*
 * Class:        PrintfFormat
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
package umontreal.ssj.util;

import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.util.Formatter;

/**
 * This class acts like a  StringBuffer which defines new types of `append`
 * methods. It defines certain functionalities of the ANSI C `printf`
 * function that also can be accessed through static methods. The information
 * given here is strongly inspired from the `man` page of the C `printf`
 * function.
 *
 * Most methods of this class format numbers for the English US locale only.
 * One can use the Java class  Formatter for performing locale-independent
 * formatting.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class PrintfFormat implements CharSequence, Appendable {
   private static NumberFormat nf = NumberFormat.getInstance (Locale.US);
   private static DecimalFormatSymbols dfs =
       new DecimalFormatSymbols (Locale.US);
   private static final int NDEC = 50;
   private static DecimalFormat[] dfe = new DecimalFormat[NDEC+1];
   private static DecimalFormat[] dfg = new DecimalFormat[NDEC+1];
   private StringBuffer sb;
/*
   private static void round (int b, StringBuffer num) {
      // round a real number up if it is >= 0.5
      // base is b

      // round up the fractional part
      int j = num.length() - 1;
      String bm1 = String.valueOf (b - 1);
      while (j >= 0 && num.charAt(j) == bm1.charAt(0)) {
         num.deleteCharAt(j);
         j--;
      }

      char c;
      if (j < 0)
         return;
      if (num.charAt(j) != '.') {
         c = num.charAt(j);
         ++c;
         num.replace(j, j + 1, Character.toString(c));
         return;
      }

      // round up the integral part
      j--;
      while (j >= 0 && num.charAt(j) == bm1.charAt(0)) {
         num.replace(j, j + 1, "0");
         j--;
      }

     if (j < 0 || num.charAt(j) == '-') {
         num.insert( j + 1, '1');
     } else if (num.charAt(j) == ' ') {
         num.replace(j, j + 1, "1");
     } else {
         c = num.charAt(j);
         ++c;
         num.replace(j, j + 1, Character.toString(c));
      }
  }*/

   /**
    * @name Constants
    * @{
    */

   /**
    * End-of-line symbol or line separator. It is "\\n" for Unix/Linux,
    * "\\r\\n" for MS-DOS/MS-Windows, and "\\r" for Apple OSX.
    */
   public static final String NEWLINE =
               System.getProperty("line.separator");

   /**
    * End-of-line symbol or line separator. Same as `NEWLINE`.
    */
   public static final String LINE_SEPARATOR =
               System.getProperty("line.separator");

   /**
    * @}
    */

   /**
    * Constructs a new buffer object containing an empty string.
    */
   public PrintfFormat() {
      sb = new StringBuffer();
   }

   /**
    * Constructs a new buffer object with an initial capacity of `length`.
    *  @param length       initial length of the buffer
    */
   public PrintfFormat (int length) {
      sb = new StringBuffer (length);
   }

   /**
    * Constructs a new buffer object containing the initial string `str`.
    *  @param str          initial contents of the buffer
    */
   public PrintfFormat (String str) {
      sb = new StringBuffer (str);
   }

   /**
    * Appends `str` to the buffer.
    *  @param str          string to append to the buffer
    *  @return this object
    */
   public PrintfFormat append (String str) {
      sb.append (str);
      return this;
   }

   /**
    * Uses the  #s(int,String) static method to append `str` to the
    * buffer. A minimum of `fieldwidth` characters will be used.
    *  @param fieldwidth   minimum number of characters that will be added
    *                      to the buffer
    *  @param str          string to append to the buffer
    *  @return this object
    */
   public PrintfFormat append (int fieldwidth, String str) {
      sb.append (s (fieldwidth, str));
      return this;
   }

   /**
    * Appends `x` to the buffer.
    *  @param x            value being added to the buffer
    *  @return this object
    */
   public PrintfFormat append (double x) {
      sb.append (x);
      return this;
   }

   /**
    * Uses the  #f(int,double) static method to append `x` to the buffer.
    * A minimum of `fieldwidth` characters will be used.
    *  @param fieldwidth   minimum length of the converted string to be
    *                      appended
    *  @param x            value to be appended to the buffer
    *  @return this object
    */
   public PrintfFormat append (int fieldwidth, double x) {
      sb.append (f (fieldwidth, x));
      return this;
   }

   /**
    * Uses the  #f(int,int,double) static method to append `x` to the
    * buffer. A minimum of `fieldwidth` characters will be used with the
    * given `precision`.
    *  @param fieldwidth   minimum length of the converted string to be
    *                      appended
    *  @param precision    number of digits after the decimal point of the
    *                      converted value
    *  @param x            value to be appended to the buffer
    *  @return this object
    */
   public PrintfFormat append (int fieldwidth, int precision, double x) {
      sb.append (f (fieldwidth, precision, x));
      return this;
   }

   /**
    * Appends `x` to the buffer.
    *  @param x            value to be appended to the buffer
    *  @return this object
    */
   public PrintfFormat append (int x) {
      sb.append (x);
      return this;
   }

   /**
    * Uses the  #d(int,long) static method to append `x` to the buffer. A
    * minimum of `fieldwidth` characters will be used.
    *  @param fieldwidth   minimum length of the converted string to be
    *                      appended
    *  @param x            value to be appended to the buffer
    *  @return this object
    */
   public PrintfFormat append (int fieldwidth, int x) {
      sb.append (d (fieldwidth, x));
      return this;
   }

   /**
    * Appends `x` to the buffer.
    *  @param x            value to be appended to the buffer
    *  @return this object
    */
   public PrintfFormat append (long x) {
      sb.append (x);
      return this;
   }

   /**
    * Uses the  #d(int,long) static method to append `x` to the buffer. A
    * minimum of `fieldwidth` characters will be used.
    *  @param fieldwidth   minimum length of the converted string to be
    *                      appended
    *  @param x            value to be appended to the buffer
    *  @return this object
    */
   public PrintfFormat append (int fieldwidth, long x) {
      sb.append (d (fieldwidth, x));
      return this;
   }

   /**
    * Uses the  #format(int,int,int,double) static method with the same
    * four arguments to append `x` to the buffer.
    *  @param fieldwidth   minimum length of the converted string to be
    *                      appended
    *  @param accuracy     number of digits after the decimal point
    *  @param precision    number of significant digits
    *  @param x            value to be appended to the buffer
    *  @return this object
    */
   public PrintfFormat append (int fieldwidth, int accuracy, int precision,
                               double x) {
      sb.append (format (fieldwidth, accuracy, precision, x));
      return this;
   }

   /**
    * Appends a single character to the buffer.
    *  @param c            character to be appended to the buffer
    *  @return this object
    */
   public PrintfFormat append (char c) {
      sb.append (c);
      return this;
   }

   /**
    * Clears the contents of the buffer.
    */
   public void clear() {
      sb.setLength (0);
   }

   /**
    * Returns the  StringBuffer associated with that object.
    *  @return the internal  StringBuffer object
    */
   public StringBuffer getBuffer() {
      return sb;
   }

   /**
    * Converts the buffer into a  String.
    *  @return the  String conversion of the internal buffer
    */
   public String toString() {
      return sb.toString();
   }

   /**
    * Same as  {@link #s(int,String) s(0, str)}. If the string `str` is
    * null, it returns the string "null".
    *  @param str          the string to process
    *  @return the same string
    */
   public static String s (String str) {
      if (str == null)
         return "null";
      else
         return str;
   }

   /**
    * Formats the string `str` like the <tt>%s</tt> in the C `printf`
    * function. The `fieldwidth` argument gives the minimum length of the
    * resulting string. If `str` is shorter than `fieldwidth`, it is
    * left-padded with spaces. If `fieldwidth` is negative, the string is
    * right-padded with spaces if necessary. The  String will never be
    * truncated. If `str` is null, it calls  {@link #s(int,String)
    * s(fieldwidth, "null")}. The `fieldwidth` argument has the same
    * effect for the other methods in this class.
    *  @param fieldwidth   minimum length of the returned string
    *  @param str          the string to process
    *  @return the same string padded with spaces if necessary
    */
   public static String s (int fieldwidth, String str) {
      if (str == null)
         return s (fieldwidth, "null");

      int fw = Math.abs (fieldwidth);
      if (str.length() < fw) {
         // We have to pad with spaces
         StringBuffer buf = new StringBuffer();
         int sl = str.length();
         for (int i = 0; i < fw-sl; i++)
            buf.append (' ');
         // Add the spaces before or after
         return fieldwidth >= 0 ? buf.toString() + str
                     : str + buf.toString();
      }
      else
         return str;
   }

   /**
    * @name Integers
    * @{
    */

   /**
    * Same as  {@link #d(int,int,long) d(0, 1, x)}.
    *  @param x            the string to process
    *  @return the same string, padded with spaces or zeros if appropriate
    */
   public static String d (long x) {
        return d (0, 1, x);
   }

   /**
    * Same as  {@link #d(int,int,long) d(fieldwidth, 1, x)}.
    *  @param fieldwidth   minimum length of the returned string
    *  @param x            the string to process
    *  @return the same string, padded with spaces or zeros if appropriate
    */
   public static String d (int fieldwidth, long x) {
        return d (fieldwidth, 1, x);
   }

   /**
    * Formats the long integer `x` into a string like <tt>%d</tt> in the C
    * `printf` function. It converts its argument to decimal notation,
    * `precision` gives the minimum number of digits that must appear; if
    * the converted value requires fewer digits, it is padded on the left
    * with zeros. When zero is printed with an explicit precision 0, the
    * output is empty.
    *  @param fieldwidth   minimum length of the returned string
    *  @param precision    number of digits in the returned string
    *  @param x            the string to process
    *  @return the same string, padded with spaces or zeros if appropriate
    */
   public static String d (int fieldwidth, int precision, long x) {
        if (precision < 0)
            throw new IllegalArgumentException ("precision must " +
                                               "not be negative.");
        if (precision == 0 && x == 0)
            return s (fieldwidth, "");

        nf.setGroupingUsed (false);
        nf.setMinimumIntegerDigits (precision);
        nf.setMaximumFractionDigits (0); // will also set the min to 0
        return s (fieldwidth, nf.format (x));
   }

   /**
    * Same as  {@link #d(int,int,long) d(0, 1, x)}.
    *  @param x            the value to be processed
    *  @return the string resulting from the conversion
    */
   public static String format (long x) {
      return d (0, 1, x);
   }

   /**
    * Converts a long integer to a  String with a minimum length of
    * `fieldwidth`, the result is right-padded with spaces if necessary
    * but it is not truncated. If only one argument is specified, a
    * `fieldwidth` of 0 is assumed.
    *  @param fieldwidth   minimum length of the returned string
    *  @param x            the value to be processed
    *  @return the string resulting from the conversion
    */
   public static String format (int fieldwidth, long x) {
      return d (fieldwidth, 1, x);
   }

   /**
    * Same as  {@link #formatBase(int,int,long) formatBase(0, b, x)}.
    *  @param b            the base used for conversion
    *  @param x            the value to be processed
    *  @return a string representing `x` in base `b`
    */
   public static String formatBase (int b, long x) {
      return formatBase (0, b, x);
   }

   /**
    * Converts the integer `x` to a  String representation in base `b`.
    * Restrictions: @f$2\le@f$ `b` @f$\le10@f$.
    *  @param fieldwidth   minimum length of the returned string
    *  @param b            the base used for conversion
    *  @param x            the value to be processed
    *  @return a string representing `x` in base `b`
    */
   public static String formatBase (int fieldwidth, int b, long x) {
      boolean neg = false;                   // insert a '-' if true
      if (b < 2 || b > 10)
         throw new IllegalArgumentException ("base must be between 2 and 10.");

      if (x < 0) {
         neg = true;
         x = -x;
      } else {
         if (x == 0)
            return "0";

         neg = false;
      }
      StringBuffer sb = new StringBuffer();
      while (x > 0) {
         sb.insert(0, x % b);
         x = x/b;
      }
      if (neg)
         sb.insert(0, '-');
      return s (fieldwidth, sb.toString());
   }

   /**
    * @}
    */

   /**
    * @name Reals
    * @{
    */

   /**
    * Same as  {@link #E(int,int,double) E(0, 6, x)}.
    *  @param x            the value to be converted to string
    *  @return the converted value as a string
    */
   public static String E (double x) {
        return E (0, 6, x);
   }

   /**
    * Same as  {@link #E(int,int,double) E(fieldwidth, 6, x)}.
    *  @param fieldwidth   minimum length of the returned string
    *  @param x            the value to be converted to string
    *  @return the converted value as a string
    */
   public static String E (int fieldwidth, double x) {
        return E (fieldwidth, 6, x);
   }

   /**
    * Formats a double-precision number `x` like <tt>%E</tt> in C
    * <tt>printf</tt>. The double argument is rounded and converted in the
    * style `[-]d.dddE+-dd` where there is one digit before the
    * decimal-point character and the number of digits after it is equal
    * to the precision; if the precision is 0, no decimal-point character
    * appears. The exponent always contains at least two digits; if the
    * value is zero, the exponent is `00`.
    *  @param fieldwidth   minimum length of the returned string
    *  @param precision    number of digits after the decimal point
    *  @param x            the value to be converted to string
    *  @return the converted value as a string
    */
   public static String E (int fieldwidth, int precision, double x) {
        if (precision < 0)
            throw new IllegalArgumentException ("precision must " +
                                               "not be negative.");
        if (Double.isNaN (x))
           return s (fieldwidth, "NaN");
        if (Double.isInfinite (x))
           return s (fieldwidth, (x < 0 ? "-" : "") + "Infinite");

        DecimalFormat df;
        if (precision >= dfe.length || dfe[precision] == null) {
          // We need to create one pattern per precision value
          StringBuffer pattern = new StringBuffer ("0.");
          for (int i = 0; i < precision; i++)
              pattern.append ("0");
          pattern.append ("E00");
          df = new DecimalFormat (pattern.toString(), dfs);
          df.setGroupingUsed (false);
          if (precision < dfe.length)
            dfe[precision] = df;
        }
        else
          df = dfe[precision];
        String res = df.format (x);
        // DecimalFormat doesn't add the + sign before the value of
        // the exponent.
        int exppos = res.indexOf ('E');
        if (exppos != -1 && res.charAt (exppos+1) != '-')
            res = res.substring (0, exppos+1) + "+" + res.substring (exppos+1);
        return s (fieldwidth, res);
   }

   /**
    * Same as  {@link #e(int,int,double) e(0, 6, x)}.
    *  @param x            the value to be converted to string
    *  @return the converted value as a string
    */
   public static String e (double x) {
        return e (0, 6, x);
   }

   /**
    * Same as  {@link #e(int,int,double) e(fieldwidth, 6, x)}.
    *  @param fieldwidth   minimum length of the returned string
    *  @param x            the value to be converted to string
    *  @return the converted value as a string
    */
   public static String e (int fieldwidth, double x) {
        return e (fieldwidth, 6, x);
   }

   /**
    * The same as `E`, except that `‘e’` is used as the exponent character
    * instead of `‘E’`.
    *  @param fieldwidth   minimum length of the returned string
    *  @param precision    number of digits after the decimal point
    *  @param x            the value to be converted to string
    *  @return the converted value as a string
    */
   public static String e (int fieldwidth, int precision, double x) {
        String res = E (fieldwidth, precision, x);
        int exppos = res.indexOf ('E');
        return exppos == -1 ? res : res.substring (0,
                               exppos) + 'e' + res.substring (exppos+1);
   }

   /**
    * Same as  {@link #f(int,int,double) f(0, 6, x)}.
    *  @param x            the value to be converted to string
    *  @return the converted value as a string
    */
   public static String f (double x) {
        return f (0, 6, x);
   }

   /**
    * Same as  {@link #f(int,int,double) f(fieldwidth, 6, x)}.
    *  @param fieldwidth   minimum length of the returned string
    *  @param x            the value to be converted to string
    *  @return the converted value as a string
    */
   public static String f (int fieldwidth, double x) {
        return f (fieldwidth, 6, x);
   }

   /**
    * Formats the double-precision `x` into a string like <tt>%f</tt> in C
    * `printf`. The argument is rounded and converted to decimal notation
    * in the style `[-]ddd.ddd`, where the number of digits after the
    * decimal-point character is equal to the precision specification. If
    * the precision is explicitly 0, no decimal-point character appears.
    * If a decimal point appears, at least one digit appears before it.
    *  @param fieldwidth   minimum length of the returned string
    *  @param precision    number of digits after the decimal point
    *  @param x            the value to be converted to string
    *  @return the converted value as a string
    */
   public static String f (int fieldwidth, int precision, double x) {
        if (precision < 0)
            throw new IllegalArgumentException ("precision must " +
                                               "not be negative.");
        if (Double.isNaN (x))
           return s (fieldwidth, "NaN");
        if (Double.isInfinite (x))
           return s (fieldwidth, (x < 0 ? "-" : "" ) + "Infinite");

        nf.setGroupingUsed (false);
        nf.setMinimumIntegerDigits (1);
        nf.setMinimumFractionDigits (precision);
        nf.setMaximumFractionDigits (precision);
        return s (fieldwidth, nf.format (x));
   }

   /**
    * Same as  {@link #G(int,int,double) G(0, 6, x)}.
    *  @param x            the value to be converted to string
    *  @return the converted value as a string
    */
   public static String G (double x) {
        return G (0, 6, x);
   }

   /**
    * Same as  {@link #G(int,int,double) G(fieldwidth, 6, x)}.
    *  @param fieldwidth   minimum length of the returned string
    *  @param x            the value to be converted to string
    *  @return the converted value as a string
    */
   public static String G (int fieldwidth, double x) {
        return G (fieldwidth, 6, x);
   }

   /**
    * Formats the double-precision `x` into a string like <tt>%G</tt> in C
    * `printf`. The argument is converted in style <tt>%f</tt> or
    * <tt>%E</tt>. `precision` specifies the number of significant digits.
    * If it is 0, it is treated as 1. Style <tt>%E</tt> is used if the
    * exponent from its conversion is less than @f$-4@f$ or greater than
    * or equal to `precision`. Trailing zeros are removed from the
    * fractional part of the result; a decimal point appears only if it is
    * followed by at least one digit.
    *  @param fieldwidth   minimum length of the returned string
    *  @param precision    number of significant digits
    *  @param x            the value to be converted to string
    *  @return the converted value as a string
    */
   public static String G (int fieldwidth, int precision, double x) {
        if (precision < 0)
            throw new IllegalArgumentException ("precision must " +
                                               "not be negative.");
        if (precision == 0)
            precision = 1;

        if (Double.isNaN (x))
           return s (fieldwidth, "NaN");
        if (Double.isInfinite (x))
           return s (fieldwidth, (x < 0 ? "-" : "" ) + "Infinite");

        // Calculate the scientific notation format.
        // We cannot use E because it can make trailing zeros
        // that must be removed afterward.
        DecimalFormat df;
        if (precision >= dfg.length || dfg[precision] == null) {
          StringBuffer pattern = new StringBuffer ("0.");
          for (int i = 0; i < (precision-1); i++)
              pattern.append ("#"); // Do not show trailing zeros
          pattern.append ("E00");
          df = new DecimalFormat (pattern.toString(), dfs);
          df.setGroupingUsed (false);
          if (precision < dfg.length)
            dfg[precision] = df;
        }
        else
          df = dfg[precision];
        String res = df.format (x);

        int exppos = res.indexOf ('E');
        if (exppos == -1)
           return res;
        int expval = Integer.parseInt (res.substring (exppos+1));
        if (expval < -4 || expval >= precision) {
           // add the plus sign for the exponent if necessary.
           if (res.charAt (exppos+1) != '-')
               return s (fieldwidth, res.substring (0, exppos+1) + "+" +
                          res.substring (exppos+1));
           else
               return s (fieldwidth, res);
        }

        // Calculate the decimal notation format
        nf.setGroupingUsed (false);
        nf.setMinimumIntegerDigits (1);
        nf.setMinimumFractionDigits (0);
        // We need the number of digits after the decimal point to
        // to get precisions significant digits.
        // The integer part of x contains at most precision digits.
        // If that was not true, expval would be greater than precision.
        // If expval=0, we have a number of the form 1.234...
        // To have precisions significant digits, we need precision-1 digits.
        // If expval<0, x<1 and we need -expval additionnal
        // decimal digits.
        // If expval>0, we need less decimal digits.
        nf.setMaximumFractionDigits (precision-expval-1);
        res = nf.format (x);
        return s (fieldwidth, res);
   }

   /**
    * Same as  {@link #g(int,int,double) g(0, 6, x)}.
    *  @param x            the value to be converted to string
    *  @return the converted value as a string
    */
   public static String g (double x) {
        return g (0, 6, x);
   }

   /**
    * Same as  {@link #g(int,int,double) g(fieldwidth, 6, x)}.
    *  @param fieldwidth   minimum length of the returned string
    *  @param x            the value to be converted to string
    *  @return the converted value as a string
    */
   public static String g (int fieldwidth, double x) {
        return g (fieldwidth, 6, x);
   }

   /**
    * The same as `G`, except that `‘e’` is used in the scientific
    * notation.
    *  @param fieldwidth   minimum length of the returned string
    *  @param precision    number of significant digits
    *  @param x            the value to be converted to string
    *  @return the converted value as a string
    */
   public static String g (int fieldwidth, int precision, double x) {
        String res = G (fieldwidth, precision, x);
        int exppos = res.indexOf ('E');
        return exppos == -1 ? res :
            res.substring (0, exppos) + 'e' + res.substring (exppos+1);
   }

   /**
    * Returns a  String containing `x`. Uses a total of at least
    * `fieldwidth` positions (including the sign and point when they
    * appear), `accuracy` digits after the decimal point and at least
    * `precision` significant digits. `accuracy` and `precision` must be
    * strictly smaller than `fieldwidth`. The number is rounded if
    * necessary. If there is not enough space to format the number in
    * decimal notation with at least `precision` significant digits
    * (<tt>accuracy</tt> or `fieldwidth` is too small), it will be
    * converted to scientific notation with at least `precision`
    * significant digits. In that case, `fieldwidth` is increased if
    * necessary.
    *  @param fieldwidth   minimum length of the returned string
    *  @param accuracy     number of digits after the decimal point
    *  @param precision    number of significant digits
    *  @param x            the value to be processed
    *  @return the converted value as a string
    */
   public static String format (int fieldwidth, int accuracy, int precision,
                                double x) {
        if (Double.isNaN (x))
           return s (fieldwidth, "NaN");
        if (Double.isInfinite (x))
           return s (fieldwidth, (x < 0 ? "-" : "" ) + "Infinite");

       if (canUseDecimalNotation (fieldwidth, accuracy, precision, x))
          return f (fieldwidth, accuracy, x);
       // Use scientific notation
       else {
          String S = E (fieldwidth, precision - 1 , x);
          return processExp (S);
       }
   }

   private static boolean canUseDecimalNotation (int fieldwidth, int accuracy,
                                                 int precision, double x) {
      // Le nombre de positions occupees par la partie entiere de x
      int PosEntier = 0;
      // Le nombre de chiffres significatifs avant le point
      int EntierSign;
      // La position de l'exposant dans le string S et la longueur de S
      int Neg = 0;

      if (x == 0.0)
         EntierSign = 1;
      else {
         EntierSign = PosEntier = (int)Math.floor (
               Math.log10 (Math.abs (x)) + 1);
         if (x < 0.0)
             Neg = 1;
      }
      if (EntierSign <= 0)
          PosEntier = 1;
      return x == 0.0 || (((EntierSign + accuracy) >= precision) &&
                          (fieldwidth >= (PosEntier + accuracy + Neg + 1)));
   }

   private static int getMinAccuracy (double x) {
      if (Math.abs (x) >= 1 || x == 0)
         return 0;
      else
         return -(int)Math.floor (Math.log10 (Math.abs (x)));
   }

   private static String processExp (String s) {
      int p = s.indexOf ("E+0");
      if (p == -1)
         p = s.indexOf ("E-0");

      // remove the 0 in E-0 and in E+0
      if (p != -1)
         s = " " + s.substring (0, p + 2) + s.substring (p + 3);

      p = s.indexOf (".E");
      if (p != -1)
         s = " " + s.substring (0, p) + s.substring (p + 1);
      return s;
   }

   /**
    * This method is equivalent to  #format(int,int,int,double), except it
    * formats the given value for the locale `locale`.
    *  @param locale       the locale being used for formatting
    *  @param fieldwidth   minimum length of the returned string
    *  @param accuracy     number of digits after the decimal point
    *  @param precision    number of significant digits
    *  @param x            the value to be processed
    *  @return the converted value as a string
    */
   public static String format (Locale locale, int fieldwidth, int accuracy,
                                int precision, double x) {
         Formatter fmt = new Formatter (locale);
        if (Double.isNaN (x))
           return fmt.format ("%" + fieldwidth + "s", "NaN").toString();
        if (Double.isInfinite (x))
           return fmt.format ("%" + fieldwidth + "s", (x < 0 ? "-" : "" ) + "Infinite").toString();

       if (canUseDecimalNotation (fieldwidth, accuracy, precision, x))
          return fmt.format ("%" + fieldwidth + "." + accuracy + "f", x).toString();
       // Use scientific notation
       else {
          String S = fmt.format ("%" + fieldwidth + "." + (precision - 1) + "E", x).toString();
          return processExp (S);
       }
   }

   /**
    * Converts @f$x@f$ to a String representation in base @f$b@f$ using
    * formatting similar to the @f$f@f$ methods. Uses a total of at least
    * `fieldwidth` positions (including the sign and point when they
    * appear) and `accuracy` digits after the decimal point. If
    * `fieldwidth` is negative, the number is printed left-justified,
    * otherwise right-justified. Restrictions: @f$2 \le b \le10@f$ and
    * @f$|x| < 2^{63}@f$.
    *  @param fieldwidth   minimum length of the returned string
    *  @param accuracy     number of digits after the decimal point
    *  @param b            base
    *  @param x            the value to be processed
    *  @return the converted value as a string
    */
   public static String formatBase (int fieldwidth, int accuracy, int b,
                                    double x) {
      if (Double.isNaN (x))
         return s (fieldwidth, "NaN");
      if (Double.isInfinite (x))
         return s (fieldwidth, (x < 0 ? "-" : "" ) + "Infinite");
      if (0. == x || -0. == x)
         return s (fieldwidth, "0");
      if (Math.abs(x) >= Num.TWOEXP[63])
         throw new UnsupportedOperationException ("   |x| >= 2^63");

      long n = (long) x;
      String mant = formatBase(-1, b, n);
      if (n == x)
         return s (fieldwidth, mant);
      if (n == 0) {
         if (x < 0) {
            mant = "-0";
         } else
            mant = "0";
      }
      // round before printing
      if (x > 0)
         x += 0.5*Math.pow(b, -accuracy - 1);
      else if (x < 0)
         x -= 0.5*Math.pow(b, -accuracy - 1);
       x -= n;
      if (x < 0)
         x = -x;

      StringBuffer frac = new StringBuffer(".");
      long y;
      int j;
      for (j = 0; j < accuracy; ++j) {
         x *= b;
         y = (long) x;
         frac.append(y);
         x -= y;
         if (x == 0.)
            break;
      }

      StringBuffer number = new StringBuffer(mant);
      number.append(frac);

      // remove trailing spaces and 0
      j = number.length() - 1;
      while (j > 0 && (number.charAt(j) == '0' || number.charAt(j) == ' ' )) {
         number.deleteCharAt(j);
         j--;
      }

      return s (fieldwidth, number.toString());
   }
   // Interface CharSequence
   public char charAt (int index) {
      return sb.charAt (index);
   }

   public int length() {
      return sb.length();
   }

   public CharSequence subSequence (int start, int end) {
      return sb.subSequence (start, end);
   }

   // Interface Appendable
   public Appendable append (CharSequence csq) {
      return sb.append (csq);
   }

   public Appendable append (CharSequence csq, int start, int end) {
      return sb.append (csq, start, end);
   }

   /**
    * @}
    */

   /**
    * @name Intervals
    * @{
    */

   /**
    * Stores a string containing `x` into `res[0]`, and a string
    * containing `error` into `res[1]`, both strings being formatted with
    * the same notation. Uses a total of at least `fieldwidth` positions
    * (including the sign and point when they appear) for `x`,
    * `fieldwidtherr` positions for `error`, `accuracy` digits after the
    * decimal point and at least `precision` significant digits.
    * `accuracy` and `precision` must be strictly smaller than
    * `fieldwidth`. The numbers are rounded if necessary. If there is not
    * enough space to format `x` in decimal notation with at least
    * `precision` significant digits (<tt>accuracy</tt> or `fieldwidth`
    * are too small), it will be converted to scientific notation with at
    * least `precision` significant digits. In that case, `fieldwidth` is
    * increased if necessary, and the error is also formatted in
    * scientific notation.
    *  @param fieldwidth   minimum length of the value string
    *  @param fieldwidtherr minimum length of the error string
    *  @param accuracy     number of digits after the decimal point for
    *                      the value and error
    *  @param precision    number of significant digits for the value
    *  @param x            the value to be processed
    *  @param error        the error on the value to be processed
    *  @param res          an array that will be filled with the formatted
    *                      value and formatted error
    */
   public static void formatWithError (int fieldwidth, int fieldwidtherr,
          int accuracy, int precision, double x, double error, String[] res) {
      if (res.length != 2)
         throw new IllegalArgumentException ("The given res array must contain two elements");
      if (Double.isNaN (x)) {
         res[0] = s (fieldwidth, "NaN");
         res[1] = s (fieldwidtherr, "");
         return;
      }
      if (Double.isInfinite (x)) {
         res[0] = s (fieldwidth, (x < 0 ? "-" : "" ) + "Infinite");
         res[1] = s (fieldwidtherr, "");
         return;
      }
      if (accuracy < 0)
         accuracy = 0;

      if (canUseDecimalNotation (fieldwidth, accuracy, precision, x)) {
         res[0] = f (fieldwidth, accuracy, x);
         res[1] = f (fieldwidtherr, accuracy, error);
      }
      // Use scientific notation
      else {
         res[0] = processExp (E (fieldwidth, precision - 1, x));
         int xExp = x == 0 ? 0 : (int)Math.floor (Math.log10 (Math.abs (x)));
         int errorExp = error == 0 ? 0 : (int)Math.floor (Math.log10 (Math.abs (error)));
         int errorPrecision = precision - 1 - (xExp - errorExp);
         if (errorPrecision < 0)
            errorPrecision = 0;
         res[1] = processExp (E (fieldwidtherr, errorPrecision, error));
      }
   }

   /**
    * Stores a string containing `x` into `res[0]`, and a string
    * containing `error` into `res[1]`, both strings being formatted with
    * the same notation. This calls
    * #formatWithError(int,int,int,int,double,double,String[]) with the
    * minimal accuracy for which the formatted string for `error` is
    * non-zero. If `error` is 0, the accuracy is 0. If this minimal
    * accuracy causes the strings to be formatted using scientific
    * notation, this method increases the accuracy until the decimal
    * notation can be used.
    *  @param fieldwidth   minimum length of the value string
    *  @param fieldwidtherr minimum length of the error string
    *  @param precision    number of significant digits for the value
    *  @param x            the value to be processed
    *  @param error        the error on the value to be processed
    *  @param res          an array that will be filled with the formatted
    *                      value and formatted error
    */
   public static void formatWithError (int fieldwidth, int fieldwidtherr,
          int precision, double x, double error, String[] res) {
      int accuracy = getMinAccuracy (error);
      if (!canUseDecimalNotation (fieldwidth, accuracy, precision, x)) {
         int posEntier = (int)Math.floor (Math.log (Math.abs (x)) / Math.log (10) + 1);
         if (posEntier < 0)
            posEntier = 1;
         int newAccuracy = precision - posEntier;
         if (canUseDecimalNotation (fieldwidth, newAccuracy, precision, x))
            accuracy = newAccuracy;
      }
      formatWithError (fieldwidth, fieldwidtherr, accuracy, precision, x, error, res);
   }

   /**
    * This method is equivalent to
    * #formatWithError(int,int,int,double,double,String[]), except that it
    * formats the given value and error for the locale `locale`.
    *  @param locale       the locale being used
    *  @param fieldwidth   minimum length of the value string
    *  @param fieldwidtherr minimum length of the error string
    *  @param accuracy     number of digits after the decimal point for
    *                      the value and error
    *  @param precision    number of significant digits for the value
    *  @param x            the value to be processed
    *  @param error        the error on the value to be processed
    *  @param res          an array that will be filled with the formatted
    *                      value and formatted error
    */
   public static void formatWithError (Locale locale, int fieldwidth,
          int fieldwidtherr, int accuracy, int precision, double x,
          double error, String[] res) {
      if (res.length != 2)
         throw new IllegalArgumentException ("The given res array must contain two elements");
      Formatter fmt = new Formatter (locale);
      Formatter fmtErr = new Formatter (locale);
      if (Double.isNaN (x)) {
         res[0] = fmt.format ("%" + fieldwidth + "s", "NaN").toString();
         res[1] = fmtErr.format ("%" + fieldwidtherr + "s", "").toString();
         return;
      }
      if (Double.isInfinite (x)) {
         res[0] = fmt.format ("%" + fieldwidth + "s", (x < 0 ? "-" : "" ) + "Infinite").toString();
         res[1] = fmtErr.format ("%" + fieldwidtherr + "s", "").toString();
         return;
      }
      if (accuracy < 0)
         accuracy = 0;

      if (canUseDecimalNotation (fieldwidth, accuracy, precision, x)) {
         res[0] = fmt.format ("%" + fieldwidth + "." + accuracy + "f", x).toString();
         res[1] = fmtErr.format ("%" + fieldwidtherr + "." + accuracy + "f", error).toString();
      }
      // Use scientific notation
      else {
         res[0] = processExp (fmt.format ("%" + fieldwidth + "." + (precision - 1) + "E", x).toString());
         int xExp = x == 0 ? 0 : (int)Math.floor (Math.log10 (Math.abs (x)));
         int errorExp = error == 0 ? 0 : (int)Math.floor (Math.log10 (Math.abs (error)));
         int errorPrecision = precision - 1 - (xExp - errorExp);
         if (errorPrecision < 0)
            errorPrecision = 0;
         res[1] = processExp (fmtErr.format
         ("%" + fieldwidtherr + "." + errorPrecision + "E", error).toString());
      }
   }

   /**
    * This method is equivalent to
    * #formatWithError(int,int,int,double,double,String[]), except that it
    * formats the given value and error for the locale `locale`.
    *  @param locale       the locale being used
    *  @param fieldwidth   minimum length of the value string
    *  @param fieldwidtherr minimum length of the error string
    *  @param precision    number of significant digits for the value
    *  @param x            the value to be processed
    *  @param error        the error on the value to be processed
    *  @param res          an array that will be filled with the formatted
    *                      value and formatted error
    */
   public static void formatWithError (Locale locale, int fieldwidth,
          int fieldwidtherr, int precision, double x, double error,
          String[] res) {
      int accuracy = getMinAccuracy (error);
      if (!canUseDecimalNotation (fieldwidth, accuracy, precision, x)) {
         int posEntier = (int)Math.floor (Math.log (Math.abs (x)) / Math.log (10) + 1);
         if (posEntier < 0)
            posEntier = 1;
         int newAccuracy = precision - posEntier;
         if (canUseDecimalNotation (fieldwidth, newAccuracy, precision, x))
            accuracy = newAccuracy;
      }
      formatWithError (locale, fieldwidth, fieldwidtherr, accuracy, precision, x, error, res);
   }

}

/**
 * @}
 */