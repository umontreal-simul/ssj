/*
 * Class:        DigitalNetFromFile
 * Description:  read the parameters defining a digital net from a file
                 or from a URL address
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
package umontreal.ssj.hups;

import java.io.*;
import java.util.*;
import java.net.URL;
import java.net.MalformedURLException;
import umontreal.ssj.util.PrintfFormat;

/**
 * This class allows us to read the parameters defining a digital net either
 * from a file, or from a URL address on the World Wide Web. The parameters
 * used in building the net are those defined in class
 * @ref umontreal.ssj.hups.DigitalNet. The format of the data files must be
 * the following:
 *  (see the format in @ref umontreal.ssj.hups)
 *
 * <center>
 *  <tt><div class="SSJ-fbox">
 * <table class="SSJ-table">
 * <tr>
 *   <td colspan="4" class="l">// Any number of comment lines starting with
 * //</td>
 * </tr><tr>
 *   <td class="l">@f$b@f$</td>
 *   <td class="l"></td>
 *   <td class="l"></td>
 *   <td class="l">// @f$\mbox{Base}@f$</td>
 * </tr><tr>
 *   <td class="l">@f$k@f$</td>
 *   <td class="l"></td>
 *   <td class="l"></td>
 *   <td class="l">// Number of columns</td>
 * </tr><tr>
 *   <td class="l">@f$r@f$</td>
 *   <td class="l"></td>
 *   <td class="l"></td>
 *   <td class="l">// Maximal number of rows</td>
 * </tr><tr>
 *   <td class="l">@f$n@f$</td>
 *   <td class="l"></td>
 *   <td class="l"></td>
 *   <td class="l">// Number of points = @f$b^k@f$</td>
 * </tr><tr>
 *   <td class="l">@f$s@f$</td>
 *   <td class="l"></td>
 *   <td class="l"></td>
 *   <td class="l">// Maximal dimension of points</td>
 * </tr><tr>
 *   <td colspan="4" class="l">// dim = 1</td>
 * </tr><tr>
 *   <td class="l">@f$c_{11}@f$</td>
 *   <td class="l">@f$c_{21}@f$</td>
 *   <td class="l">@f$\cdots@f$</td>
 *   <td class="l">@f$c_{r1}@f$</td>
 * </tr><tr>
 *   <td class="l">@f$c_{12}@f$</td>
 *   <td class="l">@f$c_{22}@f$</td>
 *   <td class="l">@f$\cdots@f$</td>
 *   <td class="l">@f$c_{r2}@f$</td>
 * </tr><tr>
 *   <td class="l"></td>
 *   <td class="l">@f$\vdots@f$</td>
 *   <td class="l"></td>
 *   <td class="l"></td>
 * </tr><tr>
 *   <td class="l">@f$c_{1k}@f$</td>
 *   <td class="l">@f$c_{2k}@f$</td>
 *   <td class="l">@f$\cdots@f$</td>
 *   <td class="l">@f$c_{rk}@f$</td>
 * </tr><tr>
 *   <td colspan="4" class="l">// dim = 2</td>
 * </tr><tr>
 *   <td class="l"></td>
 *   <td class="l">@f$\vdots@f$</td>
 *   <td class="l"></td>
 *   <td class="l"></td>
 * </tr><tr>
 *   <td colspan="4" class="l">// dim = @f$s@f$</td>
 * </tr><tr>
 *   <td class="l">@f$c_{11}@f$</td>
 *   <td class="l">@f$c_{21}@f$</td>
 *   <td class="l">@f$\cdots@f$</td>
 *   <td class="l">@f$c_{r1}@f$</td>
 * </tr><tr>
 *   <td class="l">@f$c_{12}@f$</td>
 *   <td class="l">@f$c_{22}@f$</td>
 *   <td class="l">@f$\cdots@f$</td>
 *   <td class="l">@f$c_{r2}@f$</td>
 * </tr><tr>
 *   <td class="l"></td>
 *   <td class="l">@f$\vdots@f$</td>
 *   <td class="l"></td>
 *   <td class="l"></td>
 * </tr><tr>
 *   <td class="l">@f$c_{1k}@f$</td>
 *   <td class="l">@f$c_{2k}@f$</td>
 *   <td class="l">@f$\cdots@f$</td>
 *   <td class="l">@f$c_{rk}@f$</td>
 * </tr>
 * </table>
 *  </div> </tt>
 * </center>
 *
 * The figure above gives the general format of the data file needed by
 * `DigitalNetFromFile`. The values of the parameters on the left must appear
 * in the file as integers. On the right of each parameter, there is an
 * optional comment that is disregarded by the reader program. In general,
 * the Java line comments `//` are accepted anywhere and will ensure that the
 * rest of the line is dropped by the reader. Blank lines are also
 * disregarded by the reader program. For each dimension, there must be a
 * @f$k\times r@f$ matrix of integers in @f$\{0, 1, …, b-1\}@f$ (note that
 * the matrices must appear in transposed form).
 *
 * The predefined files of parameters are kept in different directories,
 * depending on the criteria used in the searches for the parameters defining
 * the digital net. These files have all been stored at the address
 * [http://simul.iro.umontreal.ca/ssj/data](http://simul.iro.umontreal.ca/ssj/data).
 * Each file contains the parameters for a specific digital net. The name of
 * the files gives information about the main parameters of the digital net.
 * For example, the file named `Edel/OOA2/B3S13R9C9St6` contains the
 * parameters for a digital net proposed by Yves Edel (see
 * [http://www.mathi.uni-heidelberg.de/~yves/index.html](http://www.mathi.uni-heidelberg.de/~yves/index.html))
 * based on ordered orthogonal arrays; the digital net has base `B = 3`,
 * dimension `S = 13`, the generating matrices have `R = 9` rows and `C = 9`
 * columns, and the strength of the net is `St = 6`.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class DigitalNetFromFile extends DigitalNet {
   private String filename;

   private void readMatrices (StreamTokenizer st,
                              int r, int k, int dim)
      throws IOException, NumberFormatException {
      // Read dim matrices with r rows and k columns.
      // dim is the dimension of the digital net.
      genMat = new int[dim * k][r];
      for (int i = 0; i < dim; i++)
         for (int c = 0; c < k; c++) {
             for (int j = 0; j < r; j++) {
                 st.nextToken ();
                 genMat[i*numCols + c][j]  = (int) st.nval;
             }
             // If we do not use all the rows, drop the unused ones.
             for (int j = r; j < numRows; j++) {
                 st.nextToken ();
             }
         }
   }


   void readData (StreamTokenizer st) throws
                                      IOException, NumberFormatException
   {
      // Read beginning of data file, but do not read the matrices
      st.eolIsSignificant (false);
      st.slashSlashComments (true);
      int i = st.nextToken ();
      if (i != StreamTokenizer.TT_NUMBER)
         throw new NumberFormatException(" readData: cannot read base");
      b = (int) st.nval;
      st.nextToken ();   numCols = (int) st.nval;
      st.nextToken ();   numRows = (int) st.nval;
      st.nextToken ();   numPoints = (int) st.nval;
      st.nextToken ();   dim = (int) st.nval;
      if (dim < 1)
         throw new IllegalArgumentException (" dimension dim <= 0");
   }


   static BufferedReader openURL (String filename)
                                  throws MalformedURLException, IOException {
      try {
         URL url = new URL (filename);
         BufferedReader input = new BufferedReader (
                                    new InputStreamReader (
                                        url.openStream()));
         return input;

      } catch (MalformedURLException e) {
         System.err.println (e + "   Invalid URL address:   " + filename);
         throw e;

      }  catch (IOException e) {
          // This can receive a FileNotFoundException
         System.err.println (e + " in openURL with " + filename);
         throw e;
      }
   }

   static BufferedReader openFile (String filename) throws
            IOException {
      try {
         BufferedReader input;
         File f = new File (filename);

         // If file with relative path name exists, read it
         if (f.exists()) {
            if (f.isDirectory())
               throw new IOException (filename + " is a directory");
            input = new BufferedReader (new FileReader (filename));
         } else {              // else read it from ssj.jar
            String pseudo = "umontreal/ssj/hups/data/";
            StringBuffer pathname = new StringBuffer (pseudo);
            for (int ci = 0; ci < filename.length(); ci++) {
               char ch = filename.charAt (ci);
               if (ch == File.separatorChar)
                  pathname.append ('/');
               else
                  pathname.append (ch);
            }
            InputStream dataInput =
                DigitalNetFromFile.class.getClassLoader().getResourceAsStream (
                  pathname.toString());
            if (dataInput == null)
               throw new FileNotFoundException();
            input = new BufferedReader (new InputStreamReader (dataInput));
         }
         return input;

       } catch (FileNotFoundException e) {
         System.err.println (e + " *** cannot find  " + filename);
         throw e;

      } catch (IOException e) {
         // This will never catch FileNotFoundException since there
         // is a catch clause above.
         System.err.println (e + " cannot read from  " + filename);
         throw e;
      }
   }

   /**
    * Constructs a digital net after reading its parameters from file
    * `filename`. If a file named `filename` can be found relative to the
    * program’s directory, then the parameters will be read from this
    * file; otherwise, they will be read from the file named `filename` in
    * the `ssj.jar` archive. If `filename` is a URL string, it will be
    * read on the World Wide Web. For example, to construct a digital net
    * from the parameters in file `B3S13R9C9St6` in the current directory,
    * one must give the string `"B3S13R9C9St6"` as argument to the
    * constructor. As an example of a file read from the WWW, one may give
    * as argument to the constructor the string `
    * "http://simul.iro.umontreal.ca/ssj/data/Edel/OOA3/B3S13R6C6St4"`.
    * Parameter `w` gives the number of digits of resolution, `r1` is the
    * number of rows, and `s1` is the dimension. Restrictions: `s1` must
    * be less than the maximal dimension, and `r1` less than the maximal
    * number of rows in the data file. Also `w` @f$\ge@f$ `r1`.
    *  @param filename     Name of the file to be read
    *  @param r1           Number of rows for the generating matrices
    *  @param w            Number of digits of resolution
    *  @param s1           Number of dimensions
    */
   public DigitalNetFromFile (String filename, int r1, int w, int s1)
          throws MalformedURLException, IOException
   {
      super ();
      BufferedReader input = null;
      StreamTokenizer st = null;
      try {
         if (filename.startsWith("http:") || filename.startsWith("ftp:"))
            input = openURL(filename);
         else
            input = openFile(filename);
         st = new StreamTokenizer (input);
         readData (st);

      } catch (MalformedURLException e) {
         System.err.println ("   Invalid URL address:   " + filename);
         throw e;
      } catch (FileNotFoundException e) {
         System.err.println ("   Cannot find  " + filename);
         throw e;
      } catch (NumberFormatException e) {
         System.err.println ("   Cannot read number from " + filename);
         throw e;
      }  catch (IOException e) {
         System.err.println ("   IOException:   " + filename);
         throw e;
      }

      if (b == 2) {
         System.err.println ("   base = 2, use DigitalNetBase2FromFile");
         throw new IllegalArgumentException
             ("base = 2, use DigitalNetBase2FromFile");
      }
      if ((double)numCols * Math.log ((double)b) > (31.0 * Math.log (2.0)))
         throw new IllegalArgumentException
            ("DigitalNetFromFile:   too many points" + PrintfFormat.NEWLINE);
      if (r1 > numRows)
         throw new IllegalArgumentException
            ("DigitalNetFromFile:   One must have   r1 <= Max num rows" +
                PrintfFormat.NEWLINE);
      if (s1 > dim)
         throw new IllegalArgumentException
            ("DigitalNetFromFile:   One must have   s1 <= Max dimension" +
                 PrintfFormat.NEWLINE);
      if (w < 0) {
         r1 = w = numRows;
         s1 = dim;
      }
      if (w < numRows)
         throw new IllegalArgumentException
            ("DigitalNetFromFile:   One must have   w >= numRows" +
              PrintfFormat.NEWLINE);

      try {
         readMatrices (st, r1, numCols, s1);
      } catch (NumberFormatException e) {
         System.err.println (e + "   cannot read matrices from " + filename);
         throw e;
      }  catch (IOException e) {
         System.err.println (e + "   cannot read matrices from  " + filename);
         throw e;
      }
      input.close();

      this.filename = filename;
      numRows = r1;
      dim = s1;
      outDigits = w;
      int x = b;
      for (int i=1; i<numCols; i++) x *= b;
      if (x != numPoints) {
         System.out.println ("DigitalNetFromFile:   numPoints != b^k");
         throw new IllegalArgumentException (" numPoints != b^k");
      }

      // Compute the normalization factors.
      normFactor = 1.0 / Math.pow ((double) b, (double) outDigits);
      double invb = 1.0 / b;
      factor = new double[outDigits];
      factor[0] = invb;
      for (int j = 1; j < outDigits; j++)
         factor[j] = factor[j-1] * invb;
  }

   /**
    * Same as  {@link #DigitalNetFromFile() DigitalNetFromFile(filename,
    * r, r, s)} where `s` is the dimension and `r` is given in data file
    * `filename`.
    *  @param filename     Name of the file to be read
    *  @param s            Number of dimensions
    */
   public DigitalNetFromFile (String filename, int s)
          throws MalformedURLException, IOException
   {
       this (filename, -1, -1, s);
   }

   DigitalNetFromFile ()
   {
       super ();
   }


   public String toString() {
      StringBuffer sb = new StringBuffer ("File:   " + filename +
         PrintfFormat.NEWLINE);
      sb.append (super.toString());
      return sb.toString();
   }

/**
 * Writes the parameters and the generating matrices of this digital net to a
 * string. This is useful to check that the file parameters have been read
 * correctly.
 */
public String toStringDetailed() {
      StringBuffer sb = new StringBuffer (toString());
      sb.append (PrintfFormat.NEWLINE + "n = " + numPoints  +
                 PrintfFormat.NEWLINE);
      sb.append ("dim = " + dim  + PrintfFormat.NEWLINE);
      for (int i = 0; i < dim; i++) {
         sb.append (PrintfFormat.NEWLINE + " // dim = " + (1 + i) +
                    PrintfFormat.NEWLINE);
         for (int c = 0; c < numCols; c++) {
            for (int r = 0; r < numRows; r++)
                sb.append (genMat[i*numCols + c][r] + " ");
            sb.append (PrintfFormat.NEWLINE);
         }
      }
      return sb.toString ();
   }


   static class NetComparator implements Comparator {
      // Used to sort list of nets. Sort first by base, then by dimension,
      // then by the number of rows. Don't forget that base = 4 are in files
      // named B4_2* and that the computations are done in base 2.
      public int compare (Object o1, Object o2) {
         DigitalNetFromFile net1 = (DigitalNetFromFile) o1;
         DigitalNetFromFile net2 = (DigitalNetFromFile) o2;
         if (net1.b < net2.b)
            return -1;
         if (net1.b > net2.b)
            return 1;
         if (net1.filename.indexOf("_") >= 0 &&
             net2.filename.indexOf("_") < 0 )
            return 1;
         if (net2.filename.indexOf("_") >= 0 &&
             net1.filename.indexOf("_") < 0 )
            return -1;
         if (net1.dim < net2.dim)
            return -1;
         if (net1.dim > net2.dim)
            return 1;
         if (net1.numRows < net2.numRows)
            return -1;
         if (net1.numRows > net2.numRows)
            return 1;
         return 0;
      }
   }


   private static List getListDir (String dirname) throws IOException {
      try {
         String pseudo = "umontreal/ssj/hups/data/";
         StringBuffer pathname = new StringBuffer (pseudo);
         for (int ci = 0; ci < dirname.length(); ci++) {
            char ch = dirname.charAt (ci);
            if (ch == File.separatorChar)
               pathname.append ('/');
            else
               pathname.append (ch);
         }
         URL url = DigitalNetFromFile.class.getClassLoader().getResource (
                      pathname.toString());
         File dir = new File (url.getPath());
         if (!dir.isDirectory())
            throw new IllegalArgumentException (
               dirname + " is not a directory");
         File[] files = dir.listFiles();
         List alist = new ArrayList (200);
         if (!dirname.endsWith (File.separator))
            dirname += File.separator;
         for (int i = 0; i < files.length; i++) {
            if (files[i].isDirectory())
               continue;
            if (files[i].getName().endsWith ("gz") ||
                files[i].getName().endsWith ("zip"))
               continue;
            DigitalNetFromFile net = new DigitalNetFromFile();
            BufferedReader input = net.openFile(dirname + files[i].getName());
            StreamTokenizer st = new StreamTokenizer (input);
            net.readData (st);
            net.filename = files[i].getName();
            alist.add (net);
         }
         if (alist != null && !files[0].isDirectory())
            Collections.sort (alist, new NetComparator ());
         return alist;

      } catch (NullPointerException e) {
         System.err.println ("getListDir: cannot find directory   " + dirname);
         throw e;

      } catch (NumberFormatException e) {
         System.err.println (e + "***   cannot read number ");
         throw e;

      }  catch (IOException e) {
         System.err.println (e);
         throw e;
      }
   }


   private static String listFiles (String dirname) {
      try {
         String pseudo = "umontreal/ssj/hups/data/";
         StringBuffer pathname = new StringBuffer (pseudo);
         for (int ci = 0; ci < dirname.length(); ci++) {
            char ch = dirname.charAt (ci);
            if (ch == File.separatorChar)
               pathname.append ('/');
            else
               pathname.append (ch);
         }
         URL url = DigitalNetFromFile.class.getClassLoader().getResource (
                      pathname.toString());
         File dir = new File (url.getPath());
         File[] list = dir.listFiles();
         List alist = new ArrayList (200);
         final int NPRI = 3;
         StringBuffer sb = new StringBuffer(1000);
         for (int i = 0; i < list.length; i++) {
            if (list[i].isDirectory()) {
               sb.append (PrintfFormat.s(-2, list[i].getName()));
               sb.append (File.separator + PrintfFormat.NEWLINE);
            } else {
               sb.append (PrintfFormat.s(-25, list[i].getName()));
               if (i % NPRI == 2)
                  sb.append (PrintfFormat.NEWLINE);
            }
         }
         if (list.length % NPRI > 0)
            sb.append (PrintfFormat.NEWLINE);
         return sb.toString();

      } catch (NullPointerException e) {
         System.err.println ("listFiles: cannot find directory   " + dirname);
         throw e;
      }
   }

/**
 * Lists all files (or directories) in directory `dirname`. Only relative
 * pathnames should be used. The files are parameter files used in defining
 * digital nets. For example, calling `listDir("")` will give the list of the
 * main data directory in SSJ, while calling `listDir("Edel/OOA2")` will give
 * the list of all files in directory `Edel/OOA2`.
 */
public static String listDir (String dirname) throws IOException {
      try {
         List list = getListDir (dirname);
         if (list == null || list.size() == 0)
            return listFiles (dirname);
         StringBuffer sb = new StringBuffer(1000);

         sb.append ("Directory:   " + dirname  + PrintfFormat.NEWLINE +
                    PrintfFormat.NEWLINE);
         sb.append (PrintfFormat.s(-25, "     File") +
                    PrintfFormat.s(-15, "       Base") +
                    PrintfFormat.s(-10, "Dimension") +
                    PrintfFormat.s(-10, " numRows") +
                    PrintfFormat.s(-10, "numColumns" +
                    PrintfFormat.NEWLINE));
         int base = 0;
         for (int i = 0; i < list.size(); i++) {
            DigitalNet net = (DigitalNet) list.get(i);
            int j = ((DigitalNetFromFile)net).filename.lastIndexOf
                (File.separator);
            if (net.b != base) {
               sb.append (
      "----------------------------------------------------------------------"
            + PrintfFormat.NEWLINE);
            base = net.b;
            }
            String name = ((DigitalNetFromFile)net).filename.substring(j+1);
            sb.append (PrintfFormat.s(-25, name) +
                       PrintfFormat.d(10, net.b) +
                       PrintfFormat.d(10, net.dim) +
                       PrintfFormat.d(10, net.numRows) +
                       PrintfFormat.d(10, net.numCols) +
                       PrintfFormat.NEWLINE);
         }
         return sb.toString();

      } catch (NullPointerException e) {
         System.err.println (
            "formatPlain: cannot find directory   " + dirname);
         throw e;
      }
   }

   /**
    * Creates a list of all data files in directory `dirname` and writes
    * that list in format HTML in output file `filename`. Each data file
    * contains the parameters required to build a digital net. The
    * resulting list contains a line for each data file giving the name of
    * the file, the base, the dimension, the number of rows and the number
    * of columns of the corresponding digital net.
    */
   public static void listDirHTML (String dirname, String filename)
          throws IOException {
      String list = listDir(dirname);
      StreamTokenizer st = new StreamTokenizer (new StringReader(list));
      st.eolIsSignificant(true);
      st.ordinaryChar('/');
      st.ordinaryChar('_');
      st.ordinaryChar('-');
      st.wordChars('-', '-');
      st.wordChars('_', '_');
      st.slashSlashComments(false);
      st.slashStarComments(false);
      PrintWriter out = new PrintWriter (
                            new BufferedWriter (
                               new FileWriter (filename)));
      out.println ("<html>" + PrintfFormat.NEWLINE +
          "<head>" + PrintfFormat.NEWLINE + "<title>");
      while (st.nextToken () != st.TT_EOL)
         ;
      out.println ( PrintfFormat.NEWLINE + "</title>" +
           PrintfFormat.NEWLINE + "</head>");
//      out.println ("<body background bgcolor=#e1eae8 vlink=#c00000>");
      out.println ("<body>");
      out.println ("<table border>");
      out.println ("<caption> Directory: " + dirname + "</caption>");
      st.nextToken(); st.nextToken();
      while (st.sval.compareTo ("File") != 0)
         st.nextToken();
      out.print ("<tr align=center><th>" + st.sval + "</th>");
      while (st.nextToken () != st.TT_EOL) {
         out.print ("<th>" + st.sval + "</th>" );
      }
      out.println ("</tr>" + PrintfFormat.NEWLINE);
      while (st.nextToken () != st.TT_EOF) {
          switch(st.ttype) {
          case StreamTokenizer.TT_EOL:
             out.println ("</tr>");
             break;
          case StreamTokenizer.TT_NUMBER:
             out.print ("<td>" + (int) st.nval + "</td>" );
             break;
          case StreamTokenizer.TT_WORD:
             if (st.sval.indexOf ("---") >= 0) {
                st.nextToken ();
                continue;
             }
             out.print ("<tr align=center><td>" + st.sval + "</td>");
             break;
          default:
             out.print (st.sval);
             break;
        }
      }

      out.println ("</table>");
      out.println ("</body>" + PrintfFormat.NEWLINE + "</html>");
      out.close();
}

}
