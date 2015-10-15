/*
 * Class:        TextDataWriter
 * Description:  Text data writer
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       David Munger 
 * @since        August 2009

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
package umontreal.ssj.util.io;
import java.io.*;
import java.lang.reflect.Array;

/**
 * Text data writer. Writes fields as columns or as rows in a text file.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class TextDataWriter extends CachedDataWriter {

   /**
    * @name Fields
    * @{
    */

   /**
    * Default value for the column separator.
    */
   public final String DEFAULT_COLUMN_SEPARATOR = "\t";

   /**
    * Default value for the header prefix.
    */
   public final String DEFAULT_HEADER_PREFIX = "";


   protected BufferedWriter out;

   protected Format format;
   protected boolean withHeaders;

   protected String columnSeparator = DEFAULT_COLUMN_SEPARATOR;
   protected String headerPrefix = DEFAULT_HEADER_PREFIX;
   
   protected String floatFormatString = null;


   /**
    * Returns the maximum field length.
    *
    */
   protected int getMaxFieldLength() {
      int nRows = 0;      
      for (DataField f : super.getFields()) {
         if (f.isArray())
            nRows = Math.max(nRows, f.getArrayLength());
      }
      return nRows;
   }

   /**
    * Outputs fields as columns.
    *
    */
   protected void outputAsColumns() throws IOException {
      
      if (withHeaders) {
         // output field headers
         out.write(headerPrefix);
         int iAnonymous = 0;
         boolean firstColumn = true;
         for (DataField f : super.getFields()) {
            // separator
            if (!firstColumn)
               out.write(columnSeparator);
            else
               firstColumn = false;

            if (f.getLabel() == null)
               // anonymous field
               out.write(String.format("_data%02d_", ++iAnonymous));
            else
               // named field
               out.write(f.getLabel());
         }
         out.write("\n");
      }

      int nRows = getMaxFieldLength();

      for (int iRow = 0; iRow < nRows; iRow++) {
         boolean firstColumn = true;
         for (DataField f : super.getFields()) {

            // separator
            if (!firstColumn)
               out.write(columnSeparator);
            else
               firstColumn = false;

            // output field data
            if (f.isArray()) {
               // field is an array, output its current entry
               if (iRow < f.getArrayLength())
                  writeFormat(Array.get(f.asObject(), iRow));
            }
            else {
               // field is not an array, output only in first row
               if (iRow == 0)
                  writeFormat(f.asObject());
            }
         }
         out.write("\n");
      }
   }
 
   
   /**
    * Outputs fields as rows.
    *
    */
   protected void outputAsRows() throws IOException {

      int iAnonymous = 0;

      for (DataField f : super.getFields()) {

         // output field header
         if (withHeaders) {
            if (f.getLabel() == null)
               // anonymous field
               out.write(String.format("_data%02d_", ++iAnonymous));
            else
               // named field
               out.write(f.getLabel());            

            out.write(columnSeparator);
         }
         
         // output field data

         if (f.isArray()) {

            int nCols = f.getArrayLength();

            for (int iCol = 0; iCol < nCols; iCol++) {

               // separator
               if (iCol > 0)
                  out.write(columnSeparator);
               
               writeFormat(Array.get(f.asObject(), iCol));
            }
         }
         else {
            writeFormat(f.asObject());
         }
         
         out.write("\n");

      }
   }
   
   /**
    * Formats the object in accordance with the current format strings settings.
    *
    */
   protected void writeFormat(Object o) throws IOException {
      String s = null;
      if (floatFormatString != null && (o instanceof Double || o instanceof Float))
         s = String.format((java.util.Locale)null, floatFormatString, o); // pass null to avoid localization
      else
         s = o.toString();
      out.write(s);
   }

   /**
    * @}
    */

   /**
    * Output format: organize fields as columns or as rows.
    */
   public enum Format { COLUMNS, ROWS }

   /**
    * Class constructor. Truncates any existing file with the specified
    * name.
    *  @param filename     name of the file to write to
    *  @param format       organize fields as columns if set to `COLUMNS`
    *                      or as rows if set to `ROWS`
    *  @param withHeaders  output headers or not
    */
   public TextDataWriter (String filename, Format format, boolean withHeaders)
         throws IOException {
      this.out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename)));
      this.format = format;
      this.withHeaders = withHeaders;
   }

   /**
    * Class constructor. Truncates any conflicting file.
    *  @param file         file to write to
    *  @param format       organize fields as columns if set to `COLUMNS`
    *                      or as rows if set to `ROWS`
    *  @param withHeaders  output headers or not
    */
   public TextDataWriter (File file, Format format, boolean withHeaders)
         throws IOException {
      this.out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file)));
      this.format = format;
      this.withHeaders = withHeaders;
   }

   /**
    * Class constructor.
    *  @param outputStream output stream to write to
    *  @param format       organize fields as columns if set to `COLUMNS`
    *                      or as rows if set to `ROWS`
    *  @param withHeaders  output headers or not
    */
   public TextDataWriter (OutputStream outputStream, Format format,
                          boolean withHeaders)
         throws IOException {
      this.out = new BufferedWriter(new OutputStreamWriter(outputStream));
      this.format = format;
      this.withHeaders = withHeaders;
   }

   /**
    * Changes the output format.
    *  @param format       organize fields as columns if set to `COLUMNS`
    *                      or as rows if set to `ROWS`
    */
   public void setFormat (Format format) {
      this.format = format;
   }

   /**
    * Sets the format string used to output floating point numbers.
    *  @param formatString format string (e.g., <tt>%.4g</tt>)
    */
   public void setFloatFormatString (String formatString) {
      this.floatFormatString = formatString;
   }

   /**
    * Changes the column separator.
    */
   public void setColumnSeparator (String columnSeparator) {
      this.columnSeparator = columnSeparator;
   }

   /**
    * Changes the header prefix (a string that indicates the beginning of
    * the header line for the `COLUMNS` format).
    */
   public void setHeaderPrefix (String headerPrefix) {
      this.headerPrefix = headerPrefix;
   }

   /**
    * Flushes any pending data and closes the file or stream.
    */
   public void close() throws IOException {
      if (format == Format.COLUMNS)
         outputAsColumns();
      else if (format == Format.ROWS)
         outputAsRows();
      out.close();
   }

}