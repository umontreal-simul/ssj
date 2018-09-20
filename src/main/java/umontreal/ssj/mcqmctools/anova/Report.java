package umontreal.ssj.mcqmctools.anova;

import java.util.*;

public class Report {
   
   protected static int labelWidth = 24;
   protected static int indentation = 3;
   protected static int floatPrecision = 4;
   protected static int valueWidth = 0;
   
   protected String label;
   protected String value;
   protected List<Report> subReports;

   public static final String NEWLINE = System.getProperty("line.separator");

   // static setup methods
   
   /**
    * Set to <0 for left-aligned label with given width.
    * Set to 0 for label-adjusted width.
    * Set to >0 for right-aligned value with given width.
    *
    */
   public static void setLabelWidth(int labelWidth) {
      Report.labelWidth = labelWidth;
   }

   public static void setIndentation(int indentation) {
      Report.indentation = indentation;
   }

   public static void setFloatPrecision(int floatPrecision) {
      Report.floatPrecision = floatPrecision;
   }
   
   /**
    * Set to <0 for left-aligned value with given width.
    * Set to 0 for value-adjusted width.
    *
    */
   public static void setValueWidth(int valueWidth) {
      Report.valueWidth = valueWidth;
   }

   // constructors
 
   public Report(String label, String value) {
      this.label = label;
      this.value = (value == null) ? null : String.format(
            "%" + (valueWidth == 0 ? "" : valueWidth) + "s", value);
      subReports = new ArrayList<Report>();
   }
   
   public Report(String label, int value) {
      this(label, String.format("%d", value));
   }
   
   public Report(String label, double value) {
      this(label, String.format("%." + floatPrecision + "g", value));
   }
   
   public Report(String label, int[] value) {
      this(label, formatVector(value));
   }
   
   public Report(String label, double[] value) {
      this(label, formatVector(value));
   }
   
   public Report(String label, boolean value) {
      this(label, value ? "yes" : "no");
   }
   
   public Report(String label) {
      this.label = label;
      this.value = null;
      subReports = new ArrayList<Report>();
   }
   
   public Report() {
      this(null);
   }
   
   
   // add a sub-report
   
   public void add(Report subReport) {
      subReports.add(subReport);
   }

   // shortcut methods
   
   public void add(String label, String value) {
      add(new Report(label, value));
   }
   
   public void add(String label, int value) {
      add(new Report(label, value));
   }
   
   public void add(String label, double value) {
      add(new Report(label, value));
   }
   
   public void add(String label, boolean value) {
      add(new Report(label, value));
   }

   public void add(String label, int[] value) {
      add(new Report(label, value));
   }
   
   public void add(String label, double[] value) {
      add(new Report(label, value));
   }
   
   public void add(String label) {
      add(new Report(label));
   }
   

   // stringification methods
   
   public String toString() {
      return toString(0);
   }
   
   public String toString(int baseIndentation) {
      String s = "";
      if (label != null) {
         if (baseIndentation > 0)
            s += String.format("%" + baseIndentation + "s", "");
         if (value != null)
            s += String.format("%" + (labelWidth == 0 ? "" : labelWidth) + "s %s",
                  label + ":", value);
         else
            s += "==> " + label;
         s += NEWLINE;
      }
      if (!subReports.isEmpty()) {
         for (Report report : subReports)
            s += report.toString(baseIndentation + indentation);
      }
      return s;
   }

   protected static String formatVector(int[] value) {
      String s = "";
      for (int i = 0; i < value.length; i++) {
         if (i > 0) s += " ";
         s += String.format("%d", value[i]);
      }
      return s;
   }

   protected static String formatVector(double[] value) {
      String s = "";
      for (int i = 0; i < value.length; i++) {
         if (i > 0) s += " ";
         s += String.format("%." + floatPrecision + "g", value[i]);
      }
      return s;
   }
}
