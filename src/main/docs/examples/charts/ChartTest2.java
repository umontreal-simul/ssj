package charts;
import umontreal.ssj.charts.*;
import java.awt.Color;

public class ChartTest2
{
   private static double[][] getPoints1() {
      double[][] points = new double[2][40];
      for (int i = 0; i < points[0].length; i++) {
         double x = i / 4.0;
         points[0][i] = x;
         points[1][i] = Math.sqrt(x);
      }
      return points;
   }

   private static double[][] getPoints2() {
      double[][] points = new double[2][21];
      for (int i = 0; i < points[0].length; i++) {
         double x = -Math.PI + 2 * i * Math.PI / (points[0].length - 1);
         points[0][i] = x;
         points[1][i] = Math.cos(x);
      }
      return points;
   }

   private static double[][] getPoints3() {
      double[][] points = new double[2][11];
      for (int i = 0; i < points[0].length; i++) {
         points[0][i] = -5 + i;
         points[1][i] = -3 + i;
      }
      return points;
   }

   public static void main(String[] args) {
      double[][] data1 = getPoints1();
      double[][] data2 = getPoints2();
      double[][] data3 = getPoints3();

      // Create a new chart with the previous data series.
      XYLineChart chart = new XYLineChart(null, "X", "Y", data1, data2, data3);

      // Customizing axes
      Axis xaxis = chart.getXAxis();
      Axis yaxis = chart.getYAxis();
      String[] labels = { "-9", "$-\\lambda$", "$-\\sqrt{2}$",
                          "0", "$\\frac{14}{\\pi}$", "\\LaTeX" };
      double[] values = { -9, -5, -Math.sqrt(2), 0, 14.0 / Math.PI, 9 };
      xaxis.setLabels(values, labels);
      yaxis.setLabels(1);

      // Data plots customizing
      XYListSeriesCollection collec = chart.getSeriesCollection();
      collec.setColor(0, new Color(0, 64, 128));
      collec.setName(0, "$f(x) = \\sqrt(x)$");
      collec.setMarksType(0, "");
      collec.setDashPattern(0, "dotted");
      collec.setName(1, "$f(x) = \\cos(x)$");
      collec.setMarksType(1, "");
      collec.setColor(2, Color.ORANGE);
      collec.setPlotStyle(2, "ycomb,very thick");
      collec.setMarksType(2, "*");

      // Export to LaTex format
      chart.toLatexFile("ChartTest2.tex", 12, 8);  // 12cm width, 8cm height
   }
}
