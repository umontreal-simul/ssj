package charts;
import umontreal.ssj.charts.XYLineChart;

public class ChartTest1
{
   private static double[][] getPoints1() {
      double[][] points = new double[2][200];
      for (int i = 0; i < points[0].length; i++) {
         double x = i / 25.0;
         points[0][i] = x;
         points[1][i] = Math.sqrt (x);
      }
      return points;
   }

   private static double[][] getPoints2() {
      double[][] points = new double[2][21];
      for (int i = 0; i < points[0].length; i++) {
         double x = -Math.PI + 2 * i * Math.PI / (points[0].length - 1);
         points[0][i] = x;
         points[1][i] = Math.cos (x);
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
      // Get data; data1 has length 2 and contains one array for
      // X-axis values, and one array for Y-axis values.
      double[][] data1 = getPoints1();
      double[][] data2 = getPoints2();
      double[][] data3 = getPoints3();

      // Create a new chart with the previous data series.
      XYLineChart chart = new XYLineChart(null, "X", "Y", data1, data2, data3);
      chart.toLatexFile("ChartTest1.tex", 12, 8);
   }
}
