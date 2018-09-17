package charts;
import umontreal.ssj.charts.XYLineChart;

public class NormalChart
{
   private static double[][] getPoints() {
      // The density of the standard normal probability distribution
      // points contains one array for X values and one array for Y values
      final int N = 400;
      double[][] points = new double[2][N + 1];
      final double CPI = Math.sqrt (2*Math.PI);
      for (int i = 0; i <= N; ++i) {
         double x = -3.5 + i * 7.0 / N;
         points[0][i] = x;
         points[1][i] = Math.exp (-x*x/2.0) / CPI;
      }
      return points;
   }

   public static void main(String[] args) {
      double[][] points = getPoints();
      XYLineChart chart = new XYLineChart(null, "X", null, points);
      chart.setAutoRange00(true, true);      // Axes pass through (0,0)
      chart.toLatexFile("NormalChart.tex", 12, 8);
      chart.view(800,500);
   }
}
