package charts;
import umontreal.ssj.charts.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;
import java.awt.Color;

public class HistogramChartTest
{
   private static double[] getPoints1() {
      NormalGen gen = new NormalGen(new MRG32k3a(), 0, 2);
      final int N = 100000;
      double[] ad = new double[N];
      for (int i = 0; i < N; i++)
         ad[i] = gen.nextDouble();
      return ad;
   }

   private static double[] getPoints2() {
      ExponentialGen gen = new ExponentialGen(new MRG32k3a(), 1);
      final int N = 100000;
      double[] ad = new double[N];
      for (int i = 0; i < N; i++)
         ad[i] = gen.nextDouble();
      return ad;
   }

   public static void main(String[] args) {
      double[] data1 = getPoints1();
      double[] data2 = getPoints2();

      // Create a new chart with the previous data series.
      HistogramChart chart = new HistogramChart(null, null, null, data1, data2);

      // Customizes the data plots
      HistogramSeriesCollection collec = chart.getSeriesCollection();
      collec.setColor(0, new Color(255, 0, 0, 128));
      collec.setColor(1, new Color(0, 255, 0, 128));
      collec.setBins(0, 40, -6, 6);

      // Define range bounds.
      double[] bounds = { -6, 6, 0, 30000 };
      chart.setManualRange00(bounds, true, true);

      chart.toLatexFile("HistogramChartTest.tex", 12, 8);
   }
}
