package charts;
import umontreal.ssj.charts.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;
import java.awt.Color;

public class HistogramTest1
{
   private static double[] getData() {
      NormalGen gen = new NormalGen(new LFSR113());
      final int N = 100000;
      double[] ad = new double[N];
      for (int i = 0; i < N; i++)
         ad[i] = gen.nextDouble();
      return ad;
   }

   public static void main(String[] args) {
      double[] data = getData();

      HistogramChart chart;
      chart = new HistogramChart("Standard Normal", null, null, data);

      // Customizes the data plot
      HistogramSeriesCollection collec = chart.getSeriesCollection();
      collec.setBins(0, 80);
      double[] bounds = { -4, 4, 0, 5000 };
      chart.setManualRange(bounds);

      chart.view(800, 500);
      chart.toLatexFile("HistogramTest1.tex", 12, 8);
   }
}
