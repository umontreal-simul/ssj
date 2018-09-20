package charts;
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;
import umontreal.ssj.charts.*;
import java.util.Arrays;
import java.awt.Color;

public class EmpiricalChartTest
{
   private static double[] getPoints1() {
      RandomVariateGen gen = new UniformGen(new LFSR113());
      final int N = 10;
      double[] data = new double[N];
      for (int i = 0; i < N; i++)
         data[i] = gen.nextDouble();
      Arrays.sort(data);
      return data;
   }

   private static double[] getPoints2() {
      RandomVariateGen gen = new BetaGen(new LFSR113(), 3, 1);
      final int N = 20;
      double[] data = new double[N];
      for (int i = 0; i < N; i++)
         data[i] = gen.nextDouble();
      Arrays.sort(data);
      return data;
   }

   public static void main(String[] args) {
      double[] data1 = getPoints1();
      double[] data2 = getPoints2();

      // Create a new chart with the previous data series.
      EmpiricalChart chart = new EmpiricalChart(null, null, null, data1, data2);

      // Data plots customizing
      EmpiricalSeriesCollection collec = chart.getSeriesCollection();
      collec.setMarksType(0, "square*");
      collec.setColor(0, Color.MAGENTA);

      chart.enableGrid(0.1, 0.1);            // Enables grid
      chart.toLatexFile("EmpiricalChartTest.tex", 12, 8);
   }
}
