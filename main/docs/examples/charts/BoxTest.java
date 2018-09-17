package charts;
import umontreal.ssj.charts.*;
import umontreal.ssj.randvar.*;
import umontreal.ssj.rng.*;
import java.io.*;

public class BoxTest
{
   public static void main (String[] args) throws IOException {
      int count = 1000;
      double[] data1 = new double[count];
      double[] data2 = new double[count];

      RandomStream stream = new LFSR113();
      RandomVariateGen log = new LognormalGen(stream);
      RandomVariateGen poi = new PoissonGen(stream, 5.0);

      for (int i = 0; i < count; i++) {
         data1[i] = log.nextDouble();
         data2[i] = poi.nextDouble();
      }

      BoxChart bc = new BoxChart("Boxplot1", "Series", "Y", data1, data2);
      bc.view(600, 400);
   }
}
