package charts;
import umontreal.ssj.probdist.*;
import umontreal.ssj.charts.*;
import java.io.*;

public class DistIntTest
{
   public static void main(String[] args) throws IOException {
      PoissonDist dist = new PoissonDist(50);
      DiscreteDistIntChart dic = new DiscreteDistIntChart(dist);

      // Export to Latex format
      String output = dic.toLatexProb(12, 8);  // 12cm width, 8cm height
      Writer file = new FileWriter("DistIntTest.tex");
      file.write(output);
      file.close();
   }
}
