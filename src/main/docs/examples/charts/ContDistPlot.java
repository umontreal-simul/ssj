package charts;
import umontreal.ssj.probdist.*;
import umontreal.ssj.charts.*;

public class ContDistPlot
{
   public static void main (String[] args) {
      ContinuousDistribution dist = new NormalDist();
      ContinuousDistChart plot = new ContinuousDistChart(dist, -3.5, 3.5, 1000);
      plot.viewDensity(600, 400);
   }
}
