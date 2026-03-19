package ift6561examples;

import java.io.*;
import umontreal.ssj.charts.HistogramChart;
// import umontreal.ssj.charts.EmpiricalChart;
// import umontreal.ssj.probdist.ContinuousDistribution;
import umontreal.ssj.rng.*;
import umontreal.ssj.mcqmctools.*;
import umontreal.ssj.stat.TallyStore;

/**
 * This class is for the specific stochastic activity network in class San13.
 * The goal is to estimate the distribution of the length T of the longest path
 * and make a histogram.
 */

public class San13Dist extends San13 {

   // The constructor reads link length distributions in a file.
   public San13Dist(String fileName) throws IOException {
      super(fileName);
   }

   public String toString() {
      String s = "SAN network with 9 nodes and 13 links, from Elmaghraby (1977)\n"
            + "Estimate distribution of length T of longest path \n";
      return s;
   }

   public static void main(String[] args) throws IOException {
      int n = 100000;
      San13Dist san = new San13Dist("src/main/docs/examples/ift6561examples/san13a.dat");
      TallyStore statT = new TallyStore("TallyStore for SAN13 example");
      System.out.println(MonteCarloExperiment.simulateRunsDefaultReportStudent(san, n, new LFSR113(), statT, 0.95, 5));
      statT.quickSort();
      Writer file;

      /*
       * TallyStore statTaggregated = statT; int groupSize = 100; if (groupSize > 1)
       * statTaggregated = statT.aggregate (groupSize); double [] aggreg =
       * statTaggregated.getArray(); EmpiricalChart cdf = new
       * EmpiricalChart("Empirical cdf of $T$", "Values of $T$", "cdf", aggreg,
       * statTaggregated.numberObs()); double[] bounds2 = { 0, 200, 0, 1.0 };
       * cdf.setManualRange(bounds2); cdf.view(800, 500); //String cdfLatex =
       * cdf.toLatex(12.0, 8.0); //file = new
       * FileWriter("src/main/docs/examples/ift6561examples/san13cdf.tex");
       * //file.write(cdfLatex); //file.close();
       */

      HistogramChart hist = new HistogramChart("Distribution of $T$", "Values of $T$", "Frequency", statT.getArray(),
            n);
      double[] bounds = { 0, 200, 0, 12000 };
      hist.setManualRange(bounds);
      (hist.getSeriesCollection()).setBins(0, 40, 0, 200);
      hist.view(800, 500);
      String histLatex = hist.toLatex(12.0, 8.0);
      file = new FileWriter("src/main/docs/examples/ift6561examples/san13chart.tex");
      file.write(histLatex);
      file.close();

      // Print p-th empirical quantile
      double p = 0.99;
      int index = (int) Math.round(p * n);
      double xip = statT.getArray()[index];
      System.out.printf("%5.3g -th quantile: %9.6g \n", p, xip);
   }
}
