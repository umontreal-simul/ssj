package ift6561examples;

import java.io.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.hups.*;
import umontreal.ssj.mcqmctools.*;

// This is question 4 of homework 2 of ift6561, A-2018 and A-2023.

public class ProductBarrierIs extends ProductBarrier {

   double maxU1; // For IS, U1 is generated in [0, maxU1)
   double maxU2 = 0.999999999999; // 1.0 - 10^{-12}
   // When u2 > maxU2, we set it to maxU2, to avoid NaN

   // Constructor.
   public ProductBarrierIs(double C, double K, double a, double b, double mu1, double sigma1, double mu2,
         double sigma2) {
      super(C, K, a, b, mu1, sigma1, mu2, sigma2);
      maxU1 = dist1.cdf(a / C);
   }

   // Generates and returns X, without IS.
   public void simulate(RandomStream stream) {
      double CW1 = C * dist1.inverseF(maxU1 * stream.nextDouble());
      double minU2 = dist2.cdf(b / CW1); // U2 must be larger than this.
      double u2 = minU2 + (1.0 - minU2) * stream.nextDouble();
      if (u2 > maxU2)
         u2 = maxU2; // To avoid NaN
      double X = CW1 * dist2.inverseF(u2);
      payoff = (X - K) * maxU1 * (1.0 - minU2); // Unbiased IS estimator.
   }

   // Descriptor of model
   public String toString() {
      return "Simplified financial option with product and barriers, with IS";
   }
   
   public String getTag() {
      return "ProductBarrierIS";
   }

   public static void main(String[] args) throws IOException {
      double C = 100, K = 102, a = 100, b = 102;
      double mu1 = 0.01, sigma1 = 0.04;
      double mu2 = 0.01, sigma2 = 0.04;
      int n = 10000;

      RandomStream stream = new MRG32k3a();
      Tally statX = new TallyStore("Option payoffs"); // To store the n observations of X.
      ProductBarrier pb = new ProductBarrier(C, K, a, b, mu1, sigma1, mu2, sigma2);
      ProductBarrier pbis = new ProductBarrierIs(C, K, a, b, mu1, sigma1, mu2, sigma2);
      System.out.println("C = 100, K = 102, a = 100, mu1 = mu2 = 0.01, sigma1 = sigma2 = 0.04 \n");

      // We first compare MC with MC+IS
      System.out.println("\n********  Comparison of MC vs MC+IS  ******** \n\n");

      pbis.b = pb.b = 102;
      System.out.println("b = 102 \n");
      System.out.println(MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX));
      System.out.println(MonteCarloExperiment.simulateRunsDefaultReportStudent(pbis, n, stream, statX));

      pbis.b = pb.b = 110;
      System.out.println("b = 110 \n");
      System.out.println(MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX));
      System.out.println(MonteCarloExperiment.simulateRunsDefaultReportStudent(pbis, n, stream, statX));

      pbis.b = pb.b = 120;
      System.out.println("b = 120 \n");
      System.out.println(MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX));
      System.out.println(MonteCarloExperiment.simulateRunsDefaultReportStudent(pbis, n, stream, statX));

      // Now we try adding RQMC with and without IS, with Sobol+LMS+DS
      System.out.println("\n********  Comparisons with RQMC, for b = 110  ******** \n\n");

      pbis.b = pb.b = 110;
      int m = 20; // Number of RQMC replicates.
      PointSetRandomization rand = new LMScrambleShift(new MRG32k3a());
      DigitalNet p = new SobolSequence(14, 31, 2); // n = 2^{14} points in 2 dim.
      n = p.getNumPoints(); // Number of RQMC points.
      System.out.println(RQMCExperiment.makeComparisonExperimentMCvsRQMC(pb, stream, p, rand, n, m));
      System.out.println(RQMCExperiment.makeComparisonExperimentMCvsRQMC(pbis, stream, p, rand, n, m));

      p = new SobolSequence(18, 31, 2); // n = 2^{14} points in 2 dim.
      n = p.getNumPoints(); // Number of RQMC points.
      System.out.println(RQMCExperiment.makeComparisonExperimentMCvsRQMC(pb, stream, p, rand, n, m));
      System.out.println(RQMCExperiment.makeComparisonExperimentMCvsRQMC(pbis, stream, p, rand, n, m));
   }
}
