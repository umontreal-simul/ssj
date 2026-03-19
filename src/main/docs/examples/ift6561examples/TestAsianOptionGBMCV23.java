package ift6561examples;

import java.io.*;
import umontreal.ssj.stochprocess.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.stat.list.lincv.ListOfTalliesWithCV;
import umontreal.ssj.mcqmctools.*;
import umontreal.ssj.util.Chrono;
import umontreal.ssj.util.PrintfFormat;

public class TestAsianOptionGBMCV23 {

   // Testing experiments for an Asian Option under a GBM process, with and without
   // a control variate.
   // The program first makes an MC experiment to estimate the distribution of the
   // payoff without the CV.
   // It computes a confidence interval on the mean, then makes a histogram of the
   // positive payoffs.
   // It also produces a standalone Latex file with the histogram.
   // After that, it makes another experiment with the CV, and produces a report
   // that compares the variances
   // and gives a 95% confidence interval.
   
   
   
   /**
    * Simulate the model n times independently using n substreams of `stream`, 
    * and collects the n observations of (payoff, CV1, CV2) in statWithCV.
    * Same as in MonteCarloExperiment. 
    */
   public static void simulateRunsCV (AsianOptionGBMCV model, int n, RandomStream stream,
         ListOfTalliesWithCV<TallyStore> statWithCV) {
      statWithCV.init();
      for (int i = 0; i < n; i++) {
         model.simulate(stream);
         statWithCV.add(model.getPerformance(), model.getValuesCV());
         stream.resetNextSubstream();
      }
   }
   
   public static String getResultsWithCV (ListOfTalliesWithCV<TallyStore> statWithCV) {
      TallyStore statX  = statWithCV.get(0);
      TallyStore statC1 = statWithCV.get(1);
      TallyStore statC2 = statWithCV.get(2);

      // 0. Without CV
      double meanMC = statX.average();
      double varMC = statX.variance();

      // 1. With geometric payoff CV1
      double varC1 = statC1.variance();
      double covC1X = statC1.covariance (statX);
      double meanC1 = statC1.average();
      double beta1 = covC1X / varC1;
      double meanWithCV1 = meanMC - beta1 * meanC1; 
      double varWithCV1 = varMC - beta1 * covC1X;

      // 2. With sum CV2
      double varC2 = statC2.variance();
      double covC2X = statC2.covariance (statX);
      double meanC2 = statC2.average();
      double beta2 = covC2X / varC2;
      double meanWithCV2 = meanMC - beta2 * meanC2; 
      double varWithCV2 = varMC - beta2 * covC2X;

      // 3. With both CV1 and CV2
      double covC1C2 = statC2.covariance (statC1);
      double detSigmaC = varC1 * varC2 - covC1C2 *  covC1C2;
      double beta21 = (varC2 * covC1X - covC1C2 * covC2X) / detSigmaC;
      double beta22 = (varC1 * covC2X - covC1C2 * covC1X) / detSigmaC;
      double meanWithBoth = meanMC - beta21 * meanC1 - beta22 * meanC2;
      double varWithBoth = 1.0 - (covC1X * covC1X * varC2 + covC2X * covC2X * varC1
                            - 2.0 * covC1X * covC2X * covC1C2) / detSigmaC;

      // Result formatting
      PrintfFormat str = new PrintfFormat();
      str.append ("Method      Average     Variance      VRF \n");
      str.append ("MC     "+ PrintfFormat.format (10, 3, 4, meanMC)
                         + "  " + PrintfFormat.format (10, 3, 4, varMC) + "  "
                         + PrintfFormat.format (10, 3, 4, 1.0) + "\n");
      str.append ("CV1    " + PrintfFormat.format (10, 3, 4, meanWithCV1)
                          + "  " + PrintfFormat.format (10, 3, 4, varWithCV1) + "  "
                          + PrintfFormat.format (10, 3, 4, varMC / varWithCV1) + "\n");
      str.append ("CV2    " + PrintfFormat.format (10, 3, 4, meanWithCV2)
                          + ",  " + PrintfFormat.format (10, 3, 4, varWithCV2) + ",  "
                          + PrintfFormat.format (10, 3, 4, varMC / varWithCV2) + "\n");
      str.append ("CV1+CV2" + PrintfFormat.format (10, 1, 4, meanWithBoth)
                          + ",  " + PrintfFormat.format (10, 1, 4, varWithBoth) +",  "
                          + PrintfFormat.format (10, 1, 4, varMC / varWithBoth) + "\n");
      double[] varCV = new double[1];
      statWithCV.varianceWithCV(varCV); 
      str.append ("with List" + PrintfFormat.format (10, 1, 4, statWithCV.averageWithCV(0))
       + ",  " + PrintfFormat.format (10, 1, 4, varCV[0]) + ",  "
       + PrintfFormat.format (10, 1, 4, varMC / varCV[0]) + "\n");
      //+ ",  " + PrintfFormat.format (10, 1, 4, statWithCV.covarianceWithCV(0,0)) +",  "
      //+ PrintfFormat.format (10, 1, 4, varMC / statWithCV.covarianceWithCV(0,0)) + "\n");
      str.append("\n");
      
      return str.toString();
   }

   public static void testOneCase (AsianOptionGBMCV asian, int n, RandomStream stream) {
      ListOfTalliesWithCV<TallyStore> statWithCV = ListOfTalliesWithCV.createWithTallyStore(1, 2);
      simulateRunsCV (asian, n, stream, statWithCV);
      System.out.println ("\n" + asian.toString() + "\n");   
      System.out.println ("Strike price K=" + (int)asian.strike); 
      System.out.println (getResultsWithCV (statWithCV));
      
      ListOfTalliesWithCV<TallyStore> statWithCV2 = ListOfTalliesWithCV.createWithTallyStore(1, 2);
      System.out.println (MonteCarloExperiment.simulateRunsDefaultReportCV 
            (asian, n, stream, statWithCV2, 0.95, 5));
   }
             
      
   public static void main(String[] args) throws IOException {
      int d;
      double T1;
      double T;
      // double K;
      double s0 = 100.0;
      double r = 0.08;
      double sigma = 0.2;
      int n = 10000;     // Number of simulation runs.

      RandomStream stream = new LFSR113();
      GeometricBrownianMotion gbm = new GeometricBrownianMotion(s0, r, sigma, stream);
      Chrono timer = new Chrono();
      
      d = 12;  T1 = 1.0 / 12.0;  T = 72.0 / 365.0;  
      testOneCase (new AsianOptionGBMCV (gbm, r, d, T1, T, 90.0), n, stream);
      testOneCase (new AsianOptionGBMCV (gbm, r, d, T1, T, 95.0), n, stream);
      testOneCase (new AsianOptionGBMCV (gbm, r, d, T1, T, 100.0), n, stream);
      testOneCase (new AsianOptionGBMCV (gbm, r, d, T1, T, 105.0), n, stream);


      System.out.println ("Total CPU time:      " + timer.format() + "\n");

   }

}
