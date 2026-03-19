package ift6561examples;

import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.stat.TallyStore;
import umontreal.ssj.util.*;

/**
 * Code inspired by the program of Luc Veillette, 2011.
 */
public class AsianCV {
   double   discount;       // exp(-r * zeta[t])
   double   K;              // Strike price
   int      t;              // Number of observations
   double[] muDelta;        // Differences * (r - sigma^2/2)
   double[] sigmaSqrtDelta; // Sqrt(Zeta_{i+1}  - Zeta_i)*sigma
   double[] logS;           // S: the GBM; logS = log(S)
   double[] alogS;          // S: the antithetic GBM; alogS = log(S)
   double   geoEC;          // Expected value of the geometric CV
   double   sumEC;          // Expected value of the sum CV

   public AsianCV (double r, double sigma, double K, double s0, int t, double[] zeta) {
      discount = Math.exp (-r * zeta[t]);
      this.K = K;
      this.t = t;
      double mu = r - 0.5 * sigma * sigma;
      muDelta = new double[t];
      sigmaSqrtDelta = new double[t];
      logS = new double[t+1];
      alogS = new double[t+1];
      for (int j = 0; j < t; j++) {
         double delta = zeta[j+1] - zeta[j];
         muDelta[j] = mu * delta;
         sigmaSqrtDelta[j] = sigma * Math.sqrt (delta);
      }   
      logS[0] = Math.log (s0);
      alogS[0] = logS [0];
      // Computation of CV expected values
      double my  = 0;
      double sy2 = 0;
      sumEC = 0;
      for (int j = 1; j <= t; j++) {
         my += zeta[j];
         sy2 += (zeta[j]-zeta[j-1])*(t- j + 1)*(t - j + 1);
         sumEC += Math.exp (r*zeta[j]);
      }
      my = logS[0] + my*mu/t;
      sy2 *= sigma*sigma/(t*t);
      double d = (-Math.log (K) + my)/Math.sqrt (sy2);
      geoEC = discount*(Math.exp (my + sy2/2)*NormalDist.cdf01 (d + Math.sqrt (sy2))
                        - K*NormalDist.cdf01 (d));
      sumEC *= s0;
   }

   // Calcul le mouvement Brownien des prix S.
   public void generatePath (RandomStream stream) {
      for (int j = 0; j < t; j++) {
         double z = NormalDist.inverseF01 (stream.nextDouble());
         logS[j+1] = logS[j] + muDelta[j] + sigmaSqrtDelta[j]*z;
         alogS[j+1] = alogS[j] + muDelta[j] - sigmaSqrtDelta[j]*z;
      }
   }

   // Computation of payoff with arithmetic mean
   public double getAriPayoff (double[] logPath) {
      double ariAverage = 0.0;
       for (int j = 1; j <= t; j++)
          ariAverage += Math.exp (logPath[j]);
       ariAverage /= t;
       if (ariAverage > K)
          return discount * (ariAverage - K);
       else 
          return 0.0;
   }

   // Computation of payoff with geometric mean
   public double getGeoControl (double[] logPath) {
       double geoAverage = 0.0;
       for (int j = 1; j <= t; j++)
          geoAverage += logPath[j];
       geoAverage = Math.exp (geoAverage/t);
       if (geoAverage > K)
          return discount * (geoAverage - K);
       else 
          return 0.0;          
   }
   
   // Compute the sum of the observations as a CV
   public double getSumControl (double[] logPath) {
      double sum = 0.0;
       for (int j = 1; j <= t; j++)
          sum += Math.exp (logPath[j]);
       return sum;
   }

   // Performs simulations.
   public void simulateRuns (int n, RandomStream stream) {
      TallyStore statX   = new TallyStore (n);
      TallyStore statC1  = new TallyStore (n);
      TallyStore statC2  = new TallyStore (n);
      TallyStore statAX  = new TallyStore (n);
      TallyStore statAC1 = new TallyStore (n);
      TallyStore statAC2 = new TallyStore (n);

      double ariPayoff, geoControl, sumControl;
      for (int i = 0; i < n; i++) {
         generatePath (stream);
         // Without antithetic variables
         ariPayoff = getAriPayoff (logS);
         geoControl = getGeoControl (logS);
         sumControl = getSumControl (logS);
         statX.add (ariPayoff);
         statC1.add (geoControl);
         statC2.add (sumControl);
         // With antithetic variates
         statAX.add ((getAriPayoff (alogS) + ariPayoff)/2);
         statAC1.add ((getGeoControl (alogS) + geoControl)/2);
         statAC2.add ((getSumControl (alogS) + sumControl)/2);
      }
      // 1. Without CV
      double meanPayoff1 = statX.average();
      double varPayoff1 = statX.variance();
      // 2. With geometric CV
      double varC1 = statC1.variance();
      double covC1X = statC1.covariance (statX);
      double meanC1 = statC1.average();
      double beta2 = covC1X / varC1;
      double meanPayoff2 = meanPayoff1 - beta2*(meanC1 - geoEC); 
      double varPayoff2 = varPayoff1 + beta2*beta2*varC1 - 2*beta2*covC1X;
      // 3. With sum CV
      double varC2 = statC2.variance();
      double covC2X = statC2.covariance (statX);
      double meanC2 = statC2.average();
      double beta3 = covC2X / varC2;
      double meanPayoff3 = meanPayoff1 - beta3*(meanC2 - sumEC); 
      double varPayoff3 = varPayoff1 + beta3*beta3*varC2 - 2*beta3*covC2X;
      // 4. With both CVs
      double covC1C2 = statC2.covariance (statC1);
      double detSigmaC = varC1*varC2 - covC1C2*covC1C2;
      double beta41 = (varC2*covC1X - covC1C2*covC2X)/detSigmaC;
      double beta42 = (varC1*covC2X - covC1C2*covC1X)/detSigmaC;
      double meanPayoff4 = meanPayoff1 - beta41*(statC1.average() - geoEC)
         - beta42*(statC2.average() - sumEC);
      double redFactor4 = 1- (covC1X*covC1X*varC2 - 2*covC1X*covC2X*covC1C2
                              + covC2X*covC2X*varC1)
         /varPayoff1/detSigmaC;
      double varPayoff4 = redFactor4*varPayoff1;
      // 5. With both CVs and AV
      double meanPayoffA1 = statAX.average();
      double varPayoffA1 = statAX.variance();
      double varAC1 = statAC1.variance();
      double covAC1X = statAC1.covariance (statAX);
      double meanAC1 = statAC1.average();
      double varAC2 = statAC2.variance();
      double covAC2X = statAC2.covariance (statAX);
      double meanAC2 = statAC2.average();
      double covAC1C2 = statAC2.covariance (statAC1);
      double detSigmaAC = varAC1*varAC2 - covAC1C2*covAC1C2;
      double beta51 = (varAC2*covAC1X - covAC1C2*covAC2X)/detSigmaAC;
      double beta52 = (varAC1*covAC2X - covAC1C2*covAC1X)/detSigmaAC;
      double meanPayoff5 = meanPayoff1 - beta51*(statC1.average() - geoEC)
         - beta52*(statC2.average() - sumEC);
      double redFactor5 = 1 - (covAC1X*covAC1X*varAC2 - 2*covAC1X*covAC2X*covAC1C2
                               + covAC2X*covAC2X*varAC1)
         /varPayoffA1/detSigmaAC;
      double varPayoff5 = redFactor5*varPayoffA1;
      // Result formatting
      System.out.println ("Strike price K=" + (int)K);
      System.out.println ("Method         Average     Variance        Ratio");
      System.out.println("1. No VRT   "+PrintfFormat.format (10, 1, 4, meanPayoff1)
                         + ",  " + PrintfFormat.format (10, 1, 4, varPayoff1) + ",  "
                         + PrintfFormat.format (10, 1, 4, varPayoff1/varPayoff1));
      System.out.println ("2. VRT Geo. " + PrintfFormat.format (10, 1, 4, meanPayoff2)
                          + ",  " + PrintfFormat.format (10, 1, 4, varPayoff2) + ",  "
                          + PrintfFormat.format (10, 1, 4, varPayoff1/varPayoff2));
      System.out.println ("3. VRT Sum. " + PrintfFormat.format (10, 1, 4, meanPayoff3)
                          + ",  " + PrintfFormat.format (10, 1, 4, varPayoff3) + ",  "
                          + PrintfFormat.format (10, 1, 4, varPayoff1/varPayoff3));
      System.out.println ("4. 2+3      " + PrintfFormat.format (10, 1, 4, meanPayoff4)
                          + ",  " + PrintfFormat.format (10, 1, 4, varPayoff4) +",  "
                          + PrintfFormat.format (10, 1, 4, varPayoff1/varPayoff4));
      System.out.println ("5. 2+3+AV   " + PrintfFormat.format (10, 1, 4, meanPayoff5)
                          + ",  "+ PrintfFormat.format (10, 1, 4, varPayoff5) + ",  "
                          + PrintfFormat.format (10, 1, 4, varPayoff1/varPayoff5));
      System.out.println();
   }

   public static void main (String[] args) { 
      int t;
      RandomStream stream = new MRG32k3a();
      double[] K = new double[] { 80, 90, 100, 110 };
      
      // (i) t=10, zeta_j=110+j, j=1..t.
      t = 10;
      double[] zeta1 = new double[t+1];
      zeta1[0] = 0;
      for (int j=1; j<=t; j++)
         zeta1[j] = (110.0 + j) / 365;
      System.out.println ("(i) t=10, zeta_j=(110+j)/365, j=1..t.");      
      for (int k=0; k<4; k++) {
         AsianCV process = new AsianCV (0.08, 0.2, K[k], 100, t, zeta1);      	
         process.simulateRuns (10000, stream);
      }
      
      // (ii) t=10, zeta_j=12j, j=1..t.
      t = 10;
      double[] zeta2 = new double[t+1];
      for (int j=0; j<=t; j++)
         zeta2[j] = 12.0*j/365;
      System.out.println ("(ii) t=10, zeta_j=12j/365, j=1..t.");
      for (int k=0; k<4; k++) {
         AsianCV process = new AsianCV (0.08, 0.2, K[k], 100, t, zeta2);      	
         process.simulateRuns (10000, stream);
      }
            
      // (iii) t=120, zeta_j=j, j=1..t.
      t = 120;
      double[] zeta3 = new double[t+1];
      for (int j=0; j<=t; j++)
         zeta3[j] = j/365.0;
      System.out.println("(iii) t=120, zeta_j=j/365, j=1..t.");
      for (int k=0; k<4; k++) {
         AsianCV process = new AsianCV (0.08, 0.2, K[k], 100, t, zeta3);
         process.simulateRuns (10000, stream);
      }  
   }
}
