package latnetbuilder;

import umontreal.ssj.rng.*;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.*;

public class SobolTestFunc {
   double c;
   int s;
   double[] vals;
   
   public SobolTestFunc (double c, int s) {
      this.c = c;
      this.s = s;
      vals = new double[s];
   }

   // Generates the process S.
   public void generatePath (RandomStream stream) {
       for (int j = 0; j < s; j++) {
    	   double x = stream.nextDouble();
    	   if (x < 0 || x > 1) {
    		   throw new RuntimeException("oups");
    	   }
    	   vals[j] = x;
       }
   }

   // Computes and returns the mean of the process.
   public double getMean () {
       double average = 1.0;
       for (int j = 0; j < s; j++) {
    	   average *= (1 + c * (vals[j] - 0.5));
       }
       return average - 1.0;
   }

   // Performs n indep. runs using stream and collects statistics in statValue.
   public void simulateRuns (int n, RandomStream stream, Tally statValue) {
      statValue.init();
      for (int i=0; i<n; i++) {
         generatePath (stream);
         statValue.add (getMean ());
         stream.resetNextSubstream();
      }
   }

   public static void main (String[] args) {
      int s = 12;
      double c = 0.5;
      SobolTestFunc process = new SobolTestFunc (c, s);
      Tally statValue = new Tally ("Stats on value of Sobol test func");

      Chrono timer = new Chrono();
      int n = 100000;
      process.simulateRuns (n, new MRG32k3a(), statValue);
      statValue.setConfidenceIntervalStudent();
      System.out.println (statValue.report (0.95, 3));
      System.out.println ("Total CPU time:      " + timer.format() + "\n");
   }
}
