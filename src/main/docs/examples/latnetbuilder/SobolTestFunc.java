package latnetbuilder;

import umontreal.ssj.rng.*;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.*;

interface Integrand {
   public void simulateRuns(int n, RandomStream stream, Tally statValue);
   public void simulateRuns(int n, double[][] points, Tally statValue);
}

public class SobolTestFunc implements Integrand {
   double c;
   int s;

   public SobolTestFunc(double c, int s) {
      this.c = c;
      this.s = s;
   }

   // Computes and returns the value of the function at a given point.
   public double getValue(double[] point) {
      assert point.length == s;
      double average = 1.0;
      for (int j = 0; j < s; j++) {
         average *= (1 + c * (point[j] - 0.5));
      }
      return average - 1.0;
   }

   // Performs n indep. runs using stream and collects statistics in statValue.
   public void simulateRuns(int n, RandomStream stream, Tally statValue) {
      statValue.init();
      double[] point = new double[s];
      for (int i = 0; i < n; i++) {
         stream.nextArrayOfDouble(point, 0, s);
         statValue.add(getValue(point));
         stream.resetNextSubstream();
      }
   }

   // Performs n indep. runs using points and collects statistics in statValue.
   public void simulateRuns(int n, double[][] points, Tally statValue) {
      assert n == points.length;
      statValue.init();
      for (int i = 0; i < n; i++) {
         statValue.add(getValue(points[i]));
      }
   }

   public static void main(String[] args) {
      int s = 12;
      double c = 0.5;
      SobolTestFunc process = new SobolTestFunc(c, s);
      Tally statValue = new Tally("Stats on value of Sobol test func");

      Chrono timer = new Chrono();
      int n = 100000;
      process.simulateRuns(n, new MRG32k3a(), statValue);
      statValue.setConfidenceIntervalStudent();
      System.out.println(statValue.report(0.95, 3));
      System.out.println("Total CPU time:      " + timer.format() + "\n");
   }
}
