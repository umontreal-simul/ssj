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
   int r = -1;
   double[] coordinates;

   public SobolTestFunc(double c, int s) {
      this.c = c;
      this.s = s;
   }

   public SobolTestFunc(double c, int s, int r) {
      this.c = c;
      this.s = s;
      this.r = r;
   }

   public SobolTestFunc(double[] coordinates, int s) {
      assert s == coordinates.length;
      this.s = s;
      this.coordinates = coordinates;
   }

   public SobolTestFunc(double[] coordinates, int s, int r) {
      assert s == coordinates.length;
      this.s = s;
      this.r = r;
      this.coordinates = coordinates;
   }

   // Computes and returns the value of the function at a given point.
   public double getValue(double[] point) {
      assert point.length == s;
      assert c == 0. || coordinates == null;
      double average = 1.0;
      for (int j = 0; j < s; j++) {
         double coordinate;
         if (c != 0.){
            coordinate = c;
         }
         else {
            coordinate = coordinates[j];
         }
         if (r == -1){
            average *= (1 + coordinate * (point[j] - 0.5));
         }
         else if (r == 1){
            average *= (1 + coordinate * (point[j] - 0.5) * Math.abs(point[j] - 0.5));
         }
         else {
            throw new RuntimeException("r value not set correctly.");
         }
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
