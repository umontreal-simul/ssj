package tutorial;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.Chrono;

// Class to simulate and compare two different (S,s) policies with CRNs.
public class InventoryCRN extends Inventory {

   Tally statDiff = new Tally ("stats on difference");

   public InventoryCRN (double lambda, double c, double h,
                        double K, double k, double p) {
      super (lambda, c, h, K, k, p);
   }

   public void simulateDiff (int n, int m, int s1, int S1, int s2, int S2) {
      statDiff.init();
      for (int i = 0; i < n; i++) {
         double value1 = simulate (m, s1, S1);
         double value2 = simulate (m, s2, S2);
         statDiff.add (value2 - value1);
      }
   }

   public void simulateDiffCRN (int n, int m, int s1, int S1, int s2, int S2) {
      statDiff.init();
      streamDemand.resetStartStream();
      streamOrder.resetStartStream();
      for (int i = 0; i < n; i++) {
         double value1 = simulate (m, s1, S1);
         streamDemand.resetStartSubstream();
         streamOrder.resetStartSubstream();
         double value2 = simulate (m, s2, S2);
         statDiff.add (value2 - value1);
         streamDemand.resetNextSubstream();
         streamOrder.resetNextSubstream();
      }
   }

   public static void main (String[] args) {
      InventoryCRN system = new InventoryCRN (100.0, 2.0, 0.1, 10.0, 1.0, 0.95);
      Chrono timer = new Chrono();

      system.simulateDiff (500, 2000, 80, 198, 80, 200);
      system.statDiff.setConfidenceIntervalStudent();
      System.out.println (system.statDiff.report (0.9, 3));
      double varianceIndep = system.statDiff.variance();
      System.out.println ("Total CPU time: " + timer.format() + "\n");
      
      timer.init();
      system.simulateDiffCRN (500, 2000, 80, 198, 80, 200);
      System.out.println (system.statDiff.report (0.9, 3));
      double varianceCRN = system.statDiff.variance();
      System.out.println ("Total CPU time: " + timer.format());
      System.out.printf ("Variance ratio:  %8.4g%n", varianceIndep/varianceCRN);
   }
}
