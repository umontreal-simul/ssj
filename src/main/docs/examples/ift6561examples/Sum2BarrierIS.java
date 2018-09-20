package ift6561examples;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.stat.Tally;

// code for example sum-barrier-is in introduction.
public class Sum2BarrierIS {

   double a;   double b;   double K;
   double phiam1;      // Phi(a-1)
   RandomStream stream  = new MRG32k3a();
 
   public Sum2BarrierIS (double a, double b, double K) {
      this.a = a;   this.b = b;   this.K = K;
      phiam1 = NormalDist.cdf01 (a - 1.0);
   }

   // Generate X with Monte Carlo
   public double simulateOneMCRun() {
      double y1 = 1.0 + NormalDist.inverseF01 (stream.nextDouble());
      double y2 = 1.0 + NormalDist.inverseF01 (stream.nextDouble());
      if (y1 < a && (y1 + y2) > b)
         return (y1 + y2 - K);
      else
         return 0;
   }

   // Generate X with Importance Sampling
   public double simulateOneISRun() {
      double y1 = 1.0 + NormalDist.inverseF01 (stream.nextDouble() * phiam1);
      double phi = NormalDist.cdf01 (b - 1.0 - y1);   // Phi (b-1-y1)
      double y2 = 1.0 + NormalDist.inverseF01 (phi + stream.nextDouble() * (1.0 - phi));
      return ((y1 + y2 - K) * phiam1 * (1.0 - phi));
   }

   public static void main (String[] args) {
      int n = 100000;
      Sum2BarrierIS system = new Sum2BarrierIS (0.5, 2.0, 1.0);

      // with MC
      Tally mcTally = new Tally ("MC estimator");
      for (int i=0; i<n; i++) mcTally.add (system.simulateOneMCRun());
      System.out.println (mcTally.reportAndCIStudent (0.95, 4));

      // with IS
      Tally isTally = new Tally ("IS estimator");
      for (int i=0; i<n; i++) isTally.add (system.simulateOneISRun());
      System.out.println (isTally.reportAndCIStudent (0.95, 4));

      System.out.println ("Variance for MC: " + mcTally.variance());
      System.out.println ("Variance for IS: " + isTally.variance());
      System.out.println ("Variance ratio:  " + mcTally.variance() / isTally.variance());
   }
}
