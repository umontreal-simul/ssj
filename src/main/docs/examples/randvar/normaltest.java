import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;

public class normaltest
{
   public static void main (String[] args) {
      // Create 3 parallel streams of random numbers
      RandomStream stream1 = new MRG31k3p();
      RandomStream stream2 = new MRG31k3p();
      RandomStream stream3 = new MRG31k3p();

      // Create 3 parallel streams of normal random variates
      RandomVariateGen gen1 = new NormalGen (stream1);
      RandomVariateGen gen2 = new NormalGen (stream2);
      RandomVariateGen gen3 = new NormalGen (stream3);

      final int n = 5;
      genere (gen1, n);
      genere (gen2, n);
      genere (gen3, n);
   }

   private static void genere (RandomVariateGen gen, int n) {
      double u;
      for (int i = 0; i < n; i++) {
         u = gen.nextDouble();
         System.out.printf ("%12.6f%n", u);
      }
      System.out.println ("----------------------");
   }
}
