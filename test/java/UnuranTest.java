import org.junit.Test;
import static org.junit.Assert.*;
import umontreal.ssj.randvar.UnuranContinuous;
import umontreal.ssj.randvar.UnuranDiscreteInt;
import umontreal.ssj.rng.MRG31k3p;

public class UnuranTest {

   @Test
   public void testUnuranContinuous() {
      final double[] expected = new double[]{
          0.6294093864066692,
         -1.2277564567077497,
         -0.3424878256069986,
          0.08287788846040008,
          1.3744521552951496,
         -0.14927140089739152,
         -0.7238407550828191,
          0.7421332531396772,
          1.657117449078283,
         -0.05304194486402798
      };
      UnuranContinuous gen = new UnuranContinuous(new MRG31k3p(), "normal()");
      for (int i = 0; i < expected.length; i++)
         assertEquals("output[" + i + "]", expected[i], gen.nextDouble(), 1e-9);
   }

   @Test
   public void testUnuranDiscreteInt() {
      final int[] expected = new int[]{ 11, 8, 9, 7, 11, 9, 10, 6, 9, 15 };
      UnuranDiscreteInt gen = new UnuranDiscreteInt(new MRG31k3p(), "poisson(10)");
      for (int i = 0; i < expected.length; i++)
         assertEquals("output[" + i + "]", expected[i], gen.nextInt());
   }
}
