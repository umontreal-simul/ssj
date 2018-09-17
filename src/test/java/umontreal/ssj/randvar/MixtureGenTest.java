package umontreal.ssj.randvar;

import static org.junit.Assert.assertTrue;
import org.junit.Test;
import umontreal.ssj.probdist.*;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;

/**
 * Test the class {@link umontreal.ssj.randvar.MixtureGen}.
 */
public class MixtureGenTest {

   @Test
   public void testMixtureMean() {
      // define the distributions
      Distribution[] dists = new Distribution[2];
      dists[0] = new NormalDist(10, 5);
      dists[1] = new ExponentialDist(0.2);
      
      double[] weights = new double[]{0.5, 0.5};
      
      RandomStream rs = new MRG32k3a();
      
      MixtureGen mg = new MixtureGen(rs, dists, weights);
      
      int n = 100000;
      
      // compute mean
      Tally t = new Tally("Mixture stats");
      for (int i = 0; i < n; i++) {
         t.add(mg.nextDouble());
      }
      
      double trueMean = 7.5;
      double eps = 0.1;
      System.out.println("Real mean: " + trueMean + " , average found: " + t.average());
      assertTrue(t.average() > (trueMean - eps) && t.average() < (trueMean + eps));
   }
   
   @Test
   public void testMixtureMean2() {
      // define the distributions
      Distribution[] dists = new Distribution[3];
      dists[0] = new NormalDist(10, 5);
      dists[1] = new ExponentialDist(0.2);
      dists[2] = new NormalDist(20,5);
      
      double[] weights = new double[]{0.3, 0.3, 0.4};
      
      RandomStream rs = new MRG32k3a();
      
      MixtureGen mg = new MixtureGen(rs, dists, weights);
      
      int n = 100000;
      
      // compute mean
      Tally t = new Tally("Mixture stats");
      for (int i = 0; i < n; i++) {
         t.add(mg.nextDouble());
      }
      
      double trueMean = 12.5;
      double eps = 0.1;
      System.out.println("Real mean: " + trueMean + " , average found: " + t.average());
      assertTrue(t.average() > (trueMean - eps) && t.average() < (trueMean + eps));
   }
   
}
