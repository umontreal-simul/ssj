package umontreal.ssj.mcqmctools;

import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.probdist.*;
import umontreal.ssj.gof.*;

import java.util.*;

/**
 * IMPORTANT: The observations passed to observe() must implement
 * equals() and hashCode() properly.
 *
 */
public class ObservationCounter<E> extends BasicObservationCollector<E> {

   protected Map<E, Integer> countsMap;
   protected int nCollisions;
   
   public ObservationCounter() {
      this.countsMap = new HashMap<E, Integer>();
      this.nCollisions = 0;
   }

   public void init() {
      super.init();
      countsMap.clear();
      nCollisions = 0;
   }

   /**
    * Do not modify the returned map!
    *
    */
   public Map<E, Integer> getCounts() {
      return countsMap;
   }
   
   public int getNumCollisions() {
      return nCollisions;
   }

   public void observe(E observation) {
      super.observe(observation);
      if (countsMap.containsKey(observation)) {
         countsMap.put(observation, countsMap.get(observation) + 1);
         nCollisions++;
      }
      else {
         countsMap.put(observation, new Integer(1));
      }
   }
   
   public Report report() {
      Report r = super.report();
      r.add("distinct observations", countsMap.size());
      r.add("average count", (double)getNumObservations() / countsMap.size());
      r.add("number of collisions", getNumCollisions());
      
      int minCount = -1;
      int maxCount = 0;
      for (Integer count : countsMap.values()) {
         if (minCount == -1) minCount = count;
         minCount = Math.min(minCount, count);
         maxCount = Math.max(maxCount, count);
      }

      r.add("minimum count", minCount);
      r.add("maximum count", maxCount);

      return r;
   }
      
   
   /**
    * Returns the p-value of the chi2 test.
    *
    * @param expectedProbs expected probabilities.
    *
    */
   public double testChi2(Map<? extends E, Double> expectedProbs) {
      int[] countsArray = new int[countsMap.size()];
      double[] expectedProbsArray = new double[countsMap.size()];
      int i = 0;
      for (E observation : countsMap.keySet()) {
         countsArray[i] = countsMap.get(observation);
         expectedProbsArray[i] = expectedProbs.get(observation);
         i++;
      }
      
      return getPValue(expectedProbsArray, countsArray);
   }

   public double testChi2(double expectedProb) {
      int[] countsArray = new int[countsMap.size()];
      double[] expectedProbsArray = new double[countsMap.size()];
      Arrays.fill(expectedProbsArray, expectedProb);
      int i = 0;
      for (Integer count : countsMap.values())
         countsArray[i++] = count;
      
      return getPValue(expectedProbsArray, countsArray);
   }
   
   /**
    * Returns the p-value of the chi2 test.
    *
    */
   protected double getPValue(double[] expectedProbs, int[] counts) {
      double y = GofStat.chi2(expectedProbs, counts, 0, counts.length-1);
      return ChiSquareDist.barF(counts.length, 5, y);
   }
   
   public String toString() {
      return "Observation Counter";
   }   
}
