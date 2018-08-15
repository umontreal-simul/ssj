package umontreal.ssj.mcqmctools.anova;

import umontreal.ssj.stat.Tally;
import umontreal.ssj.hups.*;
import umontreal.ssj.mcqmctools.*;

/**
 * QMC sampler.
 *
 * The samples produced by this sampler are deterministic.
 *
 */
public class QMCSampler implements Sampler {

   protected PointSet points;
   
   public QMCSampler (PointSet points) {
      this.points = points;
   }
   
   public PointSet getPointSet() {
      return points;
   }
   
   /** @copydoc Sampler::getNumSimulationsPerSample() */
   @Override public int getNumSimulationsPerSample() {
      return 1;
   }

   /** @copydoc Sampler::getNumSamples() */
   @Override public int getNumSamples() {
      return points.getNumPoints();
   }

   /** @copydoc Sampler::simulate(MonteCarloModel<? extends E>, ObservationCollector<E>) */
   @Override public <E> void simulateRuns (MonteCarloModel<? extends E> model, ObservationCollector<E> collector) {

      for (PointSetIterator it = points.iterator(); it.hasNextPoint(); it.resetToNextPoint()) {
    	  model.simulate(it);
          collector.observe(model.getPerformance());
      }
   }

   /** @copydoc Sampler::simulate(MonteCarloModelDouble, ObservationCollector<Double>) */
   @Override public void simulateRuns (MonteCarloModelDouble model, Tally collector) {
      for (PointSetIterator it = points.iterator(); it.hasNextPoint(); it.resetToNextPoint()) {
          model.simulate(it);
          collector.add(model.getPerformance());
       }
   }

   @Override public String toString() {
      return "QMC Sampler [nPoints=" + points.getNumPoints() + "]"
         + " [points=" + points.getClass().getSimpleName() + "]";
   }
}
