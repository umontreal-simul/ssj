package umontreal.ssj.mcqmctools.anova;

import umontreal.ssj.stat.Tally;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.hups.*;
import umontreal.ssj.mcqmctools.*;

/**
 * QMC sampler.
 *
 * The samples produced by this sampler are deterministic.
 *
 */
public class RQMCSampler implements RandomSampler {

   protected RQMCPointSet points;
   
   public RQMCSampler (RQMCPointSet points, RandomStream stream) {
      this.points = points;
      setStream(stream);
   }

   public RQMCSampler (RQMCPointSet points) {
      this.points = points;
   }
   
   /**
    * Returns the internal RQMC point set.
    *
    */
   public RQMCPointSet getRQMCPointSet() {
      return points;
   }
   
   /** @copydoc RandomSampler::getStream() */
   @Override public RandomStream getStream() {
      return points.getRandomization().getStream();
   }

   /** @copydoc RandomSampler::setStream(RandomStream) */
   @Override public void setStream (RandomStream stream) {
      points.getRandomization().setStream(stream);
   }

   /** @copydoc Sampler::getNumSimulationsPerSample() */
   @Override public int getNumSimulationsPerSample() {
      return 1;
   }

   /** @copydoc Sampler::getNumSamples() */
   @Override public int getNumSamples() {
      return points.getPointSet().getNumPoints();
   }

   /** @copydoc Sampler::simulate(MonteCarloModel<? extends E>, ObservationCollector<E>) */
   @Override public <E> void simulateRuns (MonteCarloModel<? extends E> model, ObservationCollector<E> collector) {
      randomize();
      for (PointSetIterator it = points.iterator(); it.hasNextPoint(); it.resetToNextPoint()) {
         model.simulate(it);
         collector.observe(model.getPerformance());
      }
   }

   /** @copydoc Sampler::simulate(MonteCarloModelDouble, ObservationCollector<Double>) */
   @Override public void simulateRuns (MonteCarloModelDouble model, Tally collector) {
      randomize();
      for (PointSetIterator it = points.iterator(); it.hasNextPoint(); it.resetToNextPoint()) {
          model.simulate(it);
          collector.add(model.getPerformance());
       }
   }

   @Override public String toString() {
      return "RQMC Sampler [points=" + points.getPointSet().getNumPoints() + "]"
         + " [dimension=" + points.getPointSet().getDimension() + "]"
         + " [pointset=" + points.getPointSet().getClass().getSimpleName() + "]"
         + " [randomization=" + points.getRandomization().getClass().getSimpleName() + "]";
   }

   /**
    * Randomizes the integrator by randomizing the internal point set.
    * If the randomization stream is a PointSetIterator, advance to the next point in the stream.
    *
    */
   protected void randomize() {
      points.randomize();
      if (getStream() instanceof PointSetIterator)
         ((PointSetIterator)getStream()).resetToNextPoint();
   }
}
