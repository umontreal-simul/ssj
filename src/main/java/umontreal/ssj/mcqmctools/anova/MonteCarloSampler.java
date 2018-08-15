package umontreal.ssj.mcqmctools.anova;

import umontreal.ssj.stat.Tally;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.hups.PointSetIterator;
import umontreal.ssj.mcqmctools.*;

/**
 * Monte Carlo sampler.
 *
 * The samples produced by this class are random and independent if the
 * underlying random stream produces random and independent outputs.
 *
 */
public class MonteCarloSampler implements RandomSampler {

   protected int nSamples;
   protected RandomStream stream;
   
   public MonteCarloSampler (int nSamples) {
      this.nSamples = nSamples;
   }
   
   public MonteCarloSampler (int nSamples, RandomStream stream) {
      this.nSamples = nSamples;
      this.stream = stream;
   }
   
   /** @copydoc RandomIntegrator::getStream() */
   public RandomStream getStream() {
      return stream;
   }
   
   /** @copydoc Sampler::getNumSimulationsPerSample() */
   @Override public int getNumSimulationsPerSample() {
      return 1;
   }

   /** @copydoc Sampler::getNumSamples() */
   @Override public int getNumSamples() {
      return nSamples;
   }

   /**
    * Sets the number of samples.
    *
    */
   public void setNumSamples (int nSamples) {
      this.nSamples = nSamples;
   }

   /**
    * Use \c stream to produce random points.
    *
    * Sets the internal random stream to \c stream.
    *
    */
   @Override public void setStream (RandomStream stream) {
      this.stream = stream;
   }

   /** @copydoc Sampler::simulate(MonteCarloModel<? extends E>, ObservationCollector<E>) */
   @Override public <E> void simulateRuns (MonteCarloModel<? extends E> model, ObservationCollector<E> collector) {

      PointSetIterator psit = null;
      if (stream instanceof PointSetIterator)
         psit = (PointSetIterator)stream;
      
      for (int i = 0; i < nSamples; i++) {
         model.simulate(stream);
         collector.observe(model.getPerformance());
         if (psit != null)
            psit.resetToNextPoint();
      }
   }

   /** @copydoc Sampler::simulate(MonteCarloModelDouble, ObservationCollector<Double>) */
   @Override public void simulateRuns (MonteCarloModelDouble model, Tally collector) {
 
      boolean isPointSet = (stream instanceof PointSetIterator);
      
      for (int i = 0; i < nSamples; i++) {
         model.simulate(stream);
         collector.add(model.getPerformance());
         if (isPointSet)
            ((PointSetIterator)stream).resetToNextPoint();
      }
   }
   
   @Override public String toString() {
      String s = "Monte Carlo Sampler [samples=" + getNumSamples() + "]";
      if (getStream() != null)
         s += " [stream=" + getStream().getClass().getSimpleName() + "]";
      return s;
   }   
}
