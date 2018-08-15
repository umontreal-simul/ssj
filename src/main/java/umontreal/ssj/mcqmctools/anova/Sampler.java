package umontreal.ssj.mcqmctools.anova;

import umontreal.ssj.mcqmctools.MonteCarloModel;
import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.stat.*;

public interface Sampler {
   /**
    * Returns the number of times the model is simulated each time a sample is produced.
    *
    */
   public int getNumSimulationsPerSample();

   /**
    * Returns the number of samples produced each time the simulateRuns() method is called.
    *
    */
   public int getNumSamples();

   /**
    * Simulates the model multiple times.
    *
    * The samples are added to the observation collector \c collector.
    *
    */
   public <E> void simulateRuns (MonteCarloModel<? extends E> model, ObservationCollector<E> collector);

   /**
    * Simulates the model multiple times.
    *
    * The samples are added to the collector \c collector.
    *
    */
   public void simulateRuns (MonteCarloModelDouble model, Tally collector);
}
