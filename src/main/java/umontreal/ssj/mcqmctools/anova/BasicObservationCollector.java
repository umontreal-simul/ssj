package umontreal.ssj.mcqmctools.anova;

/**
 * Does nothing but counting the total number of observations.
 * Provides useful functionality for deriving classes.
 *
 * Deriving classes overriding methods, e.g. init() or add(), must call their parent's versions.
 *
 */
public class BasicObservationCollector<E> implements ObservationCollector<E> {

   protected String label;
   protected int nObservations;
   
   public BasicObservationCollector() {
      this(null);
   }

   /**
    * @param label  label describing the collector
    *
    */
   public BasicObservationCollector(String label) {
      this.label = label;
      nObservations = 0;
   }
   
   public String getLabel() {
      return label;
   }
   
   public int getNumObservations() {
      return nObservations;
   }
   
   /**
    * Resets the collector. Forgets all previous observations.
    *
    */
   public void init() {
      nObservations = 0;
   }
   
   /**
    * Updates the number of observations.
    *
    */
   public void observe(E observation) {
      nObservations++;
   }
   
   /**
    * Returns a report of measurements on the data collected by the collector.
    *
    */
   public Report report() {
      String collectorLabel = (label == null) ? this.getClass().getSimpleName() : label;
      Report r = new Report(collectorLabel);
      r.add("number of observations", getNumObservations());
      return r;
   }
   
   public String toString() {
      return "Basic Observation Collector";
   }
}
