package umontreal.ssj.mcqmctools;

public interface ObservationCollector<E> {

   /**
    * Clears all observations.
    *
    */
   public void init();

   /**
    * Adds an observation.
    *
    */
   public void observe(E obs);
}
