package umontreal.ssj.mcqmctools.anova;

import java.util.ArrayList;

public class ObservationCollectorList<E> extends ArrayList<ObservationCollector<? super E>>
      implements ObservationCollector<E> {

   public ObservationCollectorList() {
      super();
   }
   
   public ObservationCollectorList(int initialCapacity) {
      super(initialCapacity);
   }
   
   @Override public void init() {
      for (ObservationCollector<? super E> collector : this)
         collector.init();
   }
   
   @Override public void observe(E obs) {
      for (ObservationCollector<? super E> collector : this)
         collector.observe(obs);
   }
   
   public String toString() {
      String s = "List of Observation Collectors: ";
      for (ObservationCollector<? super E> collector : this)
         s = s + collector.toString() + ", ";
      return s;
   }
}
