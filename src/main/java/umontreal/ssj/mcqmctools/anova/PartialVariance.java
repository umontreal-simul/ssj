package umontreal.ssj.mcqmctools.anova;

import umontreal.ssj.stat.Tally;


/**
 * Represents the partial variance of a function with respect to a given coordinate set.
 *
 */
public class PartialVariance extends Tally implements Comparable<PartialVariance> {

   protected CoordinateSet coords;
   protected Tally totalVar;

   protected PartialVariance(CoordinateSet coords) {
      this(coords, null);
   }
   
   protected PartialVariance(CoordinateSet coords, Tally totalVar) {
      super("variance for coordinates " + coords);
      this.coords = coords;
      this.totalVar = totalVar;
   }
   
   /**
    * Returns the coordinate set associated with the current partial variance.
    *
    */
   public CoordinateSet getCoordinates() {
      return coords;
   }
   
   /**
    * Returns the sensitivity index (fraction of the total variance) of the coordinate set under
    * consideration.
    *
    * Throws an IllegalStateException if the total variance is not set.
    * Should be overridden in deriving classes. 
    */
   public double sensitivityIndex() {

      if (totalVar == null)
         throw new IllegalStateException("trying to access the sensitivity index without"
              + " a reference to the total variance");

      return average() / totalVar.average();
   }
   
   /**
    * Returns 1, 0 or -1 if the current partial variance is larger, equal or smaller than the
    * partial variance \c var, respectively.
    *
    */
   @Override public int compareTo(PartialVariance var) {
      double v = average();
      double vx = var.average();
      return v > vx ? 1 : v < vx ? -1 : 0;
   }
   
   @Override public String toString() {
      
      String s = String.format("%30s: %9.4g", coords.toString(), average());
      
      if (numberObs() > 1) {
         double dvar = Math.sqrt(variance() / numberObs());
         s += String.format(" ± %.2g", dvar);
      }
         
      if (totalVar != null) {
         double varFrac = sensitivityIndex();
         if (varFrac >= 0)
            s += String.format("  (%.4g %%)", 100 * varFrac);
      }

      return s;
   }
}
