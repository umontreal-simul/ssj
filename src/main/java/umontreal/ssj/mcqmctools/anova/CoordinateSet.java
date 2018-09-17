package umontreal.ssj.mcqmctools.anova;

import umontreal.ssj.stat.Tally;

import java.util.*;

/**
 * Represents a set of coordinates.
 *
 * Coordinate indices start at 0.
 *
 */
public abstract class CoordinateSet {

   public abstract List<Integer> asList();

   /**
    * Returns the maximum coordinate index, starting at 0 for the first coordinate.
    *
    */
   public abstract int maxCoordinate();

   /**
    * Returns all subsets of the current coordinate set, whose cardinality is at most \c maxOrder.
    *
    */
   public abstract List<CoordinateSet> subsets(boolean includeEmptySet, int maxOrder);

   @Override public boolean equals(Object o) {
      if (o instanceof CoordinateSet)
         return ((CoordinateSet)o).containsAll(this) && this.containsAll((CoordinateSet)o);
      else
         return false;
   }


   /**
    * Returns \c true if the current set contains coordinate \c coord.
    * Subclasses should override this method.
    *
    */
   public boolean contains(int coord) {
      return asList().contains(coord);
   }
   
   /**
    * Returns \c true if the current set contains all coordinates in \c cs.
    * Subclasses should override this method. The default implementation is inefficient.
    *
    */
   public boolean containsAll(CoordinateSet cs) {
      return asList().containsAll(cs.asList());
   }

   /**
    * Returns the cardinality of the current coordinate set.
    * Subclasses should override this method. The default implementation is inefficient.
    *
    */
   public int cardinality() {
      return asList().size();
   }
   
   /**
    * Returns \c true if \c cs is a subset of the current coordinate set.
    *
    */
   public boolean isSubset(CoordinateSet cs) {
      return cs.containsAll(this);
   }

   /**
    * Returns all subsets of the current coordinate set, whose cardinality is at most \c maxOrder.
    * Includes the empty sets.
    *
    */
   public List<CoordinateSet> subsets() {
      return subsets(true, Integer.MAX_VALUE);
   }

   /**
    * Returns all subsets of the current coordinate set, whose cardinality is at most \c maxOrder.
    * Does not include the empty set.
    *
    */
   public List<CoordinateSet> subsetsNotEmpty() {
      return subsets(false, Integer.MAX_VALUE);
   }

   /**
    * Returns all subsets of the current coordinate set, whose cardinality is at most \c maxOrder.
    * Includes the empty set.
    *
    */
   public List<CoordinateSet> subsets(int maxOrder) {
      return subsets(true, maxOrder);
   }

   /**
    * Returns all subsets of the current coordinate set, whose cardinality is at most \c maxOrder.
    * Does not include the empty set.
    *
    */
   public List<CoordinateSet> subsetsNotEmpty(int maxOrder) {
      return subsets(false, maxOrder);
   }

   @Override public String toString() {
      String s = "";
      for (Integer coord : asList()) {
         if (s.length() > 0) s += ",";
         s += (coord + 1);
      }
      return "{" + s + "}";
   }
}
