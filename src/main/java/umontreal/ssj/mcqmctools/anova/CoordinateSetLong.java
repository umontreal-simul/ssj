package umontreal.ssj.mcqmctools.anova;

import java.util.*;

/**
 * Implementation of CoordinateSet using a \c long bit-mask internal representation.
 *
 */
public class CoordinateSetLong extends CoordinateSet {

   protected long mask;
   
   /**
    * Constructs a coordinate set with corresponding bit mask \c mask.
    *
    */
   public CoordinateSetLong(long mask) {
      this.mask = mask;
   }

   /**
    * Returns the bit-mask representation of the current coordinate set.
    *
    */
   public long getMask() {
      return mask;
   }

   @Override public boolean equals(Object o) {
      if (o instanceof CoordinateSetLong)
         return mask == ((CoordinateSetLong)o).mask;
      else
         return super.equals(o);
   }

   @Override public List<Integer> asList() {
      List<Integer> list = new ArrayList<Integer>();
      long x = mask;
      int coord = 0;
      while (x != 0) {
         if ((x & 1) == 1)
            list.add(coord);
         coord++;
         x = x >> 1;
      }
      return list;
   }
   
   @Override public boolean contains(int coord) {
      return ((mask >> coord) & 1) == 1;
   }

   @Override public boolean containsAll(CoordinateSet cs) {
      if (cs instanceof CoordinateSetLong)
         return (mask | ((CoordinateSetLong)cs).mask) == mask;
      else
         return super.containsAll(cs);
   }

   @Override public int cardinality() {
      // count the bits set to 1
      long x = mask;
      int count = 0;
      while (x != 0) {
         count++;
         x &= x-1;
      }
      return count;
   }

   /**
    * Returns the maximum coordinate index, starting at 0 for the first coordinate.
    *
    */
   public int maxCoordinate() {
      int c = 0;
      while ((mask >> c) != 0) c++;
      return c - 1;
   }

   /**
    * Returns all subsets of the current coordinate set, whose cardinality is at most \c maxOrder.
    *
    */
   @Override public List<CoordinateSet> subsets(boolean includeEmptySet, int maxOrder) {

      maxOrder = Math.min(maxOrder, maxCoordinate() + 1);

      long maskMax = (1L << (maxCoordinate() + 1));

      List<CoordinateSet> list = new ArrayList<CoordinateSet>();

      for (int order = includeEmptySet ? 0 : 1; order <= maxOrder; order++) {

         // enumeration with Gosper's hack
         // http://home.pipeline.com/~hbaker1/hakmem/hacks.html#item175

         long mask = (1L << order) - 1;
         while (mask < maskMax) {

            // add subset to list
            CoordinateSet cs = new CoordinateSetLong(mask);
            if (containsAll(cs))
               list.add(cs);

            // compute next mask with Gosper's hack
            long u = mask & -mask; // rightmost bit
            long v = mask + u;
            if (v == 0) break;
            mask = v + (((v ^ mask) / u) >> 2);
         }
      }

      //! // inefficient exhaustive enumeration
      //! long maskMin = includeEmptySet ? 0 : 1;
      //! for (long mask = maskMin; mask < maskMax; mask++) {
      //!    CoordinateSet cs = new CoordinateSetLong(mask);
      //!    if (cs.cardinality() <= maxOrder && containsAll(cs))
      //!       list.add(cs);
      //! }
 
      return list;
   }

   /**
    * Returns a set of all coordinates in a space of dimension \c dimension.
    *
    */
   public static CoordinateSet allCoordinates(int dimension) {
      return new CoordinateSetLong((1L << dimension) - 1);
   }
}
