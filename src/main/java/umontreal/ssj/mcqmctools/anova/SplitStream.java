package umontreal.ssj.mcqmctools.anova;

import umontreal.ssj.rng.*;


/**
 * Implements a random stream that mixes two input streams by using a coordinate mask.
 *
 */
public class SplitStream extends RandomStreamBase {

   protected int curCoordIndex;
   protected CoordinateSet coords;
   protected double[] vals;
   
   /**
    * Reads 2 * nCache values from a a stream and stores them for future use.
    *
    * When nextValue() is called, a value is popped from the 2 * nCache
    * cached values at even indices for the coordinates contained in the
    * coordinate set \c coords, and at odd indices otherwise.
    *
    */
   public SplitStream(RandomStream stream, int nCache) {
      this.vals = new double[2 * nCache];
      this.curCoordIndex = 0;
      this.coords = null;
      stream.nextArrayOfDouble(vals, 0, 2 * nCache);
   }
   
   public SplitStream(double[] vals) {
      this.curCoordIndex = 0;
      this.coords = null;
      this.vals = vals;
   }

   public SplitStream clone() {
      SplitStream s = new SplitStream(vals.clone());
      s.curCoordIndex = curCoordIndex;
      s.coords = coords;
      return s;
   }
   
   public void setCoordinates(CoordinateSet coords) {
      this.coords = coords;
   }
   
   public CoordinateSet getCoordinates() {
      return coords;
   }
   
   @Override protected double nextValue() {
      int a = (coords != null && coords.contains(curCoordIndex)) ? 0 : 1;
      return vals[2*(curCoordIndex++) + a];
      //return vals[a*vals.length/2 + (curCoordIndex++)];
   }

   @Override public void resetNextSubstream() {
      throw new UnsupportedOperationException();
   }
   
   @Override public void resetStartStream() {
      throw new UnsupportedOperationException();
   }
   
   @Override public void resetStartSubstream() {
      curCoordIndex = 0;
   }

   @Override public String toString() {
      return getClass().getSimpleName() + " [nCache=" + (vals.length/2) + "]";
   }
}

