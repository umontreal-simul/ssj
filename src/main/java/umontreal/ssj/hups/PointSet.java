/*
 * Class:        PointSet
 * Description:  Base class of all point sets
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package umontreal.ssj.hups;

import java.util.NoSuchElementException;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.util.Num;
import umontreal.ssj.util.PrintfFormat;

/**
 * This abstract class represents a general point set.
 * It defines some basic methods for accessing and
 * manipulating the points and other properties of the point set. 
 * Conceptually, a point set can be viewed as a
 * two-dimensional array, whose element @f$(i,j)@f$ contains @f$u_{i,j}@f$,
 * the *coordinate* @f$j@f$ of point @f$i@f$, for @f$i\geq 0@f$ and @f$j\geq 0@f$.
 * Each coordinate @f$u_{i,j}@f$
 * is assumed to be in the unit interval @f$[0,1]@f$. Whether the values 0
 * and 1 can occur or not may depend on the actual implementation of the point set.
 * In general, when the points are randomized, after the randomization
 * we add  (internally, when the points or coordinates are requested) 
 * @f$\epsilon/2@f$ to each coordinate for a very small 
 * @f$\epsilon > 0@f$, so that 0 never occurs in the output.
 * By default, in this `PointSet` class, we have @f$\epsilon = 2^{-54}@f$.
 *
 * All points have the same number of coordinates (their <em>dimension</em>)
 * and this number can be queried by  #getDimension. The number of points is
 * queried by  #getNumPoints. The points and coordinates are both numbered
 * starting from 0 and their number can actually be infinite.
 * 
 * The #iterator method provides a @ref PointSetIterator object which can
 * enumerate the points and their coordinates. Several iterators over the
 * same point set can coexist at any given time.  However, in the current implementation
 * they will all enumerate the same randomized points when the points are randomized,
 * because the randomizations are integrated in the point sets and not in the iterators.
 * These iterators are
 * instances of a hidden inner-class that implements the
 * @ref PointSetIterator interface. The default implementation of iterator
 * provided here relies on the method  #getCoordinate to access the
 * coordinates directly. However, this approach is rarely efficient.
 * (One exception is in a @ref CachedPointSet.) 
 * Specialized implementations that dramatically improve the performance are
 * provided in subclasses of  @ref PointSet. The  @ref PointSetIterator
 * interface actually extends the  @ref umontreal.ssj.rng.RandomStream
 * interface, so that the iterator can also be seen as a
 * @ref umontreal.ssj.rng.RandomStream and used wherever such a stream is
 * required for generating uniform random numbers. This permits one to easily
 * replace pseudorandom numbers by the coordinates of a selected set of
 * highly-uniform points, i.e., to replace Monte Carlo by QMC or RQMC
 * in a simulation program.
 *
 * The present abstract class has only one abstract method:  #getCoordinate.
 * Providing an implementation for this method is already sufficient for the
 * subclass to work. However, in most cases, efficiency can be dramatically
 * improved by overwriting  #iterator to provide a custom iterator that does
 * not necessarily rely on  #getCoordinate.  Direct use of
 * #getCoordinate to access the coordinates is generally discouraged. One
 * should access the points and coordinates via the iterators.
 * Built-in range checks when generating the points require some extra time and also assume that
 * nobody ever uses negative indices. If #getCoordinate is never accessed
 * directly by the user, it may be implemented without range checks.
 * 
 * In this abstract class, the #addRandomShift() methods generate a `double[]` array
 * to be used eventually to add a random shift modulo 1, 
 * but this random shift is not used here.
 * It can be used in subclasses when generating coordinates.   
 * In subclasses, #addRandomShift (d1, d2, stream) can also be optionally re-implemented
 * to produce something else than a `double[]` array, 
 * for example when we want the random shift to be digital. 
 * There are also some types of point sets that do not use at all such a random shift.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public abstract class PointSet {

   /** 
    * Since Java has no unsigned type, the 32nd bit cannot be used efficiently,
    * so we have only 31 bits. This mainly affects digit scrambling and bit vectors. 
    * This also limits the maximum number of columns for the generating matrices 
    * of digital nets in base 2.
    */
   protected static final int MAXBITS = 31; // Max number of usable bits.

   /** 
    * To avoid 0 for nextCoordinate when random shifting, we add this to each coordinate.
    */
   protected double EpsilonHalf = 1.0 / Num.TWOEXP[55];  // 2^{-55}

   /** 
    * Dimension of the points.
    */
   protected int dim = 0;                 

   /** 
    * Number of points.
    */
   protected int numPoints = 0;

   /** 
    * Current dimension of the shift.
    * This is useful mostly for the case where the points have an unlimited number of coordinates.
    */
   // **Pierre:** Maybe this could be defined only in `CycleBasedPointSet`.
  protected int dimShift = 0;
   
   /** 
    * Number of array elements in the shift vector, always >= dimShift.
    */
   protected int capacityShift = 0;
   
   /** 
    * This is the shift vector as a `double[]` array, which contains the current random shift
    * in case we apply a random shift modulo 1.
    * It is initially null.  
    */
   protected double[] shift;

   /** 
    * Stream used to generate the random shifts.  This stream is saved from the last
    * time we have called #addRandomShift, in case this method has to be called again 
    * internally by an iterator.  This may happen for example if the points have unbounded 
    * dimension and we need to extend the random shift to additional coordinates.
    */
   protected RandomStream shiftStream;

   /**
    * Returns the dimension (number of available coordinates) of the points.
    * If the dimension is actually infinite,
    * <tt>Integer.MAX_VALUE</tt> is returned.
    * 
    *  @return the dimension of the point set or
    * <tt>Integer.MAX_VALUE</tt> if it is infinite
    */
   public int getDimension() {
      return dim;
   }

   /**
    * Returns the number of points. If this number is actually infinite,
    * <tt>Integer.MAX_VALUE</tt> is returned.
    *  @return the number of points in the point set or
    * <tt>Integer.MAX_VALUE</tt> if the point set has infinity cardinality.
    */
   public int getNumPoints() {
      return numPoints;
   }

   /**
    * Returns @f$u_{i,j}@f$, the coordinate @f$j@f$ of the point @f$i@f$.
    * When the points are randomized (e.g., a random shift is added), the values
    * returned by this method should incorporate the randomizations.
    * **Richard:** La méthode `getCoordinate` de certaines classes
    * ne tient pas compte du random shift, contrairement à l’itérateur de
    * la même classe. Faut-il que toutes les `getCoordinate` implémentent
    * le random shift quand il existe?   Oui.
    * 
    *  @param i            index of the point to look for
    *  @param j            index of the coordinate to look for
    *  @return the value of @f$u_{i,j}@f$
    */
   public abstract double getCoordinate (int i, int j);

   /**
    * Constructs and returns a point set iterator. The default
    * implementation returns an iterator that uses the method
    * {@link #getCoordinate() getCoordinate(i,j)} to iterate over the
    * points and coordinates, but subclasses can reimplement it for better
    * efficiency.
    *  @return point set iterator for the point set
    */
   public PointSetIterator iterator() {
      return new DefaultPointSetIterator();
   }

   /**
    * Randomizes this point set using the given `rand`.
    * Use the equivalent `rand.randomize(this)` instead.
    *  @param rand          @ref PointSetRandomization to use
    */
   // @Deprecated
   public void randomize (PointSetRandomization rand) {
       rand.randomize(this);
       // Note that RandomShift.randomize(p) calls  addRandomShift !
   }

   /**
    * By default, this method generates a random shift in the protected `double[]` array `shift`,
    * to be used eventually for a random shift modulo 1.
    * This random shift will be added (modulo 1) to all the points when they are enumerated 
    * by the iterator.  
    * The random shift is generated using the `stream` to generate the uniform random numbers, 
    * one for each of the coordinates `d1` to `d2-1`.
    * The variable `shiftStream` is also set to this `stream`.
    * Allowing an arbitrary range `d1` to `d2-1` permits one to extends the random shift
    * to additional coordinates when needed, e.g., when the number of coordinates that
    * we need is unknown a priori.   There are also other situations in which we may want
    * to randomize only specific coordinates, e.g., in the Array-RQMC method.
    * 
    * This method is re-implemented in subclasses for which a different form of
    * shift is appropriate, e.g., should produce a random digital shift in base @f$b@f$ for a #DigitalNet,
    * a binary digital shift for a #DigitalNetBase2, etc.
    * 
    * These methods are used internally by randomizations.
    * For example, calling `RandomShift.randomize(PointSet)` calls `addRandomShift`.
    * Normally, one should not call them directly, but use #PointSetRandomization objects instead.
    */
   public void addRandomShift (int d1, int d2, RandomStream stream) {
		if (d1 < 0 || d1 > d2)
			throw new IllegalArgumentException("illegal parameter d1 or d2");
		if (d2 > capacityShift) {
			int d3 = Math.max(4, capacityShift);
			while (d2 > d3)
				d3 *= 2;
			double[] temp = new double[d3];
			capacityShift = d3;
			for (int i = 0; i < d1; i++)
				temp[i] = shift[i];
			shift = temp;
		}
	    shiftStream = stream;	
		dimShift = d2;
		for (int i = d1; i < d2; i++)
			shift[i] = shiftStream.nextDouble();
   }

	/**
	 * Same as {@link #addRandomShift() addRandomShift(0, dim, stream)}, where `dim` is the
	 * dimension of the point set.
	 * 
	 * @param stream
	 *            random number stream used to generate uniforms for the random shift
	 */
	public void addRandomShift(RandomStream stream) {
		addRandomShift(0, dim, stream);
	}


	/**
	 * Refreshes the random shift (generates new uniform values for the random shift coordinates)
	 * for coordinates `d1` to `d2-1`, using the saved `shiftStream`.
	 */
	public void addRandomShift(int d1, int d2) {
		addRandomShift (d1, d2, shiftStream);
	}

	/**
	 * Same as {@link #addRandomShift() addRandomShift(0, dim)}, where `dim` is the
	 * dimension of the point set.
	 */
	public void addRandomShift() {
		addRandomShift(0, dim);
	}

	/**
    * Erases the current random shift, if any.
    */
   public void clearRandomShift() {
      capacityShift = 0;
      dimShift = 0;
//      shiftStream = null;
  }

   /**
    * Formats a string that contains information about the point set.
    *  @return string representation of the point set information
    */
   public String toString() {
       StringBuffer sb = new StringBuffer ("Number of points: ");
       int x = getNumPoints();
       if (x == Integer.MAX_VALUE)
          sb.append ("infinite");
       else
          sb.append (x);
       sb.append (PrintfFormat.NEWLINE + "Point set dimension: ");
       x = getDimension();
       if (x == Integer.MAX_VALUE)
          sb.append ("infinite");
       else
          sb.append (x);
       return sb.toString();
   }

   /**
    * Same as invoking  {@link #formatPoints(int,int) formatPoints(n, d)}
    * with @f$n@f$ and @f$d@f$ equal to the number of points and the
    * dimension of this object, respectively.
    *  @return string representation of all the points in the point set
    *
    *  @exception UnsupportedOperationException if the number of points or
    * dimension of the point set is infinite
    */
   public String formatPoints() {
      PointSetIterator iter = iterator();
      return formatPoints (iter);
   }

   /**
    * Formats a string that displays the same information as returned by
    * #toString, together with the first @f$d@f$ coordinates of the first
    * @f$n@f$ points. If @f$n@f$ is larger than the number of points in
    * the point set, it is reset to that number. If @f$d@f$ is larger than
    * the dimension of the points, it is reset to that dimension. The
    * points are printed in the simplest format, separated by spaces, by
    * calling the default iterator repeatedly.
    *  @param n            number of points
    *  @param d            dimension
    *  @return string representation of first d coordinates of first n
    * points in the point set
    */
   public String formatPoints (int n, int d) {
      PointSetIterator iter = iterator();
      return formatPoints (iter, n, d);
   }

   /**
    * Same as invoking  {@link #formatPoints(PointSetIterator,int,int)
    * formatPoints(iter, n, d)} with @f$n@f$ and @f$d@f$ equal to the
    * number of points and the dimension, respectively.
    *  @param iter         iterator associated to the point set
    *  @return string representation of all the points in the point set
    *
    *  @exception UnsupportedOperationException if the number of points or
    * dimension of the point set is infinite
    */
   public String formatPoints (PointSetIterator iter) {
      int n = getNumPoints();
      if (n == Integer.MAX_VALUE)
         throw new UnsupportedOperationException (
            "Number of points is infinite");
      int d = getDimension();
      if (d == Integer.MAX_VALUE)
         throw new UnsupportedOperationException ("Dimension is infinite");
      return formatPoints (iter, n, d);
   }

   /**
    * Same as invoking  {@link #formatPoints(int,int) formatPoints(n, d)},
    * but prints the points by calling `iter` repeatedly. The order of the
    * printed points may be different than the one resulting from the
    * default iterator.
    *  @param iter         iterator associated to the point set
    *  @param n            number of points
    *  @param d            dimension
    *  @return string representation of first d coordinates of first n
    * points in the point set
    */
   public String formatPoints (PointSetIterator iter, int n, int d) {
      if (getNumPoints() < n)
         n = getNumPoints();
      if (getDimension() < d)
         d = getDimension();
      StringBuffer sb = new StringBuffer (toString());
      sb.append (PrintfFormat.NEWLINE + PrintfFormat.NEWLINE
                 + "Points of the point set:" + PrintfFormat.NEWLINE);
      for (int i=0; i<n; i++) {
        for (int j=0; j<d; j++) {
            sb.append ("  ");
            sb.append (iter.nextCoordinate());
         }
         sb.append (PrintfFormat.NEWLINE);
         iter.resetToNextPoint();
      }
      return sb.toString();
   }

   /**
    * Similar to  {@link #formatPoints() formatPoints()}, but the points
    * coordinates are printed in base @f$b@f$.
    *  @param b            base
    *  @return string representation of all the points in the point set
    *
    *  @exception UnsupportedOperationException if the number of points or
    * dimension of the point set is infinite
    */
   public String formatPointsBase (int b) {
      PointSetIterator iter = iterator();
      return formatPointsBase (iter, b);
   }

   /**
    * Similar to  {@link #formatPoints(int,int) formatPoints(n, d)}, but
    * the points coordinates are printed in base @f$b@f$.
    *  @param n            number of points
    *  @param d            dimension
    *  @param b            base
    *  @return string representation of first d coordinates of first n
    * points in the point set
    */
   public String formatPointsBase (int n, int d, int b) {
      PointSetIterator iter = iterator();
      return formatPointsBase(iter, n, d, b);
   }

   /**
    * Similar to  {@link #formatPoints(PointSetIterator)
    * formatPoints(iter)}, but the points coordinates are printed in base
    * @f$b@f$.
    *  @param iter         iterator associated to the point set
    *  @param b            base
    *  @return string representation of all the points in the point set
    *
    *  @exception UnsupportedOperationException if the number of points or
    * dimension of the point set is infinite
    */
   public String formatPointsBase (PointSetIterator iter, int b) {
      int n = getNumPoints();
      if (n == Integer.MAX_VALUE)
         throw new UnsupportedOperationException (
            "Number of points is infinite");
      int d = getDimension();
      if (d == Integer.MAX_VALUE)
         throw new UnsupportedOperationException ("Dimension is infinite");
      return formatPointsBase (iter, n, d, b);
   }

   /**
    * Similar to  {@link #formatPoints(PointSetIterator,int,int)
    * formatPoints(iter, n, d)}, but the points coordinates are printed in
    * base @f$b@f$.
    *  @param iter         iterator associated to the point set
    *  @param n            number of points
    *  @param d            dimension
    *  @param b            base
    *  @return string representation of first d coordinates of first n
    * points in the point set
    */
   public String formatPointsBase (PointSetIterator iter, int n, int d, int b) {
      if (getNumPoints() < n)
         n = getNumPoints();
      if (getDimension() < d)
         d = getDimension();
      StringBuffer sb = new StringBuffer (toString());
      sb.append (PrintfFormat.NEWLINE + PrintfFormat.NEWLINE
                 + "Points of the point set:" + PrintfFormat.NEWLINE);
      double x;
      int acc = 10;
      if (b == 2)
         acc = 20;
      else if (b == 3)
         acc = 13;
      else
         acc = 10;
      int width = acc + 3;
      String chaine;
      for (int i=0; i<n; i++) {
        for (int j=0; j<d; j++) {
            sb.append ("  ");
            x = iter.nextCoordinate();
            chaine = PrintfFormat.formatBase (-width, acc, b, x);
            sb.append (chaine);
         }
         sb.append (PrintfFormat.NEWLINE);
         iter.resetToNextPoint();
      }
      return sb.toString();
   }

   /**
    * Same as invoking  {@link #formatPointsNumbered(int,int)
    * formatPointsNumbered(n, d)} with @f$n@f$ and @f$d@f$ equal to the
    * number of points and the dimension, respectively.
    *  @return string representation of all the points in the point set
    *
    *  @exception UnsupportedOperationException if the number of points or
    * dimension of the point set is infinite
    */
   public String formatPointsNumbered() {
      int n = getNumPoints();
      if (n == Integer.MAX_VALUE)
         throw new UnsupportedOperationException (
            "Number of points is infinite");
      int d = getDimension();
      if (d == Integer.MAX_VALUE)
         throw new UnsupportedOperationException ("Dimension is infinite");
      return formatPointsNumbered (n, d);
   }

   /**
    * Same as invoking  {@link #formatPoints(int,int) formatPoints(n,d)},
    * except that the points are numbered.
    *  @param n            number of points
    *  @param d            dimension
    *  @return string representation of first d coordinates of first n
    * points in the point set
    */
   public String formatPointsNumbered (int n, int d) {
      if (getNumPoints() < n)
         n = getNumPoints();
      if (getDimension() < d)
         d = getDimension();
      StringBuffer sb = new StringBuffer (toString());
      PointSetIterator itr = iterator();
      sb.append (PrintfFormat.NEWLINE + PrintfFormat.NEWLINE
                 + "Points of the point set:");
      for (int i=0; i<n; i++) {
         sb.append (PrintfFormat.NEWLINE + "Point " +
    //                itr.getCurPointIndex() + " = (");
                                           i + "  =  (");
         boolean first = true;
         for (int j=0; j<d; j++) {
            if (first)
               first = false;
            else
               sb.append (", ");
            sb.append (itr.nextCoordinate());
         }
         sb.append (")");
         itr.resetToNextPoint();
      }
      return sb.toString();
   }



// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   /**
    * This class implements a default point set iterator.
    * Since it is inherited by subclasses, it can be used as a base class 
    * for iterators.  It is implemented as an inner class because it can then use 
    * directly the variables of the PointSet class. It would be more difficult and
    * cumbersome to access those variables if it was implemented as a
    * separate class.
    */ 
   protected class DefaultPointSetIterator implements PointSetIterator {

	   /**
	    * Index of the current point.
	    */ 
	  protected int curPointIndex = 0;

	   /**
	    * Index of the current coordinate.
	    */ 
      protected int curCoordIndex = 0;

	   /**
	    * Default constant epsilon/2 added to the points after a random shift.
	    */ 
      protected double EpsilonHalf = 1.0 / Num.TWOEXP[55];

	   /**
	    * Error message for index out of bounds.
	    */ 
      protected void outOfBounds () {
         if (getCurPointIndex() >= numPoints)
            throw new NoSuchElementException ("Not enough points available");
         else
            throw new NoSuchElementException ("Not enough coordinates available");
      }

	   /**
	    * Set current coordinate to j.
	    */ 
      public void setCurCoordIndex (int j) {
         curCoordIndex = j;
      }

	   /**
	    * Set current coordinate to 0.
	    */ 
      public void resetCurCoordIndex() {
         setCurCoordIndex (0);
      }

	   /**
	    * @return index of current coordinate.
	    */ 
      public int getCurCoordIndex() {
        return curCoordIndex;
      }

	   /**
	    * @return true of the current point has another coordinate available.
	    */ 
      public boolean hasNextCoordinate() {
        return getCurCoordIndex() < getDimension();
      }

	   /**
	    * @return the next coordinate
	    */ 
      public double nextCoordinate() {
         if (getCurPointIndex() >= numPoints || getCurCoordIndex() >= dim)
            outOfBounds();
         return getCoordinate (curPointIndex, curCoordIndex++);
      }

	   /**
	    * @return in `p` the `d` next coordinates
	    */ 
      public void nextCoordinates (double p[], int d)  {
         if (getCurCoordIndex() + d > getDimension()) outOfBounds();
         for (int j = 0; j < d; j++)
            p[j] = nextCoordinate();
      }

	   /**
	    * Resets the current point index to `i` and current coordinate to 0.
	    */ 
      // This is called with i = numPoints when nextPoint generates the
      // last point, so i = numPoints must be allowed.
      // The "no more point" error will be raised if we ask for
      // a new coordinate or point when  i = numPoints.
      public void setCurPointIndex (int i) {
         curPointIndex = i;
         resetCurCoordIndex();
      }

	   /**
	    * Resets both the current point index and the current coordinate to 0.
	    */ 
      public void resetCurPointIndex() {
         setCurPointIndex (0);
      }

	   /**
	    * Resets the current point index to the next one and current coordinate to 0,
	    * and returns the new current point index.
	    */ 
      public int resetToNextPoint() {
         setCurPointIndex (curPointIndex + 1);
         return curPointIndex;
      }

	   /**
	    * @return the current point index.
	    */ 
      public int getCurPointIndex() {
        return curPointIndex;
      }

	   /**
	    * @return `true` iff the current point is not the last one.
	    */ 
      public boolean hasNextPoint() {
        return getCurPointIndex() < getNumPoints();
      }

	   /**
	    * Returns in `p` a block of `d` successive coordinates of the current point, 
	    * starting at coordinate `fromDim`, and advances to the next point.
	    * @return the new current point index.
	    */ 
      public int nextPoint (double p[], int fromDim, int d) {
         setCurCoordIndex(fromDim);
         nextCoordinates (p, d);
         return resetToNextPoint();
      }

	   /**
	    * Same as {@link #nextPoint(double[],int,int) nextPoint(p, 0, d)}.
	    */ 
      public int nextPoint (double p[], int d) {
         resetCurCoordIndex();
         nextCoordinates (p, d);
         return resetToNextPoint();
      }

	   /**
	    * Same as #resetCurPointIndex().
	    */ 
      public void resetStartStream() {
         resetCurPointIndex();
      }

	   /**
	    * Same as #resetCurCoordIndex().
	    */ 
      public void resetStartSubstream() {
         resetCurCoordIndex();
      }

	   /**
	    * Same as #resetToNextPoint().
	    */ 
      public void resetNextSubstream() {   // Same as resetToNextPoint();
         resetToNextPoint();
      }

	   /**
	    * Not implemented here, raises an exception.  Must be here for compatibility with the 
	    * @ref RandomStream  interface.
	    */ 
      public void setAntithetic (boolean b) {
         throw new UnsupportedOperationException();
      }

	   /**
	    * Same as #nextCoordinate()
	    */ 
      public double nextDouble() {
         return nextCoordinate();
      }

	   /**
	    * Returns in `p[start..start+n-1]` a block of `n` successive coordinates of the current point, 
	    * obtained by calling #nextDouble() `n` times. 
	    */ 
      public void nextArrayOfDouble (double[] u, int start, int n) {
         if (n < 0)
            throw new IllegalArgumentException ("n must be positive.");
         for (int i = start; i < start+n; i++)
            u[i] = nextDouble();
      }

	   /**
	    * Similar to #nextDouble(), but returns an integer uniformly distributed in `[i..j]`.
	    */ 
      public int nextInt (int i, int j) {
         return (i + (int)(nextDouble() * (j - i + 1.0)));
      }

	   /**
	    * Similar to #nextArrayOfDouble but returns in `u[start..start+n-1]` a block of `n` 
	    * integers uniformly distributed in `[i..j]`.
	    */ 
      public void nextArrayOfInt (int i, int j, int[] u, int start, int n) {
         if (n < 0)
            throw new IllegalArgumentException ("n must be positive.");
         for (int k = start; k < start+n; k++)
            u[k] = nextInt (i, j);
      }

	   /**
	    * @return a printable `String` that gives the current point index and current coordinate index.
	    */ 
      public String formatState() {
         return "Current point index: " + getCurPointIndex() +
              PrintfFormat.NEWLINE + "Current coordinate index: " +
                  getCurCoordIndex();
      }
   }
}