/*
 * Class:        Discrepancy
 * Description:  Base class of all discrepancies
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        January 2009

 * SSJ is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, either version 3 of the License, or
 * any later version.

 * SSJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * A copy of the GNU General Public License is available at
   <a href="http://www.gnu.org/licenses">GPL licence site</a>.
 */
package umontreal.ssj.discrepancy;
   import umontreal.ssj.hups.*;
   import cern.colt.list.DoubleArrayList;
   import umontreal.ssj.util.PrintfFormat;

/**
 * This *abstract* class is the base class of all discrepancy classes. All
 * derived classes must implement the abstract method
 * {@link #compute(double[][],int,int) compute(points, n, s)}. <div
 * class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public abstract class Discrepancy {
   protected static double[] ONES = { 1 }; // vector of 1, dimension dim
   protected double[] gamma;        // gamma[]: weights, dimension dim
   protected double[][] Points;     // Points[n][s]:  n points in s dimensions
   protected int dim;               // Dimension of points
   protected int numPoints;         // Number of points
//   protected boolean primeF;                  // true if numPoints is prime
//   protected boolean power2F;                 // true if numPoints is a power of 2

   static final double UNSIX = 1.0/6.0;
   static final double QUARAN = 1.0/42.0;
   static final double UNTRENTE = 1.0 / 30.0;
   static final double DTIERS = 2.0 / 3.0;
   static final double STIERS = 7.0 / 3.0;
   static final double QTIERS = 14.0 / 3.0;


   static protected void setONES (int s) {
      if (s < ONES.length)
         return;
      ONES = new double[s + 1];
      for (int i = 0; i <= s; i++)
         ONES[i] = 1.0;
   }


   protected void appendGamma (StringBuffer sb, double[] gamma, int s) {
      // append the first s components of gamma to sb
      if (gamma != null) {
         for (int j = 0; j < s; ++j)
            sb.append ("  " + gamma[j]);
      }
   }


   private void initD(double[][] points, int n, int s, double[] gamma) {
      setONES (s);
      if (gamma == null)
         setGamma (ONES, s);
     else
         setGamma (gamma, s);
      dim = s;
      numPoints = n;
      this.Points = points;
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in @f$s@f$
    * dimensions. `points[i][j]` is the @f$j@f$-th coordinate of point
    * @f$i@f$. Both @f$i@f$ and @f$j@f$ start at 0. One may also choose
    * `points = null` in which case, the points must be set later.
    */
   public Discrepancy (double[][] points, int n, int s) {
      initD(points, n, s, null);
   }

   /**
    * Constructor with the @f$n@f$ points `points[i]` in @f$s@f$
    * dimensions and the @f$s@f$ weight factors
    * <tt>gamma[</tt>@f$j@f$<tt>]</tt>, @f$j = 0, 1, …, (s-1)@f$.
    * `points[i][j]` is the @f$j@f$-th coordinate of point @f$i@f$. Both
    * @f$i@f$ and @f$j@f$ start at 0. One may also choose `points = null`
    * in which case, the points must be set later.
    */
   public Discrepancy (double[][] points, int n, int s, double[] gamma) {
      initD(points, n, s, gamma);
   }

   /**
    * The number of points is @f$n@f$, the dimension @f$s@f$, and the
    * @f$s@f$ weight factors are <tt>gamma[</tt>@f$j@f$<tt>]</tt>, @f$j =
    * 0, 1, …, (s-1)@f$. The @f$n@f$ points will be chosen later.
    */
   public Discrepancy (int n, int s, double[] gamma) {
      initD(null, n, s, gamma);
   }

   /**
    * Constructor with the point set `set`. All the points are copied in
    * an internal array.
    */
   public Discrepancy (PointSet set) {
      numPoints = set.getNumPoints();
      dim = set.getDimension();
      Points = toArray (set);
      initD(Points, numPoints, dim, null);
   }

   /**
    * Empty constructor. The points and parameters *must* be defined
    * before calling methods of this or derived classes.
    */
   public Discrepancy() {
   }

   /**
    * Computes the discrepancy of all the points in maximal dimension
    * (dimension of the points).
    */
   public double compute() {
      return compute(Points, numPoints, dim, gamma);
   }

   /**
    * Computes the discrepancy of all the points in dimension @f$s@f$.
    */
   public double compute (int s) {
      return compute(Points, numPoints, s, gamma);
   }

   /**
    * Computes the discrepancy of the first `n` points of `points` in
    * dimension `s` with weights `gamma`.
    */
   public double compute (double[][] points, int n, int s, double[] gamma) {
       return -1;
  }

   /**
    * Computes the discrepancy of the first `n` points of `points` in
    * dimension `s` with weights @f$=1@f$.
    */
   public abstract double compute (double[][] points, int n, int s);

   /**
    * Computes the discrepancy of all the points of `points` in maximum
    * dimension. Calls method
    * {@link #compute(double[][],int,int,double[]) compute(points,
    * points.length, points[0].length, gamma)}.
    */
   public double compute (double[][] points) {
      return compute(points, points.length, points[0].length, gamma);
   }

   /**
    * Computes the discrepancy of the first `n` points of `T` in 1
    * dimension. Copies the points in an array of arrays and calls method
    * {@link #compute() compute(double[][], n, 1)}. It should be
    * reimplemented in subclasses for better efficiency.
    */
   public double compute (double[] T, int n) {
      double[][] p = new double[n][1];
      for(int i = 0; i < n; ++i)
         p[i][0] = T[i];
      return compute(p, n, 1);
   }

   /**
    * Computes the discrepancy of all the points of `T` in 1 dimension.
    * Calls method  {@link #compute(double[],int,double) compute(T,
    * T.length, gamma[0])}.
    */
   public double compute (double[] T) {
      return compute(T, T.length, gamma[0]);
   }

   /**
    * Computes the discrepancy of the first `n` points of `T` in 1
    * dimension with weight `gamma`.
    */
   public double compute (double[] T, int n, double gamma) {
        return -1;
   }

   /**
    * Computes the discrepancy of all the points in `set` in the same
    * dimension as the point set and with weights `gamma`.
    */
   public double compute (PointSet set, double[] gamma) {
      int n = set.getNumPoints();
      int dim = set.getDimension();
      if (dim > 1) {
         double[][] points = new double[n][dim];
         PointSetIterator it = set.iterator();
         for(int i=0; i<n; ++i) {
            it.nextPoint(points[i],dim);
         }
         return compute(points,n,dim, gamma);
      } else {
         double[] points = new double[n];
         PointSetIterator it = set.iterator();
         for(int i=0; i<n; ++i){
            points[i] = it.nextCoordinate();
            it.resetToNextPoint();
         }
         return compute(points,n, gamma[0]);
      }
   }

   /**
    * Computes the discrepancy of all the points in `set` in the same
    * dimension as the point set. All the weights @f$=1@f$.
    */
   public double compute (PointSet set) {
      int dim = set.getDimension();
      setONES (dim);
      return compute (set, ONES);
   }

   /**
    * Returns the number of points @f$n@f$.
    */
   public int getNumPoints() {
      return numPoints;
   }

   /**
    * Returns the dimension of the points @f$s@f$.
    */
   public int getDimension() {
      return dim;
   }

   /**
    * Sets the points to `points` and the dimension to @f$s@f$. The number
    * of points is @f$n@f$.
    */
   public void setPoints (double[][] points, int n, int s) {
      initD(points, n, s, null);
    }

   /**
    * Sets the points to `points`. The number of points and the dimension
    * are the same as in `points`.
    */
   public void setPoints (double[][] points) {
      setPoints (points, points.length, points[0].length);
    }

   /**
    * Sets the weight factors to `gam` for each dimension up to @f$s@f$.
    */
   public void setGamma (double[] gam, int s) {
      gamma = gam;
  }

   /**
    * Returns the weight factors `gamma` for each dimension up to @f$s@f$.
    */
   public double[] getGamma () {
      return gamma;
  }

   /**
    * Returns all the @f$n@f$ points (@f$s@f$-dimensional) of
    * @ref umontreal.ssj.hups.PointSet `set` as an array
    * <tt>points[</tt>@f$n@f$<tt>][</tt>@f$s@f$<tt>]</tt>.
    */
   public static double[][] toArray (PointSet set) {
      int n = set.getNumPoints();
      int dim = set.getDimension();
      if (dim > 1) {
         double[][] points = new double[n][dim];
         PointSetIterator it = set.iterator();
         for(int i=0; i<n; ++i)
            it.nextPoint(points[i],dim);
         return points;

      } else {
         double[][] po1 = new double[n][1];
         PointSetIterator it = set.iterator();
         for(int i=0; i<n; ++i) {
            po1[i][0] = it.nextCoordinate();
            it.resetToNextPoint();
         }
         return po1;
      }
   }

   /**
    * Sorts the first @f$n@f$ points of @f$T@f$. Returns the sorted
    * points. *Warning:* @f$T@f$ is sorted also.
    */
   public static DoubleArrayList sort (double[] T, int n) {
      DoubleArrayList sortedU = new DoubleArrayList(T);
      sortedU.sortFromTo (0, n - 1);
      return sortedU;
    }

   /**
    * Returns the parameters of this class.
    */
   public String toString()
   {
      StringBuffer sb = new StringBuffer (getName() + ":" +
                                          PrintfFormat.NEWLINE);
      sb .append ("n = " + numPoints + ",   dim = " + dim +
                  PrintfFormat.NEWLINE);
      sb .append ("gamma = [");
      appendGamma (sb, gamma, dim);
      sb.append (" ]" + PrintfFormat.NEWLINE);
      return sb.toString();
   }

   /**
    * Returns all the points of this class.
    */
   public String formatPoints() {
      StringBuffer sb = new StringBuffer ("Points = [" +
                                          PrintfFormat.NEWLINE);
      if (Points != null) {
         for (int i = 0; i < Points.length; ++i) {
            sb.append (" [ " + Points[i][0]);
            for (int j = 1; j < Points[0].length; ++j)
               sb.append (", " + Points[i][j]);
            sb.append (" ]" + PrintfFormat.NEWLINE);
         }
      }
      sb.append (" ]" + PrintfFormat.NEWLINE);
      return sb.toString();
   }

   /**
    * Returns the name of the  @ref Discrepancy.
    */
   public String getName() {
      return getClass().getSimpleName();
   }

}