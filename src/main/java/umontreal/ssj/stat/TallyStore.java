/*
 * Class:        TallyStore
 * Description:
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
package umontreal.ssj.stat;
import cern.colt.list.DoubleArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import umontreal.ssj.util.PrintfFormat;

/**
 * This class is a variant of  @ref Tally for which the individual
 * observations are stored in a list implemented as a  DoubleArrayList.
 * The class `DoubleArrayList` is imported from the COLT library and provides
 * an efficient way of storing and manipulating a list of real-valued numbers
 * in a dynamic array.
 * The  DoubleArrayList object used to store the values can be either passed
 * to the constructor or created by the constructor, and can be accessed via
 * the  #getDoubleArrayList method.
 *
 * The same counters as in  @ref Tally are maintained and are used by the
 * inherited methods. One must access the list of observations to compute
 * quantities not supported by the methods in  @ref Tally, and/or to use
 * methods provided by the COLT package.
 *
 * *Never add or remove observations directly* on the  DoubleArrayList
 * object, because this would put the counters of the `TallyStore` object in
 * an inconsistent state.
 *
 * There are two potential reasons for using a  @ref TallyStore object
 * instead of directly using a  DoubleArrayList object: (a) it can broadcast
 * observations and (b) it maintains a few additional counters that may speed
 * up some operations such as computing the average.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class TallyStore extends Tally {

   private DoubleArrayList array = null;  // Where the observations are stored.
   private Logger log = Logger.getLogger ("umontreal.ssj.stat");

   /**
    * Constructs a new `TallyStore` statistical probe.
    */
   public TallyStore() {
      super();
      array = new DoubleArrayList();
   }

   /**
    * Constructs a new `TallyStore` statistical probe with name `name`.
    *  @param name         the name of the tally.
    */
   public TallyStore (String name) {
      super (name);
      array = new DoubleArrayList();
   }

   /**
    * Constructs a new `TallyStore` statistical probe with given initial
    * capacity `capacity` for its associated array.
    *  @param capacity     initial capacity of the array of observations
    */
   public TallyStore (int capacity) {
      super();
      array = new DoubleArrayList (capacity);
   }

   /**
    * Constructs a new `TallyStore` statistical probe with name `name` and
    * given initial capacity `capacity` for its associated array.
    *  @param name         the name of the tally.
    *  @param capacity     initial capacity of the array of observations
    */
   public TallyStore (String name, int capacity) {
      super (name);
      array = new DoubleArrayList (capacity);
   }

   /**
    * Constructs a new `TallyStore` statistical probe with given
    * associated array. This array must be empty.
    *  @param a            array that will contain observations
    */
   public TallyStore (DoubleArrayList a) {
      super();
      array = a;
      array.clear();
   }

   public void init() {
       super.init();
       // We must call super before any actions inside constructors.
       // Unfortunately, the base class calls init, which would
       // result in a NullPointerException.
       if (array != null)
          array.clear();
   }
   
   /**
    * Adds one observation `x` to this probe.
    */
   public void add (double x) {
      if (collect) array.add (x);
      super.add(x);
   }

  /**
    * Returns the observations stored in this probe.
    * @return the array of observations associated with this object
    */
   public double[] getArray() {
      array.trimToSize();
      return array.elements();
   }

   /**
    * Returns the  DoubleArrayList object that contains the observations
    * for this probe. **WARNING:** In early releases, this function was
    * named `getArray`.
    *  @return the array of observations associated with this object
    */
   public DoubleArrayList getDoubleArrayList() {
      array.trimToSize();
      return array;
   }

   /**
    * Sorts the elements of this probe using the `quicksort` from Colt.
    */
   public void quickSort() {
       array.quickSort();
   }

   /**
    * Returns the sample covariance of the observations contained in this
    * tally, and the other tally `t2`. Both tallies must have the same
    * number of observations. This returns `Double.NaN` if the tallies do
    * not contain the same number of observations, or if they contain less
    * than two observations.
    *  @param t2           the other tally.
    *  @return the sample covariance.
    */
   public double covariance (TallyStore t2) {
      if (numberObs() != t2.numberObs()) {
         // System.err.println ("******* TallyStore.covariance(): " +
         // "Tally's with different number of observations");
         log.logp (Level.WARNING, "TallyStore", "covariance",
            "This tally, with name " + getName() + ", contains " + numberObs() +
            " observations while " + "the given tally, with name " +
            t2.getName() + ", contains " + t2.numberObs() + "observations");
         return Double.NaN;
      }

      if (numberObs() < 2 || t2.numberObs() < 2) {
         //System.err.println ("******* TallyStore.covariance()   with " +
         // numberObs() + " Observation");
         log.logp (Level.WARNING, "TallyStore", "covariance",
            "This tally, with name " + getName() + ", contains " + numberObs() + " observation");
         return Double.NaN;
      }

      return cern.jet.stat.Descriptive.covariance (
          getDoubleArrayList(), t2.getDoubleArrayList());
   }

   /**
    * Clones this object and the array which stores the observations.
    */
   public TallyStore clone() {
      TallyStore t = (TallyStore)super.clone();
      t.array = (DoubleArrayList)array.clone();
      return t;
   }
   
   /**
    * Returns a new `TallyStore` instance that contains all the observations of this `TallyStore`
    * than are in the interval (a, b).  This method does not sort the observations. 
    *
    * @return a new `TallyStore` object with the selected observations
    */
   public TallyStore extractSubrange (double a, double b) {
       int numObs = this.numberObs();
       double[] obs = this.getArray();
       double x;
       TallyStore t = new TallyStore ();
       for (int i = 0; i < numObs; i++) {
    	   x = obs[i];
    	   if ((x > a) & (x < b))  t.add(x);
       }
       return t;
   } 

   
    /**
     * Returns a new `TallyStore` instance that contains aggregate observations from this `TallyStore`.
     * This method divides the observations in blocks of size `gsize` successive observations,
     * compute the average of each group, and then inserts the average of each group in a new `TallyStore` object.
     *
     * Note that this method does not sort the observations. It will aggregate
     * the observations according to their actual order in the array.
     * To aggregate the sorted observations (e.g., to plot a approximation of the cdf), the user should sort
     * this `TallyStore` object before calling this method.
     *
     * @param gsize the group size to use when performing the aggregation
     *
     * @return a new `TallyStore` object with aggregated observations
     */
    public TallyStore aggregate (int gsize) {
        int numObs = this.numberObs();
        int numGroups = numObs / gsize;
        double[] obs = this.getArray();
        double sum;
        TallyStore t = new TallyStore (numGroups);
        for (int i = 0; i < numGroups; i++) {
            sum = 0.0;
            for (int j = 0; j < gsize; j++)
                sum += obs[gsize * i + j];
            sum /= gsize;
            t.add(sum);
        }
        // This is if gsize does not divide numObs.
        int rest = numObs - numGroups * gsize;
        if (rest > 0) {
            sum = 0.0;
            for (int j = 0; j < rest; j++)
                sum += obs[gsize * numGroups + j];
            sum /= rest;
            t.add(sum);
        }
        return t;
    } 

 
   /**
    * Returns the observations stored in this object as a `String`.
    */
   public String toString() {
      StringBuffer sb = new StringBuffer ();
      for (int i=0; i<numberObs(); i++)
         sb.append (i + "    " + array.getQuick(i) +
                           PrintfFormat.NEWLINE);
      return sb.toString();
   }

}
