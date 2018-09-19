/*
 * Class:        TallyHistogram
 * Description:  Histogram of a tally
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        January 2011
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
import java.util.logging.Level;
import java.util.logging.Logger;
import umontreal.ssj.util.PrintfFormat;

/**
 * This class is an extension of  @ref Tally which gives a more detailed view
 * of the observations statistics. The individual observations are assumed to
 * fall into different bins (boxes) of equal width on an interval. The total
 * number of observations falling into the bins are kept in an array of
 * counters. This is useful, for example, if one wish to build a histogram
 * from the observations. One must access the array of bin counters to
 * compute quantities not supported by the methods in  @ref Tally.
 *
 * *Never add or remove observations directly* on the array of bin counters
 * because this would put the  @ref Tally counters in an inconsistent state.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class TallyHistogram extends Tally {
   private int[] co;         // counter: num of values in bin[i]
   private int numBins;      // number of bins
   private double m_h;       // width of 1 bin
   private double m_a;       // left boundary of first bin
   private double m_b;       // right boundary of last bin
   private Logger log = Logger.getLogger ("umontreal.ssj.stat");

   /**
    * Constructs a `TallyHistogram` statistical probe. Divide the interval
    * @f$[a,b]@f$ into @f$s@f$ bins of equal width and initializes a
    * counter to 0 for each bin. Whenever an observation falls into a bin,
    * the bin counter is increased by 1. There are two extra bins (and
    * counters) that count the number of observations @f$x@f$ that fall
    * outside the interval @f$[a,b]@f$: one for those @f$x< a@f$, and the
    * other for those @f$x > b@f$.
    *  @param a            left boundary of interval
    *  @param b            right boundary of interval
    *  @param s            number of bins
    */
   public TallyHistogram(double a, double b, int s) {
      super();
      init (a, b, s);
   }

   /**
    * Constructs a new `TallyHistogram` statistical probe with name
    * `name`.
    *  @param name         the name of the tally.
    *  @param a            left boundary of interval
    *  @param b            right boundary of interval
    *  @param s            number of bins
    */
   public TallyHistogram (String name, double a, double b, int s) {
      super (name);
      init (a, b, s);
   }

   /**
    * Initializes this object. Divide the interval @f$[a,b]@f$ into
    * @f$s@f$ bins of equal width and initializes all counters to 0.
    *  @param s            number of bins
    *  @param a            left boundary of interval
    *  @param b            right boundary of interval
    */
   public void init (double a, double b, int s) {
      /* The counters co[1] to co[s] contains the number of observations
         falling in the interval [a, b].
         co[0] is the number of observations < a,
         and co[s+1] is the number of observations > b.
      */

      super.init();
      if (b <= a)
         throw new IllegalArgumentException ("   b <= a");
      co = new int[s + 2];
      numBins = s;
      m_h = (b - a) / s;
      m_a = a;
      m_b = b;
      for (int i = 0; i <= s + 1; i++)
         co[i] = 0;
   }

   /**
    * Gives a new observation @f$x@f$ to the statistical collectors.
    * Increases by 1 the bin counter in which value @f$x@f$ falls. Values
    * that fall outside the interval @f$[a,b]@f$ are added in extra bin
    * counter bin[0] if @f$x < a@f$, and in bin[@f$s+1@f$] if @f$x > b@f$.
    *  @param x            observation value
    */
   public void add (double x) {
      super.add(x);
      if (x < m_a)
        ++co[0];
      else if (x > m_b)
        ++co[1 + numBins];
      else {
         int i = 1 + (int) ((x - m_a) / m_h);
         ++co[i];
      }
   }

   /**
    * Returns the bin counters. Each counter contains the number of
    * observations that fell in its corresponding bin. The counters
    * bin[@f$i@f$], @f$i=1, 2, …, s@f$ contain the number of observations
    * that fell in each subinterval of @f$[a,b]@f$. Values that fell
    * outside the interval @f$[a,b]@f$ were added in extra bin counter
    * bin[0] if @f$x < a@f$, and in bin[@f$s+1@f$] if @f$x > b@f$. There
    * are thus @f$s+2@f$ counters.
    *  @return the array of counters
    */
   public int[] getCounters() {
      return co;
   }

   /**
    * Returns the number of bins @f$s@f$ dividing the interval
    * @f$[a,b]@f$. Does not count the two extra bins for the values of
    * @f$x<a@f$ or @f$x>b@f$.
    *  @return the number of bins
    */
   public int getNumBins() {
      return numBins;
   }

   /**
    * Returns the left boundary @f$a@f$ of interval @f$[a,b]@f$.
    *  @return left boundary of interval
    */
   public double getA() {
      return m_a;
   }

   /**
    * Returns the right boundary @f$b@f$ of interval @f$[a,b]@f$.
    *  @return right boundary of interval
    */
   public double getB() {
      return m_b;
   }

   /**
    * Clones this object and the array which stores the counters.
    */
   public TallyHistogram clone() {
      TallyHistogram image = (TallyHistogram)super.clone();
      int[] coco = new int[2 + numBins];
      System.arraycopy (co, 0, coco, 0, 2 + numBins);
      image.co = coco;
      image.m_h = m_h;
      image.m_a = m_a;
      image.m_b = m_b;
      image.numBins = numBins;
      return image;
   }

   /**
    * Returns the bin counters as a `String`.
    */
   public String toString() {
      StringBuffer sb = new StringBuffer ();
      sb.append ("---------------------------------------" +
                PrintfFormat.NEWLINE);
      sb.append (name + PrintfFormat.NEWLINE);
      sb.append ("Interval = [ " + m_a + ", " + m_b + " ]" +
                 PrintfFormat.NEWLINE);
      sb.append ("Number of bins = " + numBins + " + 2" + PrintfFormat.NEWLINE);
      sb.append (PrintfFormat.NEWLINE + "Counters = {" +
                 PrintfFormat.NEWLINE);
      sb.append ("   (-inf, " + PrintfFormat.f(6, 3, m_a)
                 + ")    " + co[0] + PrintfFormat.NEWLINE);
      for (int i = 1; i <= numBins; i++) {
         double a = m_a + (i-1)*m_h;
         double b = m_a + i*m_h;
         sb.append ("   (" +
            PrintfFormat.f(6, 3, a) + ", " +
            PrintfFormat.f(6, 3, b) + ")    " + co[i] +
                 PrintfFormat.NEWLINE);
      }
      sb.append ("   (" + PrintfFormat.f(6, 3, m_b)
                 + ", inf)    " + co[numBins + 1] +
                 PrintfFormat.NEWLINE);
      sb.append ("}" + PrintfFormat.NEWLINE);
      return sb.toString();
   }

}