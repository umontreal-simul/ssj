/*
 * Class:        StatProbe
 * Description:  statistical probe
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
import java.util.List;
import java.util.ArrayList;
import umontreal.ssj.util.PrintfFormat;

/**
 * The objects of this class are *statistical probes* or *collectors*, which
 * are elementary devices for collecting statistics. Each probe collects
 * statistics on a given variable. The subclasses  @ref Tally,
 * @ref TallyStore, and  @ref umontreal.ssj.simevents.Accumulate (from
 * package  @ref umontreal.ssj.simevents ) implement two kinds of probes, for
 * the case of successive observations @f$X_1,X_2,X_3,…@f$, and for the case
 * of a variable whose value evolves in time, respectively.
 *
 * Each instance of  @ref StatProbe contains a list of
 * @ref ObservationListener that can listen to individual observations. When
 * a probe is updated, i.e., receives a new statistical observation, it
 * broadcasts this new data to all registered observers. The broadcasting of
 * observations to registered observers can be turned ON or OFF at any time.
 * It is initially OFF by default and should stay OFF when there are no
 * registered observers, to avoid unnecessary overhead.
 *
 * The data collection by the statistical probe itself can also be turned ON
 * or OFF. By default, it is initially ON. We can turn it OFF, for example,
 * if we want to use the statistical probe only to pass data to the
 * observers, and do not need it to store any information.
 *
 * In the simplest programs, collection is ON, broadcast is OFF, and the
 * overall stats are accessed via the methods `min`, `max`, `sum`, `average`,
 * ... of the collector.
 *
 * <div class="SSJ-bigskip"></div>
 */
public abstract class StatProbe {

   private List<ObservationListener> listeners = new ArrayList<ObservationListener>();
   protected String name;
   protected double maxValue;
   protected double minValue;
   protected double sumValue;
   protected boolean collect = true;
   protected boolean broadcast = false;
   protected boolean showNobs = true;

   /**
    * Initializes the statistical collector.
    */
   abstract public void init();

   /**
    * Sets the name of this statistical collector to `name`.
    */
   public void setName (String name) {
      this.name = name;
   }

   /**
    * Returns the name associated with this probe, or `null` if no name
    * was specified upon construction.
    *  @return the name associated to this collector, or `null` if not
    * specified
    */
   public String getName() {
      return name;
   }

   /**
    * Returns the smallest value taken by the variable since the last
    * initialization of this probe. This returns
    * <tt>Double.POSITIVE_INFINITY</tt> if the probe was not updated since
    * the last initialization.
    *  @return the smallest value taken by the collector since last
    * initialization
    */
   public double min() {
      return minValue;
   }

   /**
    * Returns the largest value taken by the variable since the last
    * initialization of this probe. This returns
    * <tt>Double.NEGATIVE_INFINITY</tt> if the probe was not updated since
    * the last initialization.
    *  @return the largest value taken by the collector since last
    * initialization
    */
   public double max() {
      return maxValue;
   }

   /**
    * Returns the sum cumulated so far for this probe. The meaning of this
    * sum depends on the subclass (e.g.,  @ref Tally or
    * @ref umontreal.ssj.simevents.Accumulate ). This returns 0 if the
    * probe was not updated since the last initialization.
    *  @return the sum for this probe
    */
   public double sum() {
      return sumValue;
   }

   /**
    * Returns the average for this collector. This returns `Double.NaN` if
    * the probe was not updated since the last initialization.
    *  @return the average value of the collected observations
    */
   abstract public double average();

   /**
    * Returns a string containing a report for this statistical collector.
    * The contents of this report depends on the statistical probe as well
    * as on the parameters set by the user through probe-specific methods.
    *  @return a report for this probe, represented as a string
    */
   abstract public String report();

   /**
    * Formats and returns a short, one-line report about this statistical
    * probe. This line is composed of whitespace-separated fields which
    * must correspond to the column names given by
    * {@link #shortReportHeader() shortReportHeader()}. This report should
    * not contain any end-of-line character, and does not include the name
    * of the probe. Its contents depends on the statistical probe as well
    * as on the parameters set by the user through probe-specific methods.
    *  @return the short report for the probe.
    */
   abstract public String shortReport();

   /**
    * Returns a string containing the name of the values returned in the
    * report strings. The returned string must depend on the type of probe
    * and on the reporting options only. It must not depend on the
    * observations received by the probe. This can be used as header when
    * printing several reports. For example, <tt>
    * <pre>
    *          System.out.println (probe1.shortReportHeader());
    *          System.out.println (probe1.getName() + " " + probe1.shortReport());
    *          System.out.println (probe2.getName() + " " + probe2.shortReport());
    *          ...
    * </pre>
    * </tt> Alternatively, one can use  {@link #report()
    * report(String,StatProbe[])} to get a report with aligned probe
    * names.
    *  @return the header string for the short reports.
    */
   abstract public String shortReportHeader();

   /**
    * Formats short reports for each statistical probe in the array
    * `probes` while aligning the probes’ names. This method first formats
    * the given global name. It then determines the maximum length
    * @f$\ell@f$ of the names of probes in the given array. The first
    * line of the report is composed of @f$\ell+3@f$ spaces followed by
    * the string returned by  #shortReportHeader called on the first probe
    * in `probes`. Each remaining line corresponds to a statistical probe;
    * it contains the probe’s name followed by the contents returned by
    * #shortReport. Note that this method assumes that `probes` contains
    * no `null` element.
    *  @param globalName   the global name given to the formatted report.
    *  @param probes       the probes to include in the report.
    *  @return the formatted report.
    */
   public static String report (String globalName, StatProbe[] probes) {
      int maxn = 0;
      StatProbe firstProbe = null;
      for (StatProbe probe : probes) {
         if (firstProbe == null)
            firstProbe = probe;
         String s = probe.getName();
         if (s != null && s.length() > maxn)
            maxn = s.length();
      }
      if (firstProbe == null)
         return "";
      StringBuffer sb = new StringBuffer ("Report for ");
      sb.append (globalName).append (PrintfFormat.NEWLINE);
      for (int i = 0; i < maxn; i++)
         sb.append (' ');
      sb.append ("   ");
      sb.append (firstProbe.shortReportHeader()).append (PrintfFormat.NEWLINE);
      for (StatProbe probe : probes) {
         sb.append
            (PrintfFormat.s (-maxn, probe.getName()));
         sb.append ("   ");
         sb.append (probe.shortReport()).append (PrintfFormat.NEWLINE);
      }
      return sb.toString();
   }

   /**
    * Equivalent to  #report(String,StatProbe[]), except that `probes` is
    * an  Iterable object instead of an array. Of course, the iterator
    * returned by `probes` should enumerate the statistical probes to
    * include in the report in a consistent and sensible order.
    *  @param globalName   the global name given to the formatted report.
    *  @param probes       the probes to include in the report.
    *  @return the formatted report.
    */
   public static String report (String globalName,
                                Iterable<? extends StatProbe> probes) {
      int maxn = 0;
      StatProbe firstProbe = null;
      for (StatProbe probe : probes) {
         if (firstProbe == null)
            firstProbe = probe;
         String s = probe.getName();
         int sl = s == null ? 4 : s.length();
         if (sl > maxn)
            maxn = sl;
      }
      if (firstProbe == null)
         return "";
      StringBuffer sb = new StringBuffer ("Report for ");
      sb.append (globalName).append (PrintfFormat.NEWLINE);
      for (int i = 0; i < maxn; i++)
         sb.append (' ');
      sb.append ("   ");
      sb.append (firstProbe.shortReportHeader()).append (PrintfFormat.NEWLINE);
      for (StatProbe probe : probes) {
         sb.append
            (PrintfFormat.s (-maxn, probe.getName()));
         sb.append ("   ");
         sb.append (probe.shortReport()).append (PrintfFormat.NEWLINE);
      }
      return sb.toString();
   }

   /**
    * Determines if this statistical probe is broadcasting observations to
    * registered observers. The default is `false`.
    *  @return the status of broadcasting.
    */
   public boolean isBroadcasting() {
      return broadcast;
   }

   /**
    * Instructs the probe to turn its broadcasting ON or OFF. The default
    * value is OFF. Warning: To avoid useless overhead and performance
    * degradation, broadcasting should never be turned ON when there are
    * no registered observers.
    *  @param b            `true` to turn broadcasting ON, `false` to turn
    *                      it OFF
    */
   public void setBroadcasting (boolean b) {
      broadcast = b;
   }

   /**
    * Determines if this statistical probe is collecting values. The
    * default is `true`.
    *  @return the status of statistical collecting.
    */
   public boolean isCollecting() {
      return collect;
   }

   /**
    * Turns ON or OFF the collection of statistical observations. The
    * default value is ON. When statistical collection is turned OFF,
    * observations added to the probe are passed to the registered
    * observers if broadcasting is turned ON, but are not counted as
    * observations by the probe itself.
    *  @param b            `true` to activate statistical collection,
    *                      `false` to deactivate it
    */
   public void setCollecting (boolean b) {
      collect = b;
   }

   /**
    * Adds the observation listener `l` to the list of observers of this
    * statistical probe.
    *  @param l            the new observation listener.
    *  @exception NullPointerException if `l` is `null`.
    */
   public void addObservationListener (ObservationListener l) {
      if (l == null)
         throw new NullPointerException();
      if (!listeners.contains (l))
         listeners.add (l);
   }

   /**
    * Removes the observation listener `l` from the list of observers of
    * this statistical probe.
    *  @param l            the observation listener to be deleted.
    */
   public void removeObservationListener (ObservationListener l) {
      listeners.remove (l);
   }

   /**
    * Removes all observation listeners from the list of observers of this
    * statistical probe.
    */
   public void clearObservationListeners() {
      listeners.clear();
   }

   /**
    * Notifies the observation `x` to all registered observers if
    * broadcasting is ON. Otherwise, does nothing.
    */
   public void notifyListeners (double x) {
      if (!broadcast)
         return;
      // We could also use the enhanced for loop here, but this is less efficient.
      final int nl = listeners.size();
      for (int i = 0; i < nl; i++)
         listeners.get (i).newObservation (this, x);
   }


   public StatProbe clone() throws CloneNotSupportedException {
      StatProbe s = (StatProbe)super.clone();
      s.listeners = new ArrayList<ObservationListener>(listeners);
      return s;
   }
}