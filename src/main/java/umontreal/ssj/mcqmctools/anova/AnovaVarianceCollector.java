package umontreal.ssj.mcqmctools.anova;

import umontreal.ssj.stat.list.ListOfTallies;
import umontreal.ssj.stat.Tally;

import java.util.*;

/**
 * Extends ListOfTallies to collect ANOVA variances.
 *
 */
public class AnovaVarianceCollector extends ListOfTallies<PartialVariance> {

   protected List<AnovaObserver> observers;
   protected Tally meanCorrection;
   protected Tally totalVar;
   protected boolean sorted;

   /**
    * Constructs an ANOVA collector for the coordinate sets \c coordSets.
    *
    */
   public AnovaVarianceCollector(Iterable<CoordinateSet> coordSets) {

      sorted = false;
      totalVar = new Tally("Total variance");
      meanCorrection = new Tally("Correction to the mean");

      for (CoordinateSet cs : coordSets)
         add(new PartialVariance(cs, totalVar));

      observers = new ArrayList<AnovaObserver>();
   }

   /**
    * Add an observer, which is notified when the ANOVA variance collector is updated.
    *
    */
   public void addObserver(AnovaObserver observer) {
      observers.add(observer);
   }

   /**
    * Resets all collectors.
    *
    */
   @Override public void init() {
      super.init();
      meanCorrection.init();
      totalVar.init();
   }

   /**
    * Adds an observation of ANOVA variances.
    * The first \c size() items in \c vars are the observed partial variances.
    * The last two items are respectively the observed correction to the mean and total variance.
    *
    */
   @Override public void add(double[] vars) {
      if (sorted)
         throw new IllegalStateException("AnovaVarianceCollector cannot collect data once sorted");
      for (int i = 0; i < size(); i++)
         get(i).add(vars[i]);
      meanCorrection.add(vars[size()]);
      totalVar.add(vars[size() + 1]);

      // notify observers
      for (AnovaObserver obs : observers)
         obs.anovaUpdated(this);
   }

   
   /**
    * Sorts the computed variances.
    *
    * It is no longer possible to gather data after sorting.
    *
    */
   public void sort() {
      Collections.sort(this, Collections.reverseOrder());
      sorted = true;
   }

   /**
    * Returns the mean correction collector.
    *
    */
   public Tally getMeanCorrection() {
      return meanCorrection;
   }

   /**
    * Returns the total variance collector.
    *
    */
   public Tally getTotalVariance() {
      return totalVar;
   }

   /**
    * Returns the total variance for all projections of dimension \c order.
    *
    */
   public double getVarianceForOrder(int order) {
      double sum = 0;
      for (PartialVariance c : this)
         if (c.getCoordinates().cardinality() == order)
            sum += c.average();
      return sum;
   }

   /**
    * Returns the total variance fraction for all projections of dimension \c order.
    *
    */
   public double getVarianceFractionForOrder(int order) {
      double sum = 0;
      for (PartialVariance c : this)
         if (c.getCoordinates().cardinality() == order)
            sum += c.sensitivityIndex();
      return sum;
   }

   /**
    * Returns the total variance for a specific projection or a negative number if the
    * projection \c coords have not been analysed.
    *
    */
   public double getTotalVarianceForCoordinate(CoordinateSet coords) {
      for (PartialVariance c : this)
         if (c.getCoordinates().equals(coords))
            return c.average();
      return -1;
   }
   
   /**
    * Returns the total variance for all projections involving coordinate \c coord.
    *
    */
   public double getTotalVarianceForCoordinate(int coord) {
      double sum = 0;
      for (PartialVariance c : this)
         if (c.getCoordinates().contains(coord))
            sum += c.average();
      return sum;
   }
   
   /**
    * Returns the maximum dimension of the projections contained in the list.
    *
    */
   public int getMaxOrder() {
      int x = 0;
      for (PartialVariance c : this)
         x = Math.max(x, c.getCoordinates().cardinality());
      return x;
   }

   /**
    * Returns a description of the ANOVA collector.
    *
    */
   @Override public String toString() {
      String s = "ANOVA Collector";
      s += String.format(" [maxOrder=%d]", getMaxOrder());
      return s;
   }
   

   /**
    * Returns a report of the current state of estimation of the ANOVA variances.
    * Values are printed under the form: value ± standard deviation.
    */
   public String report() {

      String report = "";

      report += String.format("%30s: %9.4g", "Correction to mean", getMeanCorrection().average());
      if (getMeanCorrection().numberObs() > 1) {
         double dvar = Math.sqrt(getMeanCorrection().variance() / getMeanCorrection().numberObs());
         report += String.format(" ± %.2g", dvar);
      }

      report += umontreal.ssj.util.PrintfFormat.NEWLINE;
      report += String.format("%30s: %9.4g", "Total variance", getTotalVariance().average());
      if (getTotalVariance().numberObs() > 1) {
         double dvar = Math.sqrt(getTotalVariance().variance() / getTotalVariance().numberObs());
         report += String.format(" ± %.2g", dvar);
      }

      report += umontreal.ssj.util.PrintfFormat.NEWLINE;

      report += umontreal.ssj.util.PrintfFormat.NEWLINE;

      for (PartialVariance c : this)
         report += c + umontreal.ssj.util.PrintfFormat.NEWLINE;
      
      for (int i = 1; i <= getMaxOrder(); i++)
         report += umontreal.ssj.util.PrintfFormat.NEWLINE
               + String.format("%30s: %9.4g  (%.4g %%)",
               String.format("order %d", i),
               getVarianceForOrder(i),
               getVarianceFractionForOrder(i) * 100);

      return report;
   }

}
