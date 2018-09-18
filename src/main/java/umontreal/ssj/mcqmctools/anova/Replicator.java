package umontreal.ssj.mcqmctools.anova;

import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.list.ListOfTallies;
import umontreal.ssj.rng.*;
// import umontreal.ssj.hups.*;
import umontreal.ssj.mcqmctools.*;

/**
 * Replicator class.
 * Replicates mutliple samples of an integral computed by an internal random integrator.
 *
 */
public class Replicator implements RandomIntegrator {

   protected int nReplicates;
   protected RandomIntegrator integrator;

   // internal tallies, to be allocated only once, and only if needed
   protected Tally statValue = null;
   protected ListOfTallies<Tally> statValueList = null;
   
   /**
    * Creates a replicator that samples \c nReplicates replicates using \c integrator.
    * 
    */
   public Replicator(int nReplicates, RandomIntegrator integrator) {
      this.nReplicates = nReplicates;
      this.integrator = integrator;
   }

   /**
    * Returns the internal integrator.
    *
    */
   public RandomIntegrator getIntegrator() {
      return integrator;
   }

   /**
    * Returns the number of replicates.
    *
    */
   public int getNumReplicates() {
      return nReplicates;
   }
   
   /**
    * Returns the number of points (or simulations) per sample integral.
    *
    */
   public int getNumPoints() {
      return integrator.getTotalSimulations();
   }
   
   /**
    * Returns the total number of times the model is simulated per call to an integrate()
    * method: the number of replications multiplied by the number of points.
    *
    */
   @Override public int getTotalSimulations() {
      return getNumReplicates() * getNumPoints();
   }

   /**
    * Use \c stream as the source of randomness.
    *
    * Sets the random stream of the point set randomization to \c stream.
    *
    */
   @Override public void setStream(RandomStream stream) {
      integrator.setStream(stream);
   }

   /**
    * Returns the currently used random stream. May be \c null.
    *
    */
   public RandomStream getStream() {
      return integrator.getStream();
   }

   /**
    * Integrates a model by means of simulation.
    *
    * The output values are added to the statistical collector \c statValue.
    *
    */
   @Override public void integrate (MonteCarloModelDouble model, Tally statValue) {
	  statValue.init();
      for (int i = 0; i < nReplicates; i++)
         statValue.add(integrator.integrate(model));
   }

   /**
    * Shorthand to integrate without having to pass a Tally object.
    *
    */
   @Override public double integrate (MonteCarloModelDouble model) {
      if (statValue == null)
         this.statValue = new Tally();
      else
         statValue.init();
      integrate(model, statValue);
      return statValue.average();
   }

   /**
    * Integrates a model by means of simulation.
    *
    * The output values are added to the statistical collector \c stat.
    *
    */
   @Override public void integrate (MonteCarloModel<double[]> model, ListOfTallies<? extends Tally> stat) {
      double[] val = new double[stat.size()];
      ListOfTallies<Tally> innerStat = ListOfTallies.createWithTally(stat.size());

      for (int i = 0; i < nReplicates; i++) {
         innerStat.init();
         integrator.integrate(model, innerStat);
         innerStat.average(val);
         stat.add(val);
      }
   }

   /**
    * Shorthand to integrate without having to pass a ListOfTallies object.
    *
    */
   public void integrate (MonteCarloModel<double[]> model, double[] values) {
      if (statValueList == null || statValueList.size() != values.length)
         this.statValueList = ListOfTallies.createWithTally(values.length);
      else
         statValueList.init();
      integrate(model, statValueList);
      statValueList.average(values);
   }

   @Override public String toString() {
      return "Replicator"
         + " [replicates=" + nReplicates + "]"
         + " [integrator=" + integrator + "]";
   }
}
