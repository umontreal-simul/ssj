package umontreal.ssj.mcqmctools.anova;

import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.hups.PointSetIterator;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.list.ListOfTallies;
import umontreal.ssj.mcqmctools.*;

public class MonteCarloIntegrator extends MonteCarloSampler implements RandomIntegrator {

   // internal tallies, to be allocated only once, and only if needed
   protected Tally statValue = null;
   protected ListOfTallies<Tally> statValueList = null;

   public MonteCarloIntegrator (int samples) {
      super(samples);
   }
   
   public MonteCarloIntegrator (int samples, RandomStream stream) {
      super(samples, stream);
   }
   
   /**
    * Returns the number of samples.
    */
   @Override public int getNumPoints() {
      return getNumSamples();
   }

   /** @copydoc #getNumPoints() */
   @Override public int getTotalSimulations() {
      return getNumSamples();
   }

   /** @copydoc Integrator::integrate(MonteCarloModelDouble, Tally) */
   @Override public void integrate (MonteCarloModelDouble model, Tally statValue) {
 
      boolean isPointSet = (stream instanceof PointSetIterator);
      
      for (int i = 0; i < nSamples; i++) {
         model.simulate(stream);
         statValue.add(model.getPerformance());
         if (isPointSet)
            ((PointSetIterator)stream).resetToNextPoint();
      }
   }
   
   /** @copydoc Integrator::integrate(MonteCarloModelDouble) */
   @Override public double integrate (MonteCarloModelDouble model) {
      if (statValue == null)
         this.statValue = new Tally();
      else
         statValue.init();
      integrate(model, statValue);
      return statValue.average();
   }

   /** @copydoc Integrator::integrate(MonteCarloModel<double[]>, ListOfTallies<? extends Tally>) */
   @Override public void integrate (MonteCarloModel<double[]> model, ListOfTallies<? extends Tally> statValue) {

      boolean isPointSet = (stream instanceof PointSetIterator);
      
      for (int i = 0; i < nSamples; i++) {
         model.simulate(stream);
         statValue.add(model.getPerformance());
         if (isPointSet)
            ((PointSetIterator)stream).resetToNextPoint();
      }
   }

   /** @copydoc Integrator::integrate(MonteCarloModel<double[]>, double[]) */
   @Override public void integrate (MonteCarloModel<double[]> model, double[] values) {
      if (statValueList == null || statValueList.size() != values.length)
         this.statValueList = ListOfTallies.createWithTally(values.length);
      else
         statValueList.init();
      integrate(model, statValueList);
      statValueList.average(values);
   }

   @Override public String toString() {
      String s = "Monte Carlo Integrator [samples=" + getNumSamples() + "]";
      if (getStream() != null)
         s += " [stream=" + getStream().getClass().getSimpleName() + "]";
      return s;
   }   
}
