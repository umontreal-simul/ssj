package umontreal.ssj.mcqmctools.anova;

import umontreal.ssj.mcqmctools.MonteCarloModel;
import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.list.ListOfTallies;

public interface Integrator {

   /**
    * Returns the number of points per integral.
    *
    */
   public int getNumPoints();

   /**
    * Returns the total number of times the model is simulated per call to an integrate()
    * method.
    *
    */
   public int getTotalSimulations();

   /**
    * Integrates a model by means of simulation.
    *
    * The output values are added to the statistical collector \c statValue.
    *
    */
   public void integrate (MonteCarloModelDouble model, Tally statValue);

   /**
    * Shorthand to integrate without having to pass a Tally object.
    *
    */
   public double integrate (MonteCarloModelDouble model);

   /**
    * Integrates a model by means of simulation.
    *
    * The output values are added to the statistical collector \c statValue.
    *
    */
   public void integrate (MonteCarloModel<double[]> model, ListOfTallies<? extends Tally> stat);

   /**
    * Shorthand to integrate without having to pass a ListOfTallies object.
    *
    */
   public void integrate (MonteCarloModel<double[]> model, double[] values);
}
