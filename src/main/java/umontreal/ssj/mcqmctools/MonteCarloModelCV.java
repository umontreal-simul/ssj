package umontreal.ssj.mcqmctools;

/**
 * An extension of @ref MonteCarloModelDouble that also implements a vector of
 * control variates.
 */

public interface MonteCarloModelCV extends MonteCarloModelDouble {

   /**
    * Recovers the realizations of the control variates for the last run.
    */
   public double[] getValuesCV();

   /**
    * TO DO: Recovers the realizations of the control variates from the last run.
    */
   // public void getValuesCV (double[] valCV);

   /**
    * Returns the number of control variates.
    */
   public int getNumberCV();

}
