package umontreal.ssj.mcqmctools.anova;

import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.hups.PointSetIterator;
import umontreal.ssj.mcqmctools.*;
import umontreal.ssj.stat.Tally;

import java.util.*;


/**
 * ANOVA variance estimator.
 *
 * Estimates the partial variances of multiple coordinate sets using the PartialVarianceEstimator
 * class.
 *
 */
public class AnovaVarianceEstimator implements MonteCarloModel<double[]> {

   // partial variance estimator
   protected PartialVarianceEstimator varEstimator;

   // partial variance integrator
   protected RandomIntegrator integrator;

   // variance storage
   double[] vars = null;
   

   public AnovaVarianceEstimator() {
      this.integrator = null;
      this.varEstimator = new PartialVarianceEstimator();
   }

   public MonteCarloModelDouble getModel() {
      return varEstimator.getModel();
   }
   
   public double getApproximateMean() {
      return varEstimator.getApproximateMean();
   }

   public void setModel (MonteCarloModelDoubleRQMC model, double approxMean) {
      varEstimator.setModel(model, approxMean);
   }

   public Integrator getIntegrator() {
      return integrator;
   }

   /**
    * Sets the integrator.
    * The integrator must provide twice the dimension of the model.
    *
    */
   public void setIntegrator (RandomIntegrator integrator) {
      this.integrator = integrator;
   }

   /**
    * Returns the list of coordinate sets under consideration.
    *
    */
   public List<CoordinateSet> getCoordinateSets() {
      return varEstimator.getCoordinateSets();
   }

   /**
    * Set the coordinate sets to consider to \c coordSets.
    *
    */
   public void setCoordinates (List<CoordinateSet> coordSets) {
      varEstimator.setCoordinateSets(coordSets);
   }

   /**
    * Set the coordinate sets to consider to all non-empty subsets of \c coords.
    *
    */
   public void setCoordinates (CoordinateSet coords) {
      varEstimator.setCoordinateSets(coords.subsetsNotEmpty());
   }

   /**
    * Set the coordinate sets to consider to all non-empty subsets of \c coords, up to
    * cardinality \c maxOrder.
    *
    */
   public void setCoordinates (CoordinateSet coords, int maxOrder) {
      varEstimator.setCoordinateSets(coords.subsetsNotEmpty(maxOrder));
   }

   
   /**
    * Simulates the model once.
    *
    * Returns \c vars such that:
    * \li \c vars.length is \c coordSets.size() + 2.
    * \li \c vars[\c nSets] contains the correction to the approximate mean.
    * \li \c vars[\c nSets+1] contains the total variance.
    *
    */
   @Override public void simulate (RandomStream stream) {

      List<CoordinateSet> coordSets = varEstimator.getCoordinateSets();
      int nSets = coordSets.size();

      if (vars == null || vars.length != nSets + 2)
         vars = new double[nSets + 2];

      if (integrator == null)
         throw new IllegalStateException("integrator has not been initialized");

      integrator.setStream(stream);
      integrator.integrate(varEstimator, vars);

      vars[nSets+1] -= vars[nSets] * vars[nSets];

      for (int i = 0; i < nSets; i++) {

         CoordinateSet cs = coordSets.get(i);

         // the components are ordered by coordinate mask
         // so there is no need to iterate over all of them
         for (int i2 = 0; i2 < i; i2++) {

            CoordinateSet cs2 = coordSets.get(i2);

            if (cs2.isSubset(cs))
               vars[i] -= vars[i2];
         }
      }
   }


   @Override public double[] getPerformance () {
	   return vars;
	   }
	   
	   
	   /**
    * Returns the number of dimensions for the input.
    *
    */
   public int getDimension() {
      return varEstimator.getDimension();
   }
   

   /**
    * Returns a description of the partial variance estimator.
    *
    */
   @Override public String toString() {
      String s = String.format("ANOVA Variance Estimator"
            + " [model=%s]", varEstimator.getModel());
      return s;
   }
   

}
