package umontreal.ssj.mcqmctools.anova;

import umontreal.ssj.rng.RandomStream;
// import umontreal.ssj.hups.PointSetIterator;
import umontreal.ssj.mcqmctools.*;
// import umontreal.ssj.stat.Tally;

import java.util.*;


/**
 * Partial variance estimator.
 *
 * Estimates partial variances of a model with respect to multiple coordinate sets.
 *
 * Reference:
 * Monte Carlo estimators for small sensitivity indices
 * I. M. Sobol' and E. E. Myshetskaya
 * Monte Carlo Methods Appl. Vol. 13 No. 5-6 (2007), pp. 455-465
 *
 */
public class PartialVarianceEstimator implements MonteCarloModel<double[]> {

   // model whose variance is to be estimated
   protected MonteCarloModelDoubleRQMC model;

   // approximation to the average
   // the better it is, the better the variance estimator is
   protected double approxMean;

   // coordinate sets to be considered
   protected List<CoordinateSet> coordSets;
   protected CoordinateSet noCoordinate;
   protected CoordinateSet allCoordinates;

   
   // placeholder for the computed partial variances
   protected double[] vars;
   
   public PartialVarianceEstimator() {
      this.model = null;
      this.coordSets = null;
      this.approxMean = 0;
      this.vars = null;
   }

   public PartialVarianceEstimator (MonteCarloModelDoubleRQMC model, double approxMean,
         List<CoordinateSet> coordSets) {
      setModel(model, approxMean);
      setCoordinateSets(coordSets);
      this.vars = null;
   }

   public MonteCarloModelDoubleRQMC getModel() {
      return model;
   }

   /**
    * Sets the model whose partial variances are to be estimated.
    *
    * Best precision is achieved when the mean value of the model is close to \c approxMean.
    *
    */
   public void setModel (MonteCarloModelDoubleRQMC model, double approxMean) {
      this.model = model;
      this.approxMean = approxMean;
   }

   /**
    * Returns the list of coordinate sets under consideration.
    *
    */
   public List<CoordinateSet> getCoordinateSets() {
      return coordSets;
   }

   /**
    * Set the coordinate sets to consider to \c coordSets.
    *
    */
   public void setCoordinateSets (List<CoordinateSet> coordSets) {
      this.coordSets = coordSets;
      noCoordinate = new CoordinateSetLong(0);
      allCoordinates = new CoordinateSetLong(-1);
   }

   public double getApproximateMean() {
      return approxMean;
   }

   

   /**
    * Simulates the estimator once.
    *
    * Returns \c vars such that:
    * \li \c vars.length is \c coordSets.size() + 2.
    * \li \c vars[\c nSets] contains the correction to the approximate mean.
    * \li \c vars[\c nSets+1] contains the square of the above correction.
    *
    */
   @Override public void simulate (RandomStream stream) {

      // initialize storage
      if (vars == null || vars.length != coordSets.size() + 2)
         vars = new double[coordSets.size() + 2];

      
      if (model == null)
         throw new IllegalArgumentException("model has not been initialized");

      if (coordSets == null)
         throw new IllegalArgumentException("the coordinate sets have not been initialized");

      int nSets = coordSets.size();

      if (vars.length < nSets + 2)
         throw new IllegalArgumentException("vars[] must contain one more element than the"
               + " number of coordinate sets");


      SplitStream s = new SplitStream(stream, model.getDimension());

      s.setCoordinates(allCoordinates);
      s.resetStartSubstream();
      model.simulate(s);
      double valAll = model.getPerformance() - approxMean;

      s.setCoordinates(noCoordinate);
      s.resetStartSubstream();
      model.simulate(s);
      double valNone = model.getPerformance() - approxMean;
      
      // correction to the approximate mean
      vars[nSets] = valAll;

      // square correction
      vars[nSets+1] = valAll * valAll;

      for (int j = 0; j < nSets; j++) {

         // FIXME: check if the current coordinate set contains all coordinates
         // and reuse var[nSets] if applicable.

         s.setCoordinates(coordSets.get(j));
         s.resetStartSubstream();
         model.simulate(s);
         double valPartial = model.getPerformance() - approxMean;

         vars[j] = valAll * (valPartial - valNone);
      }
   }
   

   @Override public double[] getPerformance () {
	   return vars;
	   }
	   

   /**
    * Returns the number of input dimensions.
    *
    */
   public int getDimension() {
      return (model == null) ? 0 : 2 * model.getDimension();
   }
   

   /**
    * Returns a description of the partial variance estimator.
    *
    */
   @Override public String toString() {
      String s = String.format("Partial Variance Estimator"
            + " [model=%s]", model.toString());
      return s;
   }
   

}
