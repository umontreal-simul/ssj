package umontreal.ssj.mcqmctools.anova;

// import umontreal.ssj.rng.RandomStream;
// import umontreal.ssj.hups.RQMCPointSet;
// import umontreal.ssj.mcqmctools.*;
// import umontreal.ssj.mcqmctools.anova.*;

/**
 * This class automates the process of replicating estimators of the ANOVA variances.
 *
 * For more flexibility, use the AnovaVarianceEstimator class.
 *
 */
public class Anova {

   protected int maxOrder = Integer.MAX_VALUE;
   protected int maxCoordinate = Integer.MAX_VALUE;

   protected RandomIntegrator innerIntegrator;
   protected Integrator outerIntegrator;

   public Anova (Integrator outerIntegrator, RandomIntegrator innerIntegrator) {
      this.outerIntegrator = outerIntegrator;
      this.innerIntegrator = innerIntegrator;
   }

   public Anova () {
      this(null, null);
   }
   
   /**
    * Sets the outer integrator from which provides a RandomStream to the randomization of the inner
    * integrator.
    * The number of points in the outer integrator corresponds to the number of replications of the
    * ANOVA variance estimators.
    *
    */
   public void setOuterIntegrator (Integrator integrator) {
      this.outerIntegrator = integrator;
   }

   /**
    * Sets the inner integrator which is used to generate one random estimation of the ANOVA variances
    * at a time.
    *
    */
   public void setInnerIntegrator (RandomIntegrator integrator) {
      this.innerIntegrator = integrator;
   }

   /**
    * Sets the maximum coordinate index to consider.
    *
    */
   public void setMaxCoordinate (int maxCoordinate) {
      this.maxCoordinate = maxCoordinate;
   }

   /**
    * Sets the maximum projection order to consider.
    *
    */
   public void setMaxOrder (int maxOrder) {
      this.maxOrder = maxOrder;
   }

   /**
    * Produces multiple replicates of the ANOVA variance estimators.
    *
    * Equivalent to estimate(model, approxMean, null)
    *
    */
   public AnovaVarianceCollector estimate (MonteCarloModelDoubleRQMC model, double approxMean) {
      return estimate(model, approxMean, null);
   }

   /**
    * Produces multiple replicates of the ANOVA variance estimators.
    *
    * The number of replicates is determined by the number of points in the outer integrator.
    * Best precision is achieved when the mean value of the model is close to \c approxMean.
    *
    */
   public AnovaVarianceCollector estimate (MonteCarloModelDoubleRQMC model, double approxMean, AnovaObserver observer) {

      int dimension = Math.min(maxCoordinate + 1, model.getDimension());

      AnovaVarianceEstimator varEstimator = new AnovaVarianceEstimator();
      varEstimator.setModel(model, approxMean);
      varEstimator.setCoordinates(CoordinateSetLong.allCoordinates(dimension), maxOrder);
      varEstimator.setIntegrator(innerIntegrator);

      AnovaVarianceCollector varCollector = new AnovaVarianceCollector(varEstimator.getCoordinateSets());
      if (observer != null)
         varCollector.addObserver(observer);

      outerIntegrator.integrate(varEstimator, varCollector);
      varCollector.sort();

      return varCollector;
   }

   @Override public String toString() {
      String s = "ANOVA";
      if (maxCoordinate < Integer.MAX_VALUE)
         s += " [maxCoordinate=" + maxCoordinate + "]";
      if (maxOrder < Integer.MAX_VALUE)
         s += " [maxOrder=" + maxOrder + "]";
      if (outerIntegrator != null) {
         s += umontreal.ssj.util.PrintfFormat.NEWLINE;
         s += "  Outer Integrator:";
         s += umontreal.ssj.util.PrintfFormat.NEWLINE;
         s += "    " + outerIntegrator;
      }
      if (innerIntegrator != null) {
         s += umontreal.ssj.util.PrintfFormat.NEWLINE;
         s += "  Inner Integrator:";
         s += umontreal.ssj.util.PrintfFormat.NEWLINE;
         s += "    " + innerIntegrator;
      }

      return s;
   }
}
