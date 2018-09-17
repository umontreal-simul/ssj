package umontreal.ssj.mcqmctools.anova;

import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.Chrono;

public class MeanVarExperiment {
   protected String name;
   protected MonteCarloModelDouble model;
   protected RandomIntegrator integrator;
   protected double cpuSeconds;
   protected double average;
   protected double variance;
   protected int nObservations;
   protected String statReport;
   
   public MeanVarExperiment(String name, MonteCarloModelDouble model, RandomIntegrator integrator) {
      this.name = name;
      this.model = model;
      this.integrator = integrator;
      this.cpuSeconds = 0;
      this.average = 0;
      this.variance = 0;
      this.nObservations = 0;
   }

   public MeanVarExperiment(MonteCarloModelDouble model, RandomIntegrator integrator) {
      this(null, model, integrator);
   }

   public String getName() {
      return name;
   }

   public void setName(String name) {
      this.name = name;
   }

   public MonteCarloModelDouble getModel() {
      return model;
   }
   
   public void setModel(MonteCarloModelDouble model) {
      this.model = model;
   }
   
   public Integrator getIntegrator() {
      return integrator;
   }
   
   public void setIntegrator(RandomIntegrator integrator) {
      this.integrator = integrator;
   }
   
   public double getCPUSeconds() {
      return cpuSeconds;
   }
   
   public double getSecondsPerSimulation() {
      return cpuSeconds / integrator.getTotalSimulations();
   }

   public double getAverage() {
      return average;
   }

   public double getVariance() {
      return variance;
   }
   
   /**
    * Returns the variance scaled with the number of function evaluations per
    * observation.
    *
    */
   public double getScaledVariance() {
      return variance * integrator.getTotalSimulations() / nObservations;
   }
   
   public void simulate() {
      simulate(new Tally(model.getClass().getSimpleName()));
   }

   public void simulate(Tally stat) {
      Chrono timer = new Chrono();
      stat.init();
      integrator.integrate(model, stat);
      cpuSeconds = timer.getSeconds();
      average = stat.average();
      variance = stat.variance();
      nObservations = stat.numberObs();
      stat.setConfidenceIntervalStudent();
      statReport = stat.report(0.95, 4);
   }
   
   public String report() {
      return statReport
         + String.format("Variance:     %9.4g\n", getVariance())
         + String.format("Scaled Var.:  %9.4g\n", getScaledVariance())
         + String.format("CPU time:     %9s\n", Chrono.format(cpuSeconds));
   }
   
   public String toString() {
      return "MeanVarExperiment"  + (name == null ? "" : " " + name) + ":\n" +
            "  Model: " + model + "\n" +
            "  Integrator: " + integrator;
   }
}
