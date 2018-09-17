package umontreal.ssj.mcqmctools.anova;

/**
 * ANOVA observer.
 * An AnovaVarianceCollector object notifies its observers when it is updated.
 *
 */
public interface AnovaObserver {
   public void anovaUpdated(AnovaVarianceCollector anova);
}
