package umontreal.ssj.mcqmctools.anova;

import umontreal.ssj.rng.RandomStream;

public interface RandomSampler extends Sampler {
   /**
    * Use \c stream as the source of randomness.
    *
    */
   public void setStream(RandomStream stream);

   /**
    * Returns the currently used random stream. May be \c null.
    *
    */
   public RandomStream getStream();
}
