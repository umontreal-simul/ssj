package markovchainrqmc;

import umontreal.ssj.markovchainrqmc.MarkovChainComparable;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.probdist.NormalDist;

class Brownian extends MarkovChainComparable  {
   final double x0;                   // Initial position.
   final double dt, sqrtDt;           // Time interval between observations.
   double x;                           // Position.              

   public Brownian (double x0, double dt) {
       this.x0 = x0;
       this.dt = dt;
       if (dt < 0)
          throw new IllegalArgumentException("dt must be positive");
       sqrtDt = Math.sqrt (dt);        // Just for faster computation
       stateDim = 1;                   // Dimension of state.
       initialState();
   }

   // Sets initial position
   public void initialState () {
      x = x0;                          
   }

   // Simulates the next step.
   public void nextStep (RandomStream stream) {
      x += sqrtDt * NormalDist.inverseF01 (stream.nextDouble());  
   }

   // Returns performance mesure.
   public double getPerformance () {
      return Math.abs(x-x0);
   }
   
   // Compares value of x between two chains.
   public int compareTo (MarkovChainComparable m, int i) {
      if (!(m instanceof Brownian))     
         throw new IllegalArgumentException("Can't compare a " +
           "Brownian Markov chain with other types of Markov chains."); 
      switch(i) {
         case 0: 
            double mx = ((Brownian)m).x;
            return (x>mx ? 1 : (x<mx ? -1 : 0));
         default: 
            throw new AssertionError("Invalid state index" );
      }
   }
}
