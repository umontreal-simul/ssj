package markovchainrqmc;

import umontreal.ssj.markovchainrqmc.*;
import umontreal.ssj.hups.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.*;
import umontreal.ssj.util.sort.*;

public class BrownianTest {
   Brownian brownian;
   
   public BrownianTest (double x0, double dt) {
      brownian = new Brownian(x0,dt);
      tests();
   }
   
   public void tests() {
      RandomStream stream = new MRG32k3a();
      
      //1- Print trajectory and performance for 20 steps of the chain.
      System.out.println("1- Print trajectory and performance");
      brownian.initialState ();
      System.out.println("step =  0,   position =  " + brownian.x);
      for (int step = 1; step < 21; ++step){
          brownian.nextStep (stream);
          System.out.printf("step = %2d,   position = %8.5f%n", step, brownian.x);
      }
      System.out.printf("Performance = %8.5f%n", brownian.getPerformance());
      
      //2- Monte Carlo, 100 replications, 2^12 trajectories, 20 steps 
      Tally performance = new Tally("Performance of Brownian Motion");
      Tally replicates = new Tally("Replicates of Brownian Motion");
      for(int i=0; i<100; ++i){
         brownian.simulRuns(4096,20,stream,performance);
         stream.resetNextSubstream();
         replicates.add(performance.average());
      }
      System.out.println("\n 2- MC, 100 reps, 2^12 trajetories, 20 steps ");
      System.out.println(replicates.report());
      
      //3- RQMC, 100 replications, 2^12 trajectories, 20 steps 
      PointSetRandomization rand = new RandomShift(stream);
      PointSet p = new SobolSequence(12,31,20);
      brownian.simulRQMC(p,100,20,rand,replicates);
      System.out.println("\n 3- RQMC, 100 reps, 2^12 trajectories, 20 steps");
      System.out.println(replicates.report());
      
      //4- Array-RQMC, 100 replications, 2^12 trajectories, 20 steps
      PointSet p2 = new SobolSequence(12,31,1);
      MultiDimSort sort = new OneDimSort(0);
      ArrayOfComparableChains array = 
                              new ArrayOfComparableChains(brownian, rand, sort);
      array.makeCopies(4096);
      array.simulReplicatesArrayRQMC (p2, rand, sort, 0, 100, 20, replicates);
      System.out.println("\n 4- Array-RQMC, 100 replications" + 
                                              ", 2^12 trajectories, 20 steps");
      System.out.println(replicates.report());
      
      //5- Array-RQMC Simulate and use 2^12 trajectories, 20 steps 
      array.makeCopies(4096);
      array.initialStates();
      for (int step = 1; step < 21; ++step){
          array.sortChains();
          array.simulOneStepArrayRQMC (p2);
          for(MarkovChainComparable mc: array.getChains()){
             //Do something with mc
          }
      }
   }

   public static void main (String[] args) {
      double x0 = 0.0;
      double dt = 0.1;       
      BrownianTest testBrownian = new BrownianTest(x0,dt);
   }
}
