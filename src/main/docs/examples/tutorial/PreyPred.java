package tutorial;
import umontreal.ssj.simevents.*;

public class PreyPred {
   double r  = 0.005,      c  = 0.00001,
          s  = 0.01,       d  = 0.000005,     h = 5.0;
   double x0 = 2000.0,     z0 = 150.0;
   double horizon = 501.0;
   Simulator sim = new Simulator();
   Continuous x;
   Continuous z;

   public static void main (String[] args) { new PreyPred(); }

   public PreyPred() {
      x = new Preys(sim);
      z = new Preds(sim);
      sim.init();
      new EndOfSim(sim).schedule (horizon);
      new PrintPoint(sim).schedule (h);
      (sim.continuousState()).selectRungeKutta4 (h);
      x.startInteg (x0);
      z.startInteg (z0);
      sim.start();
   }

   public class Preys extends Continuous {
      public Preys(Simulator sim) { super(sim); }

      public double derivative (double t) {
         return (r * value() - c * value() * z.value());
      }
   }

   public class Preds extends Continuous {
      public Preds(Simulator sim) { super(sim); }

      public double derivative (double t) {
         return (-s * value() + d * x.value() * value());
      }
   }

   class PrintPoint extends Event {
      public PrintPoint(Simulator sim) { super(sim); }
      public void actions() {
         System.out.println (sim.time() + "  " + x.value() + "  " + z.value());
         this.schedule (h);
      }
   }

   class EndOfSim extends Event {
      public EndOfSim(Simulator sim) { super(sim); }
      public void actions() { sim.stop(); }
   }
}
