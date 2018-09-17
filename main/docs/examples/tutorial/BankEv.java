package tutorial;
import umontreal.ssj.simevents.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;
import umontreal.ssj.stat.*;

public class BankEv {
   static final double minute = 1.0 / 60.0;
   int      nbTellers;         // Number of tellers.
   int      nbBusy;            // Number of tellers busy.
   int      nbWait;            // Queue length.
   int      nbServed;          // Number of customers served so far
   double   meanDelay;         // Mean time between arrivals.
   Event    nextArriv       = new Arrival();   // The next arrival.
   RandomStream  streamArr  = new MRG32k3a();   // Customer's arrivals
   ErlangGen genServ = new ErlangConvolutionGen (new MRG32k3a(), 2, 1.0/minute);
   RandomStream  streamTeller = new MRG32k3a(); // Number of tellers
   RandomStream  streamBalk   = new MRG32k3a(); // Balking decisions
   Tally statServed = new Tally ("Nb. served per day");
   Tally avWait     = new Tally ("Average wait per day (hours)");
   Accumulate wait  = new Accumulate ("cumulated wait for this day");

   Event e9h45 = new Event() {
      public void actions() {
         meanDelay = 2.0*minute;
         nextArriv.schedule
            (ExponentialGen.nextDouble (streamArr, 1.0/meanDelay));
      }
   };

   Event e10h = new Event() {
      public void actions() {
         double u = streamTeller.nextDouble();
         if (u >= 0.2) nbTellers = 3;
         else if (u < 0.05) nbTellers = 1;
         else nbTellers = 2;
         while (nbWait > 0 && nbBusy < nbTellers) {
            nbBusy++;  nbWait--;
            new Departure().schedule (genServ.nextDouble());
         }
         wait.update (nbWait);
      }
   };

   Event e11h = new Event() {
      public void actions() {
         nextArriv.reschedule ((nextArriv.time() - Sim.time())/2.0);
         meanDelay = minute;
      }
   };

   Event e14h = new Event() {
      public void actions() {
         nextArriv.reschedule ((nextArriv.time() - Sim.time())*2.0);
         meanDelay = 2.0*minute;
      }
   };

   Event e15h = new Event() {
      public void actions() { nextArriv.cancel(); }
   };

   private boolean balk() {
      return (nbWait > 9) ||
             (nbWait > 5 && (5.0*streamBalk.nextDouble() < nbWait-5));
   }

   class Arrival extends Event {
      public void actions() {
         nextArriv.schedule
            (ExponentialGen.nextDouble (streamArr, 1.0/meanDelay));
         if (nbBusy < nbTellers) {
            nbBusy++;
            new Departure().schedule (genServ.nextDouble());
         } else if (!balk())
            { nbWait++;  wait.update (nbWait); }
      }
   }

   class Departure extends  Event {
      public void actions() {
         nbServed++;
         if (nbWait > 0) {
            new Departure().schedule (genServ.nextDouble());
            nbWait--;   wait.update (nbWait);
         }
         else nbBusy--;
      }
   };

   public void simulOneDay() {
      Sim.init();       wait.init();
      nbTellers = 0;    nbBusy    = 0;
      nbWait    = 0;    nbServed  = 0;
      e9h45.schedule (9.75);
      e10h.schedule (10.0);
      e11h.schedule (11.0);
      e14h.schedule (14.0);
      e15h.schedule (15.0);
      Sim.start();
      statServed.add (nbServed);
      wait.update();
      avWait.add (wait.sum());
   }

   public void simulateDays (int numDays) {
      for (int i=1; i<=numDays; i++)  simulOneDay();
      System.out.println (statServed.report());
      System.out.println (avWait.report());
   }

   public static void main (String[] args) {
       new BankEv().simulateDays (100);
   }
}
