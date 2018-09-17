package ift6561examples;
import umontreal.ssj.stat.Tally;
import java.io.*;

public class CallCenterCRN extends CallCenter {
   Tally statQoS1      = new Tally ("stats on QoS for config.1");
   Tally statQoS2      = new Tally ("stats on QoS for config.2");
   Tally statQoS3      = new Tally ("stats on QoS for config.3");
   Tally statLR1       = new Tally ("stats on loss rate for config.1");
   Tally statLR2       = new Tally ("stats on loss rate for config.2");
   Tally statLR3       = new Tally ("stats on loss rate for config.3");
   Tally statWait1      = new Tally ("stats on waiting time for config.1");
   Tally statWait2      = new Tally ("stats on waiting time for config.2");
   Tally statWait3      = new Tally ("stats on waiting time for config.3");
   Tally statDiffQoS21 = new Tally
      ("stats on difference between config.2 and config.1 for QoS");
   Tally statDiffQoS31 = new Tally
      ("stats on difference between config.3 and config.1 for QoS");
   Tally statDiffLR21 = new Tally
      ("stats on difference between config.2 and config.1 for loss rate");
   Tally statDiffLR31 = new Tally
      ("stats on difference between config.3 and config.1 for loss rate");
   Tally statDiffWait21 = new Tally
      ("stats on difference between config.2 and config.1 for waiting time");
   Tally statDiffWait31 = new Tally
      ("stats on difference between config.3 and config.1 for waiting time");

   double p1, p2, p3;
   double nu1, nu2, nu3;

   public CallCenterCRN (String fileName, double p1, double nu1,
                         double p2, double nu2,
                         double p3, double nu3) throws IOException {
      super (fileName);
      this.p1 = p1;     this.p2 = p2;      this.p3 = p3;
      this.nu1 = nu1;     this.nu2 = nu2;     this.nu3 = nu3;
   }

   // Overridden to support the case of no abandonment
   public double generPatience() {
      if (nu <= 0)
         // ExponentialDist.inverseF crashes with nu = 0.
         return Double.POSITIVE_INFINITY;
      else
         return super.generPatience();
   }

   public void setPatience (double p, double nu) {
      this.p = p;     this.nu = nu;
   }

   public void simulateDiffCRN (int n) {
      double qos1, qos2, qos3,  lr1, lr2, lr3, w1, w2, w3;
      statQoS1.init();     statQoS2.init();     statQoS3.init();
      statLR1.init();      statLR2.init();      statLR3.init();
      statWait1.init();    statWait2.init();    statWait3.init();
      statDiffQoS21.init();   statDiffQoS31.init();
      statDiffLR21.init();     statDiffLR31.init();
      statDiffWait21.init();   statDiffWait31.init();
      for (int i=0; i<n; i++) {
         setPatience (p1, nu1);
         streamB.resetNextSubstream();
         streamArr.resetNextSubstream();
         streamPatience.resetNextSubstream();
         (genServ.getStream()).resetNextSubstream();
         simulateOneDay();
         qos1 = (double)nGoodQoS / nCallsExpected;
         lr1 = (double)nAbandon / nCallsExpected;
         w1 = statWaitsDay.sum() / nCallsExpected;

         setPatience (p2, nu2);
         streamB.resetStartSubstream();
         streamArr.resetStartSubstream();
         streamPatience.resetStartSubstream();
         (genServ.getStream()).resetStartSubstream();
         waitList.clear();
         simulateOneDay();
         qos2 = (double)nGoodQoS / nCallsExpected;
         lr2 = (double)nAbandon / nCallsExpected;
         w2 = statWaitsDay.sum() / nCallsExpected;

         setPatience (p3, nu3);
         streamB.resetStartSubstream();
         streamArr.resetStartSubstream();
         streamPatience.resetStartSubstream();
         (genServ.getStream()).resetStartSubstream();
         waitList.clear();
         simulateOneDay();
         qos3 = (double)nGoodQoS / nCallsExpected;
         lr3 = (double)nAbandon / nCallsExpected;
         w3 = statWaitsDay.sum() / nCallsExpected;

         statQoS1.add (qos1);   statQoS2.add (qos2);   statQoS3.add (qos3);
         statLR1.add (lr1);     statLR2.add (lr2);     statLR3.add (lr3);
         statWait1.add (w1);    statWait2.add (w2);    statWait3.add (w3);
         statDiffQoS21.add (qos2 - qos1);   statDiffQoS31.add (qos3 - qos1);
         statDiffLR21.add (lr2 - lr1);   statDiffLR31.add (lr3 - lr1);
         statDiffWait21.add (w2 - w1);   statDiffWait31.add (w3 - w1);
      }
   }
            
   static public void main (String[] args) throws IOException { 
      int n = 1000;   // Number of replications.

      CallCenterCRN cc = new CallCenterCRN
         ("CallCenter.dat", 0.1, 0.001,  0, 0,  0.2, 0.04); 

      cc.simulateDiffCRN (n);
      System.out.println (
          cc.statQoS1.reportAndCIStudent (0.9) +
          cc.statQoS2.reportAndCIStudent (0.9) +
          cc.statQoS3.reportAndCIStudent (0.9) +
          cc.statDiffQoS21.reportAndCIStudent (0.9) +
          cc.statDiffQoS31.reportAndCIStudent (0.9) +
          cc.statLR1.reportAndCIStudent (0.9) +
          cc.statLR2.reportAndCIStudent (0.9) +
          cc.statLR3.reportAndCIStudent (0.9) +
          cc.statDiffLR21.reportAndCIStudent (0.9) +
          cc.statDiffLR31.reportAndCIStudent (0.9) +
          cc.statWait1.reportAndCIStudent (0.9) +
          cc.statWait2.reportAndCIStudent (0.9) +
          cc.statWait3.reportAndCIStudent (0.9) +
          cc.statDiffWait21.reportAndCIStudent (0.9) +
          cc.statDiffWait31.reportAndCIStudent (0.9));
   }
}
