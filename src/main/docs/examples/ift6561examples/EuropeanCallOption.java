package ift6561examples;
// package umontreal.ssj.finance;

import umontreal.ssj.stat.Tally;
import umontreal.ssj.stochprocess.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.mcqmctools.*;

/**
 * This class represents an <SPAN CLASS="textit">European call</SPAN>
 * option with European exercise type. The payoff of this option at the time of
 * expiration is given by the formula
 * 
 * <P>
 * </P>
 * <DIV ALIGN="CENTER" CLASS="mathdisplay"><A NAME="as:AsianPayoff"></A> payoff
 * = max(0, bar(S)<SUB>T</SUB> - <I>K</I>) </DIV>
 * <P>
 * </P>
 * where <SPAN CLASS="MATH"><I>K</I></SPAN> is the strike price, and
 * <SPAN CLASS="MATH">bar(S)<SUB>T</SUB></SPAN> is the
 * <SPAN CLASS="textit">arithmetic</SPAN> average
 * 
 * <P>
 * </P>
 * <DIV ALIGN="CENTER" CLASS="mathdisplay"> bar(S)<SUB>T</SUB> =
 * <IMG ALIGN="MIDDLE" BORDER="0" SRC="AsianOptionimg1.png" ALT="$\displaystyle
 * {\frac{1}{n}}$">&sum;<SUB>i=1</SUB><SUP>n</SUP><I>S</I>(<I>t</I><SUB>i</SUB>)
 * </DIV>
 * <P>
 * </P>
 * of the option's underlying asset price at the observation times
 * <SPAN CLASS="MATH"><I>t</I><SUB>i</SUB></SPAN> ( <SPAN CLASS="MATH"><I>i</I>
 * = 1,&#8230;, <I>n</I></SPAN>), with <SPAN CLASS="MATH"><I>t</I><SUB>n</SUB> =
 * <I>T</I></SPAN>, the time of expiration of the option.
 * 
 * <P>
 * The initial value <SPAN CLASS="MATH"><I>S</I>(<I>t</I><SUB>0</SUB>)
 * = <IMG ALIGN="BOTTOM" BORDER="0" SRC="AsianOptionimg2.png" ALT="$ \tt
 * s0$"></SPAN> of the price process is not included in the calculation of the
 * average. However it will be included if the user sets the first observation
 * time <SPAN CLASS="MATH"><I>t</I><SUB>1</SUB></SPAN> to be the same as the
 * initial time <SPAN CLASS="MATH"><I>t</I><SUB>0</SUB></SPAN>.
 * 
 */
public class EuropeanCallOption implements MonteCarloModelDouble {

   StochasticProcess priceProcess; // Underlying process for the price.
   double r;          // Interest rate.
   double expirDate;  // Expoiration date. 
   double[] path;     // Sample path of the process.
   double strike;     // Strike price.
   double discount;   // Discount factor exp(-r * obsTimes[d]).

   /**
    * Array <TT>obsTimes[0..d+1]</TT> must contain <TT>obsTimes[0] = 0</TT>, plus
    * the <SPAN CLASS="MATH"><I>d</I></SPAN> observation times.
    * 
    */
   public EuropeanCallOption(double r, double expirDate, double strike) {
      this.r = r;
      this.expirDate = expirDate;
      this.strike = strike;
      discount = Math.exp(-r * expirDate);
   }

   // This constructor also specifies the underlying process.
   public EuropeanCallOption(StochasticProcess sp, double r, double expirDate, double strike) {
      this(r, expirDate, strike);
      setProcess(sp);
   }

   /**
    * Reset the process to <TT>sp</TT>. Assumes that <TT>obsTimes</TT> have been
    * set.
    */
   public void setProcess(StochasticProcess sp) {
      // Reset the process to sp.
      priceProcess = sp;
      sp.setObservationTimes(expirDate, 1);
   }

   /**
    * Computes and returns discounted payoff. Assumes path has been generated.
    */
   public double getPerformance() {
      if (path[1] > strike)
         return discount * (path[1] - strike);
      else
         return 0.0;
   }

   /**
    * Returns the number of observation times <SPAN CLASS="MATH"><I>d</I></SPAN>.
    * 
    */
   public int getNumObsTimes() {
      return 1;
   }

   /**
    * Generate a sample path of the process using <TT>stream</TT>
    */
   public void simulate(RandomStream stream) {
      path = priceProcess.generatePath(stream);
      // Note: We cannot pre-generate RQMC points here and call
      // generatePath(points), because this is not defined for all process types.
   }

   /**
    * Performs <SPAN CLASS="MATH"><I>n</I></SPAN> independent runs using
    * <TT>stream</TT> and collects statistics in <TT>statValue</TT>. The collector
    * <TT>statValue</TT> collects only the positive payoffs.
    */
   public void simulateRuns(int n, RandomStream stream, Tally statValue, Tally statValuePos) {
      statValue.init();
      statValuePos.init();
      double x;
      for (int i = 0; i < n; i++) {
         simulate(stream);
         x = getPerformance();
         statValue.add(x);
         if (x > 0.0000000001)
            statValuePos.add(x);
         stream.resetNextSubstream();
      }
   }

   public String toString() {
      return "European call option model";
   }
   
   public String getTag() {
      return "EuropeanCallOption";
   }

}
