package ift6561examples;

import umontreal.ssj.stat.TallyStore;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.charts.XYLineChart;
import umontreal.ssj.mcqmctools.MonteCarloExperiment;

/**
 * Estimates the value of a callable bond where the coupons are paid at
 * the same time as in {@code CallableBond}, but where there is a single
 * call back opportunity at 10.0054 years (the 11th coupon).
 */
public class CallableBondCMC extends CallableBond {
    /** The number of coupons after the last decision (excluding the coupon right after that decision). */
    int nPostDecisionCoupons; 
    /** All coupon times after the decision (includes time 0.0 and the coupon right after the decision). */
    double[] postDecisionPaymentTimes; 
    /** As above, but in the case when the bond is called back, such that only the times right after 
        the last decision are needed. */
    double[] postDecisionPaymentTimesCalled; 


    public CallableBondCMC(RandomStream randomStream) {
        super(randomStream);

        nDecisions = 1;
        nPostDecisionCoupons = nCoupons - nProtectionCoupons - nDecisions;
        thresholdRates = null; // not needed for this class
        callValues = new double[]{1.025}; 

        double[] decisionTimes = new double[nDecisions + 1]; // includes time 0.0, where there is no decision
        decisionTimes[0] = 0.0;
        for (int iTime = 1; iTime <= nDecisions; iTime++)
            decisionTimes[iTime] = timeToFirstCoupon + (double) nProtectionCoupons + (double) (iTime - 1) - dtDecision;

        ornUhl = new OrnsteinUhlenbeckWithIntegratedProcess(r0, alpha, beta, sigma, new NormalGen(randomStream));
        ornUhl.setObservationTimes(decisionTimes, nDecisions);
        
        postDecisionPaymentTimes = new double[nDecisions + 1 + nPostDecisionCoupons]; // includes time 0.0
        postDecisionPaymentTimes[0] = 0.0;
        for (int iTime = 1; iTime <= (nDecisions + nPostDecisionCoupons); iTime++)
            postDecisionPaymentTimes[iTime] = couponTimes[nProtectionCoupons + iTime];
        postDecisionPaymentTimesCalled = new double[]{0.0, couponTimes[nProtectionCoupons + 1]};
    }


    /** {@code resetNextSubstream} is called after each simulation such that common random numbers can be used. */
    public TallyStore simulateCMC(int nSimulations, double thresholdRate) {
        TallyStore tally = new TallyStore("callable bond value with threshold " + thresholdRate);
        for (int iSim =0; iSim < nSimulations; iSim++) {
            double bondValue = protectionPeriodValue;
            double[] rates = ornUhl.generatePath();
            if (rates[1] < thresholdRate) { // call back the bond
                double[] expectedDiscounts = ornUhl.getExpectedFutureDiscount(postDecisionPaymentTimesCalled);
                bondValue += expectedDiscounts[1] * (callValues[0] + coupon);
            }
            else { // all subsequent coupons and the principal are paid
                double[] expectedDiscounts = ornUhl.getExpectedFutureDiscount(postDecisionPaymentTimes);
                for (int iTime = 1; iTime <= (nDecisions + nPostDecisionCoupons); iTime++) {
                    bondValue += expectedDiscounts[iTime] * coupon;
                }
                bondValue += expectedDiscounts[nDecisions + nPostDecisionCoupons] * principal;
            }

            tally.add(bondValue);
            ornUhl.getStream().resetNextSubstream();
        }
        return tally;
    }


    public void evaluateBondAtVariousThresholds(double[] thresholdRateTrials, int nSimulations, 
                                                boolean doUseCommonRandomNumbers) {
        System.out.println("using common random numbers: " + doUseCommonRandomNumbers);
        double[] results = new double[thresholdRateTrials.length];
        for (int iTrial = 0; iTrial < thresholdRateTrials.length; iTrial++) {
            if (doUseCommonRandomNumbers)
                ornUhl.getStream().resetStartStream();
            TallyStore tallyTrial = simulateCMC(nSimulations, thresholdRateTrials[iTrial]);
            System.out.println(tallyTrial.reportAndCIStudent(0.95, 4));
            results[iTrial] = tallyTrial.average();
        }
        XYLineChart chart = new XYLineChart(doUseCommonRandomNumbers ? "CRN" : "IRN", "threshold rate", "bond value", 
                                            new double[][]{thresholdRateTrials, results});
        chart.view(500, 500);
        System.out.println(chart.toLatex(4.0, 4.0));
    }

    
    public static void main(String[] args) {
        int n = 100000;

        RandomStream randomStream = new MRG32k3a();
        CallableBond bond = new CallableBond(randomStream);
        for (int iThreshold = 1; iThreshold < bond.getThresholdRates().length; iThreshold++)
            bond.getThresholdRates()[iThreshold] = Double.NEGATIVE_INFINITY;
		TallyStore statValue = new TallyStore ("Callable bond value with single threshold (-0.124) and no CMC");
		System.out.println (MonteCarloExperiment.simulateRunsDefaultReport 
				(bond, n, new MRG32k3a(), statValue));

        CallableBondCMC bondCMC = new CallableBondCMC(randomStream);
        double thresholdRate = -0.124;
        TallyStore tallyCMC = bondCMC.simulateCMC(n, thresholdRate);
        System.out.println(tallyCMC.reportAndCIStudent(0.95, 4));
        System.out.printf("%n*** variance reduction factor from CMC = %.1f%n%n%n", (statValue.variance() / tallyCMC.variance()));

        boolean doUseCommonRandomNumbers = true;
        double[] thresholdRateTrials = {-0.20, -0.12, -0.08, -0.06, -0.05, -0.04, -0.03, -0.02, 0.00, 0.08};

        bondCMC.evaluateBondAtVariousThresholds(thresholdRateTrials, n, doUseCommonRandomNumbers);

        doUseCommonRandomNumbers = false;
        bondCMC.evaluateBondAtVariousThresholds(thresholdRateTrials, n, doUseCommonRandomNumbers);
    }
}