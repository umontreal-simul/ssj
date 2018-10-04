package ift6561examples;

import umontreal.ssj.stat.TallyStore;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.stochprocess.*;
import umontreal.ssj.mcqmctools.*;
import umontreal.ssj.charts.HistogramChart;

/**
 * Estimates the value of a callable bond.  Model from the paper of Ben Ameur, Breton, et al.
 * 
 * NOTE: This program needs to be improved!
 */
public class CallableBond implements MonteCarloModelDouble {
    double coupon = 0.0425;
    double principal = 1.0;
    int nCoupons = 21;
    double[] couponTimes;
    int nDecisions = 11;
    int nProtectionCoupons = nCoupons - nDecisions;
    double[] postDecisionPaymentTimes;
    double dtDecision = 0.1666;
    double timeToFirstCoupon = 0.172;
    double r0 = 0.05;
    double alpha = 0.44178462;
    double beta = 0.098397028;
    double sigma = 0.13264223;

    double[] thresholdRates = 
        {-0.124, -0.116, -0.106, -0.095, -0.082, -0.066, -0.051, -0.033, -0.010, 0.021, Double.NEGATIVE_INFINITY};
    double[] callValues = 
    	{1.025, 1.020, 1.015, 1.010, 1.005, 1, 1, 1, 1, 1, principal};

    OrnsteinUhlenbeckWithIntegratedProcess ornUhl;
    double protectionPeriodValue;
    double expectedUncallableValue;
    double[] rates;
    double[] expectedDiscounts;


    public CallableBond (RandomStream randomStream) {
        couponTimes = new double[nCoupons + 1]; // includes time 0.0, where there is no coupon.
        couponTimes[0] = 0.0;
        for (int iTime = 1; iTime <= nCoupons; iTime++)
            couponTimes[iTime] = timeToFirstCoupon + (double) (iTime - 1);

        double[] decisionTimes = new double[nDecisions + 1]; // includes time 0.0, where there is no decision
        decisionTimes[0] = 0.0;
        for (int iTime = 1; iTime <= nDecisions; iTime++)
            decisionTimes[iTime] = timeToFirstCoupon + (double) nProtectionCoupons + (double) (iTime - 1) - dtDecision;

        ornUhl = new OrnsteinUhlenbeckWithIntegratedProcess(r0, alpha, beta, sigma, new NormalGen(randomStream));
        ornUhl.setObservationTimes(decisionTimes, nDecisions);
        
        double[] expectedCouponDiscounts = ornUhl.getTotalAnalyticDiscount(couponTimes);
        protectionPeriodValue = 0.0;
        for (int iTime = 1; iTime <= nProtectionCoupons; iTime++)
            protectionPeriodValue += coupon * expectedCouponDiscounts[iTime];
        System.out.println("Protection period expected value = " + protectionPeriodValue);

        postDecisionPaymentTimes = new double[nDecisions + 1]; // includes 0.0
        postDecisionPaymentTimes[0] = 0.0;
        for (int iTime = 1; iTime <= nDecisions; iTime++)
            postDecisionPaymentTimes[iTime] = couponTimes[nProtectionCoupons + iTime];
    }

    public void simulate(RandomStream stream) {
       rates = ornUhl.generatePath();
       expectedDiscounts = ornUhl.getExpectedFutureDiscount(postDecisionPaymentTimes);
    }

    public double getPerformance() {
    	return priceCallableBond (rates, expectedDiscounts);
     }

    public double priceCallableBond (double[] rates, double[] expectedDiscounts) {
        double callableValue = protectionPeriodValue;
        boolean didCallBack = false;
        for (int iTime = 1; iTime <= nDecisions; iTime++) {
            callableValue += coupon * expectedDiscounts[iTime];
            if (rates[iTime] < thresholdRates[iTime - 1]) {
                callableValue += callValues[iTime - 1] * expectedDiscounts[iTime];
                    didCallBack = true;
                    break;
            }
        }
        if (!didCallBack) 
            callableValue += principal * expectedDiscounts[nDecisions];
        return callableValue;
    }


    public double[] getThresholdRates() {
        return thresholdRates;
    }


    public static void main(String[] args) {
        int n = 100000;
		TallyStore statValue = new TallyStore ("Stats on bond payoffs");
        CallableBond bond = new CallableBond (new MRG32k3a());
		System.out.println (MonteCarloExperiment.simulateRunsDefaultReport 
				(bond, n, new MRG32k3a(), statValue));
        HistogramChart chart = new HistogramChart("Callable bond", "bond value", "frequency", statValue);
        chart.toLatexFile ("CallableBondHist.tex", 4.0, 4.0);
    }
}