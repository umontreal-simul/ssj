package ift6561examples;
import umontreal.ssj.stat.TallyStore;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.charts.HistogramChart;


/**
 * Estimates the value of a callable bond.
 */
public class CallableBond {
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
    double[] callValues = {1.025, 1.020, 1.015, 1.010, 1.005, 1, 1, 1, 1, 1, principal};

    OrnsteinUhlenbeckWithIntegratedProcess ornUhl;
    double protectionPeriodValue;
    double expectedUncallableValue;


    public CallableBond(RandomStream randomStream) {
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


    public TallyStore simulate(int nSimulations) {
        TallyStore tally = new TallyStore("callable bond value");
        for (int iSim =0; iSim < nSimulations; iSim++) {
            double[] rates = ornUhl.generatePath();
            double[] expectedDiscounts = ornUhl.getExpectedFutureDiscount(postDecisionPaymentTimes);

            tally.add(priceCallableBondForOneSimulation(rates, expectedDiscounts));
        }
        return tally;
    }


    public double priceCallableBondForOneSimulation(double[] rates, double[] expectedDiscounts) {
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
        int nSimulations = 100000;
        TallyStore tally = new CallableBond(new MRG32k3a()).simulate(nSimulations);

        System.out.println(tally.reportAndCIStudent(0.95, 4));

        HistogramChart chart = new HistogramChart("Callable bond", "bond value", "frequency", tally);
        System.out.println(chart.toLatex(4.0, 4.0));
    }
}