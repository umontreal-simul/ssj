package ift6561examples;
import umontreal.ssj.stat.TallyStore;
import umontreal.ssj.stat.list.lincv.ListOfTalliesWithCV;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.rng.MRG32k3a;


/**
 * Estimates the value of a callable bond using the uncallable bond as a control variate.
 */
public class CallableBondCV extends CallableBond {
    double expectedUncallableValue;


    public CallableBondCV(RandomStream randomStream) {
        super(randomStream);

        expectedUncallableValue = 0.0;
        double[] expectedCouponDiscounts = ornUhl.getTotalAnalyticDiscount(couponTimes);
        for (int iTime = 1; iTime <= nCoupons; iTime++)
            expectedUncallableValue += coupon * expectedCouponDiscounts[iTime];
        expectedUncallableValue += principal * expectedCouponDiscounts[nCoupons];
        System.out.println("Uncallable bond value = " + expectedUncallableValue);
    }


    /**
     * @return a {@code ListOfTalliesWithCV} which has its beta parameter set from the simulations,
     * without trial runs.
     */
    public ListOfTalliesWithCV<TallyStore> simulateWithCV(int nSimulations) {
        ListOfTalliesWithCV<TallyStore> tallyWithCV = ListOfTalliesWithCV.createWithTallyStore(1, 1);
        tallyWithCV.setExpectedValue(0, expectedUncallableValue);
        for (int iSim =0; iSim < nSimulations; iSim++) {
            double[] rates = ornUhl.generatePath();
            double[] expectedDiscount = ornUhl.getExpectedFutureDiscount(postDecisionPaymentTimes);
            
            double callableValue = priceCallableBond (rates, expectedDiscount);
            
            double uncallableValue = protectionPeriodValue;
            for (int iTime = 1; iTime <= nDecisions; iTime++) 
                uncallableValue += coupon * expectedDiscount[iTime];
            uncallableValue += principal * expectedDiscount[nDecisions];

            tallyWithCV.add(callableValue, uncallableValue);
        }
        tallyWithCV.estimateBeta();
        return tallyWithCV;
    }


    public static void main(String[] args) {
        int nSimulations = 1000000;
        ListOfTalliesWithCV<TallyStore> tallyWithCV = new CallableBondCV(new MRG32k3a()).simulateWithCV(nSimulations);

        TallyStore tallyBond = tallyWithCV.get(0);
        tallyBond.setName("Callable bond (no CV)");
        System.out.println(tallyBond.reportAndCIStudent(0.95));
        TallyStore tallyUncallableBond = tallyWithCV.get(1);
        tallyUncallableBond.setName("Uncallable bond (the CV)");
        System.out.println(tallyUncallableBond.reportAndCIStudent(0.95));

        System.out.println("callable bond average with CV = " + tallyWithCV.averageWithCV(0));
        double[] stdDevWithCV = new double[1];
        tallyWithCV.standardDeviationWithCV(stdDevWithCV);
        System.out.println("std. dev. with CV = " + stdDevWithCV[0]);
        System.out.println("corresponding 95% CI half-width = " + (1.96 * stdDevWithCV[0] / Math.sqrt((double) nSimulations)));
        System.out.println("beta for CV = " + tallyWithCV.getBeta().get(0,0));
        System.out.println("\nvariance reduction from CV = " + (tallyBond.variance() / stdDevWithCV[0] / stdDevWithCV[0]));
    }
}