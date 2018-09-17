package ift6561examples;

import umontreal.ssj.randvar.*;
import umontreal.ssj.randvarmulti.*;
import umontreal.ssj.stochprocess.OrnsteinUhlenbeckProcess;
import cern.colt.matrix.impl.*;

/** 
 * This class is custom-made for the `CallableBond` example.  
 * It generates an Ornstein-Uhlenbeck process as in @ref umontreal.ssj.stochprocess.OrnsteinUhlenbeckProcess,
 * but it also computes the integral of the process,  $\int_{t_0}^{t} X(\tau) d\tau$. 
 * The 2 x 2 correlation matrix for the process R(t) and its integral is decomposed with PCA.
 */
public class OrnsteinUhlenbeckWithIntegratedProcess extends OrnsteinUhlenbeckProcess {
    double[] integratedPath; // length [d + 1]
    MultinormalGen[] normalsCorrGen; // length [d]
    double[] bdt; // length [d]:  (beta * dt)
    double[] oneMExpMADtOverAlpha; // length [d]:  (1.0 - exp(-alpha * dt)) / alpha
    double c22; // sigma * sigma / 2.0 / alpha^3


    /**
     * @param x0 The starting value of the process.
     * @param alpha The force pushing the process back to its average value, b.
     * @param b The average value of the process.    
     * @param sigma The volatility.
     */
    public OrnsteinUhlenbeckWithIntegratedProcess(double x0, double alpha, double b, double sigma, NormalGen gen)  {
        super(x0, alpha, b, sigma, gen);
    }


    /**
     * Sets up various constants and decomposes the correlation matrix at each time step.
     */
    @Override protected void initArrays(int d) {
        super.initArrays(d); //[Not sure why the "d" is there, but was there in super.initArrays...]
        // this is called by init() which is called in StochasticProcess.setObservationTimes().
        integratedPath = new double[d + 1];
        bdt = new double[d];
        oneMExpMADtOverAlpha = new double[d];
        // iTime == d not included
        for (int iTime = 0; iTime < d; iTime++) {
            double dt = t[iTime + 1] - t[iTime];
            bdt[iTime] = beta * dt;
            double oneMExpMADt = -Math.expm1(-alpha * dt); // (1.0 - exp(-alpha * dt);
            oneMExpMADtOverAlpha[iTime] = oneMExpMADt / alpha;
        }

        double c11 = sigma * sigma / 2.0 / alpha;
        double c12 = c11 / alpha;
        c22 = c12 / alpha;
        normalsCorrGen = new MultinormalGen[d];
        for (int iTime = 0; iTime < d; iTime++) {
            double dt = t[iTime + 1] - t[iTime];
            double expMADt = alphadt[iTime];
            double oneMExpMADt = -Math.expm1(-alpha * dt); // (1.0 - exp(-alpha * dt);

            DenseDoubleMatrix2D covar = new DenseDoubleMatrix2D(2, 2);
            covar.set(0, 0, c11 * oneMExpMADt * (1.0 + expMADt)); // sigma^2/2/alpha * (1 - exp(-2 * alpha * dt))
            covar.set(0, 1, c12 * oneMExpMADt * oneMExpMADt); // sigma^2/2/alpha^2 * (1 - exp(-alpha*dt))^2
            covar.set(1, 0, covar.get(0, 1));
            covar.set(1, 1, c22 * (-3.0 + 2.0 * alpha * dt + 4.0 * expMADt - expMADt * expMADt));
            
            double[] mu = {0.0, 0.0}; // the expectation terms depend on the previous path values and are added elsewhere.
            normalsCorrGen[iTime] = new MultinormalPCAGen(gen, mu, covar);
        }
    }

    
    /**
     * Same as the method from the parent class, but also generates at each time
     * step a stochastic variable that corresponds to the integral of
     * the process over the time interval: $\int_{t_{i-1}}^{t_{0}} X(t) dt$.
     * Those integrated values are not returned by this method, but are stored
     * in memory and can be retrieved with {@link #getIntegratedPath() getIntegratedPath()}.
     * <p>
     * Two random uniforms (from the inner RandomStream) are required to generate each time step instead of just
     * one for the parent class.
     * @see #getIntegratedPath()
     */
    @Override public double[] generatePath() {
        path[0] = x0;
        integratedPath[0] = 0.0;
        double[] normalsCorr = new double[2];
        for (int iTime = 1; iTime <= d; iTime++) {
            normalsCorrGen[iTime - 1].nextPoint(normalsCorr);
            path[iTime] = badt[iTime - 1] + path[iTime - 1] * alphadt[iTime - 1] + normalsCorr[0];
            
            integratedPath[iTime] = integratedPath[iTime - 1] + 
                bdt[iTime - 1] + oneMExpMADtOverAlpha[iTime - 1] * (path[iTime - 1] - beta) + normalsCorr[1];
        }
        observationIndex = d;
        return path;
    }


    /**
     * This method must be called after {@link #generatePath() generatePath()} or 
     * {@link #nextObservation() nextObservation()}.
     * @return $\int_{\t_{0}}^{\t_{i}} X(\tau) d\tau$ for all (d + 1) time steps.
     */
    public double[] getIntegratedPath() {
        return integratedPath;
    }


    /**
     * {@code generatePath} must be called before this method is called since it depends on the underlying {@code path}.
     *
     * @param couponTimes The times at which the coupons are disbursed, which
     *  are later (by noticePeriod) than the time at which the decisions are taken (the times used
     *  to generate the values in this class). {@code couponTimes[0]} is meaningless (no coupon).
     * The order should correspond to the observation times of the process.  If {@code couponTimes}
     * has more elements than the number of time steps of the process, the expectation of those
     * additional terms is computed relative to the last time step of the process.
     * 
     * @return $E[exp(-\int_{\tau_m}^{_0} R(t) dt) | R(\tau_m),\int_{\tau_m}^{_0} R(t) dt) ]$ 
     * where $\tau_m$ is one of the {@code observationTimes}
     * of this {@code StochasticProcess} and $t_m$ is the corresponding {@code couponTimes}.
     * $t_m$ should be larger than $\tau_m$ by one notice period.
     * The first element of {@code couponTimes} and the returned array are meaningless; 
     * the returned value at {@code [0]} is {@code 0.0}.
     */
    public double[] getExpectedFutureDiscount(double[] couponTimes) {
        double[] analyticNotificationDiscounts = new double[couponTimes.length];
        for (int iTime = 1; iTime < couponTimes.length; iTime++) {
            int iTimePath = (iTime < d) ? iTime : d;
            double dt = couponTimes[iTime] - t[iTimePath]; // noticePeriod
            double expMADt = Math.exp(-alpha * dt);
            double mu_2 = beta * dt + (1.0 - expMADt) / alpha * (path[iTimePath] - beta);
            double sigma22 = c22 * (-3.0 + 2.0 * alpha * dt + 4.0 * expMADt - expMADt * expMADt);
            analyticNotificationDiscounts[iTime] = Math.exp(-mu_2  + sigma22 / 2.0 - integratedPath[iTimePath]);
        }
        return analyticNotificationDiscounts;
    }


    /**
     * This method is independent of the observation times of the object's instance, but it
     * is dependent on the other parameters of the process. 
     * It returns an analytic discount, which is the averaging discount at each time step..
     * @param times The times for which we want the total discounting (back to start time).
     * @return $E[exp(-\int_{t_0}^{t_m} R(t) dt) | R(t_0)]$ where $t_m$ is one of the observation times.
     *  The first element of the returned array is 1.0.
     */
    public double[] getTotalAnalyticDiscount(double[] times) {
        double[] analyticDiscount = new double[times.length];
        analyticDiscount[0] = 1.0;
        for (int iTime = 1; iTime < times.length; iTime++) {
            double dt = times[iTime] - times[0]; 
            double expMADt = Math.exp(-alpha * dt);
            double mu_2 = beta * dt + (1.0 - expMADt) / alpha * (x0 - beta);
            double sigma22 = c22 * (-3.0 + 2.0 * alpha * dt + 4.0 * expMADt - expMADt * expMADt);
            analyticDiscount[iTime] = Math.exp(-mu_2  + sigma22 / 2.0);
        }
        return analyticDiscount;
    }
}