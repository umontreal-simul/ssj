/*
 * Class:        MultivariateGeometricBrownianMotion
 * Description:
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Jean-Sébastien Parent & Clément Teule
 * @since        2008

 * SSJ is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, either version 3 of the License, or
 * any later version.

 * SSJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * A copy of the GNU General Public License is available at
   <a href="http://www.gnu.org/licenses">GPL licence site</a>.
 */
package umontreal.ssj.stochprocess;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.randvar.*;

/**
 * This class is a multivariate version of  @ref GeometricBrownianMotion. It
 * represents a multivariate GBM process @f$\{\mathbf{S}(t) =
 * (S_1(t),…,S_c(t)),  t\ge0\}@f$, which can be written as
 * @anchor REF_stochprocess_MultivariateGeometricBrownianMotion_eq_GBM2
 * @f[
 *   S_i(t) = S_i(0) \exp\left[ X_i(t) \right], \tag{GBM2}
 * @f]
 * where @f$\mathbf{X}(t) = (X_1(t),…,X_c(t))@f$ is a multivariate Brownian
 * motion. The GBM process is simulated by simulating the BM process
 * @f$\mathbf{X}@f$ (which is stored internally) and taking the exponential.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class MultivariateGeometricBrownianMotion
                           extends MultivariateStochasticProcess {

    protected NormalGen      gen;
    protected MultivariateBrownianMotion mbm;   // The underlying BM process X.
    protected double[]       mu,
                             sigma;
    protected double[]       mudt;

    /**
     * Constructs a new `MultivariateGeometricBrownianMotion` with
     * parameters @f$\mu= \mathtt{mu}@f$, @f$\sigma= \mathtt{sigma}@f$,
     * and @f$S(t_0) = \mathtt{x0}@f$, using `mbm` as the underlying
     * @ref MultivariateBrownianMotion. The parameters of `mbm` are
     * automatically reset to @f$\mu-\sigma^2/2@f$ and @f$\sigma@f$,
     * regardless of the original parameters of `mbm`. The correlation
     * structure is determined by the underlying
     * @ref MultivariateBrownianMotion. The observation times are the
     * same as those of `mbm`. The generation method depends on that of
     * `mbm` (sequential, bridge sampling, PCA, etc.).
     */
    public MultivariateGeometricBrownianMotion (int c, double[] x0,
                                                double[] mu, double[] sigma,
                                                MultivariateBrownianMotion mbm) {
        this.mbm = mbm;
        setParams (c, x0, mu, sigma);
    }

    /**
     * Sets the observation times of the
     * `MultivariateGeometricBrownianMotion`, but also those of the inner
     * @ref MultivariateBrownianMotion.
     */
    public void setObservationTimes (double[] t, int d) {
        this.d = d;
        super.setObservationTimes (t, d);
        mbm.setObservationTimes (t, d);
    }
public double[] nextObservationVector() {
        // Note : this implementation is general, to deal with
        // the possibility of generating bm with bridge sampling, for example.  ???

         // Note Clément :
         // Problème car x0 est un double[], pas un double
//          double s = x0 * Math.exp (mbm.nextObservationVector());
//          path[observationCounter] = s;
//
//         observationCounter = mbm.getCurrentObservationIndex();
//         // Could be different than simply 'observationCounter++' because of the
//         // possibility of Brownian bridge
//
//         return s;
       throw new UnsupportedOperationException ("nextObservationVector is not implemented ");
    }

/**
 * Generates and returns the vector of next observations
 */
public void nextObservationVector(double[] obs) {
       throw new UnsupportedOperationException ("nextObservationVector is not implemented ");
    }


    public double[] generatePath()  {
        for(int i = 0; i < c; i++)
            path[i] = x0[i];
        mbm.generatePath();
        for (int j = 1; j <= d; j++)
            for (int i = 0; i < c; i++)
                path[c * j + i] = x0[i] * Math.exp(mbm.getObservation(j, i));
        observationCounter = d;
        return path;
    }

/**
 * Same as in `StochasticProcess`, but also invokes `resetStartProcess` for
 * the underlying `BrownianMotion` object.
 */
public void resetStartProcess() {
        observationCounter = 0;
        mbm.resetStartProcess();
    }

    /**
     * Sets the parameters @f$S(t_0) = \mathtt{x0}@f$, @f$\mu=
     * \mathtt{mu}@f$ and @f$\sigma= \mathtt{sigma}@f$ of the process.
     * *Warning*: This method will recompute some quantities stored
     * internally, which may be slow if called repeatedly.
     */
    public void setParams (int c, double[] x0, double[] mu, double[] sigma) {
        if (x0.length < c)
            throw new IllegalArgumentException (
                     "x0 dimension :  "+ x0.length + " is smaller than the process dimension : " + c);
        if (mu.length < c)
            throw new IllegalArgumentException (
                     "mu dimension :  "+ mu.length + " is smaller than the process dimension : " + c);
        if (sigma.length < c)
            throw new IllegalArgumentException (
                     "sigma dimension :  "+ sigma.length + " is smaller than the process dimension : " + c);

        double[] zero = new double[c];
        this.c = c;
        this.x0 = x0;
        this.mu = new double[c];
        for(int i = 0; i < c; i++){
            this.mu[i]    = mu[i] - 0.5 * sigma[i] * sigma[i];
            zero[i] = 0.0;
        }
        mbm.setParams (zero, this.mu, sigma);  // reconstructs the MBM with the proper parameters
        if (observationTimesSet) init(); // Otherwise not needed.
    }

    /**
     * Resets the random stream for the underlying Brownian motion to
     * `stream`.
     */
    public void setStream (RandomStream stream) { (mbm.gen).setStream (stream); }

    /**
     * Returns the random stream for the underlying Brownian motion.
     */
    public RandomStream getStream () { return (mbm.gen).getStream (); }

    /**
     * Returns the normal random variate generator used.
     */
    public NormalGen getGen() { return gen; }

    /**
     * Returns a reference to the  @ref MultivariateBrownianMotion object
     * used to generate the process.
     */
    public MultivariateBrownianMotion getBrownianMotion() {
        return mbm;
    }


    protected void init() {
        super.init();
    }

}