/*
 * Class:        Palpha
 * Description:  computes the P_alpha figure of merit for a lattice
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        January 2009

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
package umontreal.ssj.discrepancy;
   import umontreal.ssj.hups.Rank1Lattice;
   import umontreal.ssj.util.Num;
   import umontreal.ssj.util.PrintfFormat;

/**
 * Extends the class  @ref Discrepancy and implements the methods required to
 * compute the @f$P_{\alpha}@f$ figure of merit for a lattice point set
 * @f$\Psi_s@f$ which is the intersection of a lattice @f$L@f$ and the unit
 * hypercube @f$[0, 1)^s@f$ in @f$s@f$ dimensions. @f$\Psi_s@f$ contains
 * @f$n@f$ points. For an arbitrary integer @f$\alpha> 1@f$, it is defined
 * as @cite mSLO94a, @cite vHIC98c, @cite vHIC01a&thinsp;
 * @f[
 *   P_{\alpha}(s) = \sum_{\mathbf{0\neq h}\in L_s^*} \|\mathbf{h}\|^{-\alpha},
 * @f]
 * where @f$L_s^*@f$ is the lattice dual to @f$L_s@f$, and the norm is
 * defined as @f$\|\mathbf{h}\| = \prod^s_{j=1} \max\{1, |h_j|\}@f$. When
 * @f$\alpha@f$ is even, @f$P_{\alpha}@f$ can be evaluated explicitly as
 * @anchor REF_discrepancy_Palpha_palpha_1
 * @f[
 *   \qquad P_{\alpha}(s) = -1  +  \frac{1}{n}\sum_{i=1}^n \prod_{j=1}^s \left[1 - \frac{(-1)^{\alpha/2}(2\pi)^{\alpha}}{\alpha!} B_{\alpha}(u_{ij})\right], \tag{palpha.1}
 * @f]
 * where @f$u_{ij}@f$ is the @f$j@f$-th coordinate of point @f$i@f$, and
 * @f$B_{\alpha}(x)@f$ is the Bernoulli polynomial of degree @f$\alpha@f$
 * (see  umontreal.ssj.util.Num.bernoulliPoly in class <tt>util/Num</tt>).
 *
 * One may generalize the @f$P_{\alpha}@f$ by introducing a weight for each
 * dimension to give the *weighted* @f$P_{\alpha}@f$ defined by
 * @cite vHIC98c&thinsp;
 * @f[
 *   P_{\alpha}(s) = \sum_{\mathbf{0\neq h}\in L_s^*}\beta_I^2 \|\mathbf{h}\|^{-\alpha},
 * @f]
 * where the weights are such that @f$\beta_I
 * =\beta_0\prod_{j=1}^s\beta_j^{\alpha}@f$, and for even @f$\alpha@f$
 * @anchor REF_discrepancy_Palpha_palpha_2
 * @f[
 *   \qquad P_{\alpha}(s) = \beta_0 \left\{-1  +  \frac{1}{n}\sum_{i=1}^n \prod_{j=1}^s \left[1 - \frac{(-1)^{\alpha/2}(2\pi\beta_j)^{\alpha}}{\alpha!} B_{\alpha}(u_{ij})\right] \right\}. \tag{palpha.2}
 * @f]
 * One recovers the original criterion (
 * {@link REF_discrepancy_Palpha_palpha_1 palpha.1} ) for
 * @f$P_{\alpha}@f$ by choosing all @f$\beta_j = 1@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class Palpha extends Discrepancy {
   // J'utilise les gamma de la classe mère pour beta[]
   private double[] C;         // beta * factors
   private int alpha;

   private void initFactor (int s, double[] Beta, int alpha) {
      C = new double[s];
      for (int j = 1; j <= s; j++)
         C[j-1] = Math.pow(2.0*Math.PI*Beta[j], (double) alpha) /
                  Num.factorial(alpha);
   }


   private void initA (int s, double[] bet, int alp) {
      if (null == bet) {
         setONES(s);
         gamma = ONES;
      } else
         gamma = bet;
      alpha = alp;
   }

   /**
    * Constructor with @f$n@f$ points in @f$s@f$ dimensions and with
    * `alpha` @f$= \alpha@f$. `points[i][j]` is the @f$j@f$-th coordinate
    * of point @f$i@f$. Both `i` and `j` start at 0. The weights
    * `beta[j]`, @f$j=0, 1,…, s@f$ are as in eq. (
    * {@link REF_discrepancy_Palpha_palpha_2 palpha.2} ).
    * The points and dimensions in (
    * {@link REF_discrepancy_Palpha_palpha_2 palpha.2} ) are
    * @f$u_{ij} = @f$ `points[i-1][j-1]`, but @f$\beta_j = @f$ `beta[j]`.
    * Restriction: `alpha` @f$ \in\{2, 4, 6, 8\}@f$.
    */
   public Palpha (double[][] points, int n, int s, double[] beta, int alpha) {
      super (points, n, s, beta);
      initA (s, beta, alpha);
   }

   /**
    * Constructor with all @f$\beta_j = 1@f$ (see eq. (
    * {@link REF_discrepancy_Palpha_palpha_1 palpha.1} )).
    */
   public Palpha (double[][] points, int n, int s, int alpha) {
      super(points, n, s);
      initA (s, null, alpha);
  }

   /**
    * Constructor with @f$n@f$ points in @f$s@f$ dimensions and with
    * `alpha` @f$= \alpha@f$. The @f$n@f$ points will be chosen later.
    * The weights `beta[j]`, @f$j=0, 1,…, s@f$ are as in eq. (
    * {@link REF_discrepancy_Palpha_palpha_2 palpha.2} ),
    * with @f$\beta_j = @f$ `beta[j]`. Restriction: `alpha` @f$ \in\{2,
    * 4, 6, 8\}@f$.
    */
   public Palpha (int n, int s, double[] beta, int alpha) {
      super (n, s, beta);
      initA (s, beta, alpha);
   }

   /**
    * Constructor with parameter `alpha` @f$= \alpha@f$. The points and
    * parameters *must* be defined before calling methods of this class.
    * Restriction: `alpha` @f$ \in\{2, 4, 6, 8\}@f$.
    */
   public Palpha (int alpha) {
      this.alpha = alpha;
   }

   /**
    * Constructor with the lattice `set` with weights
    * <tt>beta[j]</tt>@f$=\beta_j@f$ and parameter `alpha` @f$=
    * \alpha@f$. All the points are copied in an internal array.
    * Restriction: `alpha` @f$ \in\{2, 4, 6, 8\}@f$.
    */
   public Palpha (Rank1Lattice set, double[] beta, int alpha) {
      super(set);
      initA (set.getDimension(), beta, alpha);
   }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShift1Lattice_shift1lat
    * shift1lat} ) for the @f$s@f$-dimensional points of lattice `points`,
    * containing @f$n@f$ points. All weights @f$\beta_j = 1@f$.
    */
   public double compute (double[][] points, int n, int s) {
      return compute (points, n, s, null, alpha);
   }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShift1Lattice_shift1lat
    * shift1lat} ) for the @f$s@f$-dimensional points of lattice `points`,
    * containing @f$n@f$ points, with weights @f$\beta_j=@f$ `beta[j]`.
    */
   public double compute (double[][] points, int n, int s, double[] beta) {
      return compute (points, n, s, beta, alpha);
   }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShift1Lattice_shift1lat
    * shift1lat} ) for the @f$s@f$-dimensional points of lattice `points`,
    * containing @f$n@f$ points, with all weights @f$\beta_j=1@f$ and
    * @f$\alpha=@f$ `alpha`. Restriction: `alpha` @f$ \in\{2, 4, 6,
    * 8\}@f$.
    */
   public double compute (double[][] points, int n, int s, int alpha) {
      return compute (points, n, s, null, alpha);
   }

   /**
    * Computes the discrepancy (
    * {@link REF_discrepancy_DiscShift1Lattice_shift1lat
    * shift1lat} ) for the @f$s@f$-dimensional points of lattice `points`,
    * containing @f$n@f$ points, with weights @f$\beta_j=@f$ `beta[j]`
    * and with @f$\alpha=@f$ `alpha`. Restriction: `alpha` @f$ \in\{2,
    * 4, 6, 8\}@f$.
    */
   public double compute (double[][] points, int n, int s, double[] beta,
                          int alpha) {
      initA (s, beta, alpha);
      initFactor (s, gamma, alpha);

      switch (alpha) {
      case 2:
         return calcPalpha2 (points, n, s, gamma);
      case 4:
         return calcPalpha4 (points, n, s, gamma);
      case 6:
         return calcPalpha6 (points, n, s, gamma);
      case 8:
         return calcPalpha8 (points, n, s, gamma);
      default:
         throw new IllegalArgumentException ("alpha must be one of {2, 4, 6, 8}");
      }
   }

   /**
    * This method computes the value of @f$P_{\alpha}@f$ with @f$\alpha=
    * 2@f$ for the set of points in dimension @f$s@f$.
    */
   private double calcPalpha2 (double[][] points, int N, int s,
                                                double[] beta) {
      double sum = 0.0;
      for (int i = 0; i < N; ++i) {
         double prod = 1.0;
         for (int j = 0; j < s; ++j) {
            double u = Points[i][j];
            prod *= 1.0 + C[j] * (u*(u - 1.0) + UNSIX);
         }
         sum += prod;
      }
      sum /= N;
      return beta[0] * (sum - 1.0);
    }

   /**
    * This method computes the value of @f$P_{\alpha}@f$ with @f$\alpha=
    * 4@f$ for the set of points in dimension @f$s@f$.
    */
   private double calcPalpha4 (double[][] points, int N, int s,
                                      double[] beta) {
      double sum = 0.0;
      for (int i = 0; i < N; ++i) {
         double prod = 1.0;
         for (int j = 0; j < s; ++j) {
            double u = Points[i][j];
            prod *= 1.0 - C[j] * (((u - 2.0) * u + 1.0)*u*u - UNTRENTE);
         }
         sum += prod;
      }
      sum /= N;
      return beta[0] * (sum - 1.0);
    }

   /**
    * This method computes the value of @f$P_{\alpha}@f$ with @f$\alpha=
    * 6@f$ for the set of points in dimension @f$s@f$.
    */
   private double calcPalpha6 (double[][] points, int N, int s,
                                        double[] beta) {
      double sum = 0.0;
      for (int i = 0; i < N; ++i) {
         double prod = 1.0;
         for (int j = 0; j < s; ++j) {
            double u = Points[i][j];
            prod *= 1.0 + C[j] *
                    ((((u - 3.0) * u + 2.5) * u*u - 0.5) * u*u + QUARAN);
         }
         sum += prod;
      }
      sum /= N;
      return beta[0] * (sum - 1.0);
    }

   /**
    * This method computes the value of @f$P_{\alpha}@f$ with @f$\alpha=
    * 8@f$ for the set of points in dimension @f$s@f$.
    */
   private double calcPalpha8 (double[][] points, int N, int s,
                                  double[] beta) {
      double sum = 0.0;
      for (int i = 0; i < N; ++i) {
         double prod = 1.0;
         for (int j = 0; j < s; ++j) {
            double u = Points[i][j];
            prod *= 1.0 - C[j] * (((((u - 4.0) * u +
                  QTIERS) * u*u - STIERS) * u*u + DTIERS) * u*u - UNTRENTE);
         }
         sum += prod;
      }
      sum /= N;
      return beta[0] * (sum - 1.0);
    }
   public String toString() {
      StringBuffer sb = new StringBuffer (getName() + ":" +
                                          PrintfFormat.NEWLINE);
      sb .append ("n = " + numPoints + ",   dim = " + dim +
                  PrintfFormat.NEWLINE);
      sb .append ("gamma = [");
      appendGamma (sb, gamma, dim + 1);
      sb.append (" ]" + PrintfFormat.NEWLINE);
      sb.append ("alpha = " + alpha + PrintfFormat.NEWLINE);
      return sb.toString();
   }

   /**
    * Sets the values of @f$\beta_j = \mathtt{beta[j]}, j = 0, …, s@f$,
    * where @f$s@f$ is the dimension of the points.
    */
   public void setBeta (double[] beta) {
     gamma = beta;
    }

}