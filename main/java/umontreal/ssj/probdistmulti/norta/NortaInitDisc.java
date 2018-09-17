/*
 * Class:        NortaInitDisc
 * Description:
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Université de Montréal
 * Organization: DIRO, Université de Montréal
 * @author       Nabil Channouf
 * @since
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package umontreal.ssj.probdistmulti.norta;
import umontreal.ssj.probdistmulti.*;
import umontreal.ssj.probdist.*;

/**
 * This abstract class defines the algorithms used for NORTA initialization
 * when the marginal distributions are discrete. Four algorithms are
 * supported for now, and they are defined as subclasses of the class
 * `NortaInitDisc`.
 *
 * For two random variables @f$X_1@f$ and @f$X_2@f$, and their two marginal
 * distributions @f$F_1@f$ and @f$F_2@f$, respectively, we specify the rank
 * correlation @f$r_X=\mbox{Corr}(F_1(X_1),F_2(X_2))@f$, the parameters of
 * the marginal distributions and a parameter for truncation @f$tr@f$. For
 * the correlation matching, we must have finite supports for the two
 * distributions. Then if the support of each marginal is infinite, we have
 * to upper-bound it at the quantile of order @f$tr@f$. For example, if the
 * marginals have their support points in @f$[0,+\infty)@f$, the software
 * will truncate to @f$[0,F_l^{-1}(tr)]@f$, for @f$l=1,2@f$. The parameter
 * @f$tr@f$ must to be given by the user, depending on the type of the two
 * distributions. If the marginals have finite supports, one can simply give
 * @f$tr=1@f$.
 *
 * Each algorithm  @ref NI1,  @ref NI2a,  @ref NI2b and  @ref NI3 can be used
 * to calculate the corresponding correlation
 * @f$\rho_Z=\mbox{Corr}(Z_1,Z_2)@f$, where @f$Z_1@f$ and @f$Z_2@f$ are
 * standard normal random variables. These subclasses implement the specific
 * methods for NORTA initialization presented in @cite sAVR06a&thinsp;.
 *
 * Each type of algorithm should be defined as a subclass of
 * @ref NortaInitDisc. Each subclass must implement the method  #computeCorr
 * which returns the solution @f$\rho_Z@f$. When executing this method, the
 * subclass may call the methods  #integ and  #deriv, depending on the type
 * of algorithm. For example, the subclass  @ref NI1 calls only the method
 * #integ, since the algorithm do not use the derivative
 * @cite sAVR06a&thinsp;. Each subclass must also call the method
 * #computeParams(,) which is executed immediately before the beginning of
 * the root-finder algorithm.
 *
 * When creating a class representing an algorithm, the  #toString method can
 * be called to display information about the inputs.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public abstract class NortaInitDisc {
   protected double rX; // Target rank correlation
   // Marginal disttributions of r.variables X_1 and X_2
   protected DiscreteDistributionInt dist1;
   protected DiscreteDistributionInt dist2;
   protected double tr; /* Parameter for the quantile upper limit;
   		    used for truncation in the unbounded case. */
   protected double mu1, mu2, sd1, sd2; /* Means and standard deviations
   					    of F_1(X_1) and F_2(X_2). */

   private int m1, m2;    // Number of support points of the 2 marginals.
   private double[] p1;   // Probability masses of the marginal 1.
   private double[] p2;   // Probability masses of the marginal 2.
   /* Quantiles of the cumulative probability masses by the standard
      normal C.D.F for the 2 marginals */
   private double[] z1;
   private double[] z2;


   private String tabToString(double[] tab , String message)
   {
      String desc = message + "\n [";
      for (int i = 0; i < tab.length; i++)
         if (i == tab.length - 1)
            desc += "]\n";
         else
            desc += tab[i] + ",";
      return desc;
   }

   /**
    * Constructor with the target rank correlation `rX`, the two discrete
    * marginals `dist1` and `dist2` and the parameter for the truncation
    * `tr`. This constructor can be called only by the constructors of the
    * subclasses.
    */
   public NortaInitDisc (double rX, 
                         DiscreteDistributionInt dist1, 
                         DiscreteDistributionInt dist2,
                         double tr) {
      this.rX = rX;
      this.dist1 = dist1;
      this.dist2 = dist2;
      this.tr = tr;
   }

   /**
    * Computes and returns the correlation @f$\rho_Z@f$. Every subclass
    * of  @ref NortaInitDisc must implement this method.
    */
   public abstract double computeCorr();

   /**
    * Computes the following inputs of each marginal distribution:
    * <ul><li>
    * The number of support points @f$m_1@f$ and @f$m_2@f$ for the two
    * distributions.
    * </li>
    * <li>
    * The means and standard deviations of @f$F_1(X_1)@f$ and
    * @f$F_2(X_2)@f$, respectively.
    * </li>
    * <li>
    * The vectors @f$p_1[i]@f$, @f$p_2[j]@f$,
    * @f$z_1[i]=\Phi^{-1}(f_1[i])@f$ and @f$z_2[j]=\Phi^{-1}(f_2[j])@f$,
    * where @f$f_1[i]@f$ and @f$f_2[j]@f$, for @f$i=0,…,m_1-1@f$;
    * @f$j=0,…,m_2-1@f$, are the cumulative probability functions, and
    * @f$\Phi@f$ is the standard normal distribution function.
    * </li>
    * </ul> Every subclass of  @ref NortaInitDisc must call this method.
    */
   public void computeParams ()
   {
      m1 = dist1.inverseFInt (tr) + 1;
      m2 = dist2.inverseFInt (tr) + 1;
      // Support points of X_1 and X_2
      int[] y1 = new int[m1];
      int[] y2 = new int[m2];
      p1 = new double[m1];
      p2 = new double[m2];
      // Cumulative probability masses of X_1 and X_2
      double[] f1 = new double[m1];
      double[] f2 = new double[m2];
      z1 = new double[m1];
      z2 = new double[m2];
      double u11 = 0.0, u22 = 0.0;

      // Compute mu1, sd1, p1 and z1 of X_1
      for (int i = 0; i < m1; i++) {
         y1[i] = i;
         p1[i] = dist1.prob (y1[i]);
         f1[i] = dist1.cdf (y1[i]);
         z1[i] = NormalDist.inverseF01 (f1[i]);
         if (z1[i] == Double.NEGATIVE_INFINITY)
            z1[i] = NormalDist.inverseF01 (2.2e-308);
         if (z1[i] == Double.POSITIVE_INFINITY)
            z1[i] = NormalDist.inverseF01 (1.0 - Math.ulp (1.0));
         mu1 += f1[i] * p1[i];
         u11 += f1[i] * f1[i] * p1[i];
      }
      sd1 = Math.sqrt (u11 - mu1 * mu1);

      // Compute mu2, sd2, p2 and z2 of X_2
      for (int i = 0; i < m2; i++) {
         y2[i] = i;
         p2[i] = dist2.prob (y2[i]);
         f2[i] = dist2.cdf (y2[i]);
         z2[i] = NormalDist.inverseF01 (f2[i]);
         if (z2[i] == Double.NEGATIVE_INFINITY)
            z2[i] = NormalDist.inverseF01 (2.2e-308);
         if (z2[i] == Double.POSITIVE_INFINITY)
            z2[i] = NormalDist.inverseF01 (1.0 - Math.ulp (1.0));
         mu2 += f2[i] * p2[i];
         u22 += f2[i] * f2[i] * p2[i];
      }
      sd2 = Math.sqrt (u22 - mu2 * mu2);
   }

   /**
    * Computes the function
    * @anchor REF_probdistmulti_norta_NortaInitDisc_gr_M
    * @f{align}{
    *    \tag{gr_M} g_r(r)=\sum_{i=0}^{m_1-2} p_{1,i+1} \sum_{j=0}^{m_2-2} p_{2,j+1} \bar{\Phi}_r(z_{1,i },z_{2,j}),
    * @f}
    * which involves the bivariate normal integral @f$ \bar{\Phi}_r(x,y)=
    * \int_x^{\infty}\int_y^{\infty}\phi_r(z_1,z_2) dz_1 dz_2@f$.
    * Method  #barF of class
    * @ref umontreal.ssj.probdistmulti.BiNormalDonnellyDist (from package
    * probdistmulti of SSJ @cite iLEC04j&thinsp;) is used to compute
    * @f$\bar{\Phi}_r(x,y)@f$, with @f$m_1@f$, @f$m_2@f$, and the vectors
    * @f$p_1[i], i=1,…,m_1-1@f$; @f$z_1[i], i=0,…,m_1-2@f$; @f$p_2[j],
    * j=1,…,m_2-1@f$; @f$z_2[j], j=0,…,m_2-2@f$. The correlation parameter
    * @f$r@f$ must be in @f$[-1,1]@f$. This method may be called by
    * subclasses of  @ref NortaInitDisc.
    *  @param r            initial value of the correlation
    */
   public double integ (double r)
   {
      double gr = 0.0; // The returned value.
      for (int i = 0; i < m1 - 1; i++) {
         double sum = 0.0;
         for (int j = 0; j < m2 - 1; j++) {
            sum += p2[j + 1]
                   * BiNormalDonnellyDist.barF (z1[i], z2[j], r);
         }
         gr += p1[i + 1] * sum;
      }
      return gr;
   }

   /**
    * Computes the derivative of @f$g_r@f$, given by
    * @anchor REF_probdistmulti_norta_NortaInitDisc_grp_M
    * @f{align}{
    *    \tag{grp_M} g’_r(r)=\sum_{i=0}^{m_1-2} p_{1,i+1} \sum_{j=0}^{m_2-2} p_{2,j+1} \phi_r(z_{1,i },z_{2,j}),
    * @f}
    * where @f$\phi_r@f$ is the bivariate standard binormal density. The
    * method uses @f$m_1@f$, @f$m_2@f$, and the vectors @f$p_1[i],
    * i=1,…,m_1-1@f$; @f$z_1[i], i=0,…,m_1-2@f$; @f$p_2[j],
    * j=1,…,m_2-1@f$; @f$z_2[j], j=0,…,m_2-2@f$. The correlation parameter
    * @f$r@f$ must be in @f$[-1,1]@f$. This method may be called by
    * subclasses of  @ref NortaInitDisc.
    *  @param r            initial value of the correlation
    */
   public double deriv (double r)
   {
      double c = Math.sqrt (1.0 - r * r);
      double c1 = 2 * c * c;
      double gp = 0.0; // The returned value

      for (int i = 0; i < m1 - 1; i++) {
         double z1sq = z1[i] * z1[i];
         double t1 = 2 * r * z1[i];
         double sum = 0.0;
         for (int j = 0; j < m2 - 1; j++) {
            sum += p2[j + 1]
                   * Math.exp ((t1 * z2[j] - z1sq - z2[j] * z2[j]) / c1);
         }
         gp += p1[i + 1] * sum;
      }
      gp = gp / (2.0 * Math.PI * c);
      return gp;
   }
   public String toString()
   {
      String desc = "";
      desc += "rX = " + rX + "\n";
      desc += "tr = " + tr + "\n";
      desc += "m1 = " + m1 + "\n";
      desc += "m2 = " + m2 + "\n";
      desc += "mu1 = " + mu1 + "\n";
      desc += "mu2 = " + mu2 + "\n";
      desc += "sd1 = " + sd1 + "\n";
      desc += "sd2 = " + sd2 + "\n";
      desc += tabToString(p1 , "Table p1 : ");
      desc += tabToString(z1 , "Table z1 : ");
      desc += tabToString(p2 , "Table p2 : ");
      desc += tabToString( z2, "Table z2 : ");
      return desc;
   }

}
