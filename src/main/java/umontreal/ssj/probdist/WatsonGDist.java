/*
 * Class:        WatsonGDist
 * Description:  Watson G distribution
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
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
package umontreal.ssj.probdist;
import umontreal.ssj.util.*;
import umontreal.ssj.functions.MathFunction;

/**
 * Extends the class  @ref ContinuousDistribution for the Watson @f$G@f$
 * distribution (see @cite tDAR83a, @cite tWAT76a&thinsp;). Given a sample of
 * @f$n@f$ independent uniforms @f$U_i@f$ over @f$[0,1]@f$, the @f$G@f$
 * statistic is defined by
 * @anchor REF_probdist_WatsonGDist_eq_WatsonG
 * @f{align}{
 *    G_n 
 *    & 
 *   =
 *    \sqrt{n} \max_{\Rule{0.0pt}{7.0pt}{0.0pt} 1\le j \le n} \left\{ j/n - U_{(j)} + \bar{U}_n - 1/2 \right\} \tag{WatsonG} 
 *    \\  & 
 *   =
 *    \sqrt{n}\left(D_n^+ + \bar{U}_n - 1/2\right), \nonumber
 * @f}
 * where the @f$U_{(j)}@f$ are the @f$U_i@f$ sorted in increasing order,
 * @f$\bar{U}_n@f$ is the average of the observations @f$U_i@f$, and
 * @f$D_n^+@f$ is the Kolmogorov-Smirnov+ statistic. The distribution
 * function (the cumulative probabilities) is defined as @f$F_n(x) = P[G_n
 * \le x]@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_edf
 */
public class WatsonGDist extends ContinuousDistribution {
   protected int n;

   private static class Function implements MathFunction {
      protected int n;
      protected double u;

      public Function (int n, double u) {
         this.n = n;
         this.u = u;
      }

      public double evaluate (double x) {
         return u - cdf(n,x);
      }
   }

   /**
    * Constructs a *Watson* distribution for a sample of size @f$n@f$.
    */
   public WatsonGDist (int n) {
      setN (n);
   }


   public double density (double x) {
      return density (n, x);
   }

   public double cdf (double x) {
      return cdf (n, x);
   }

   public double barF (double x) {
      return barF (n, x);
   }

   public double inverseF (double u) {
      return inverseF (n, u);
   }

/**
 * Computes the density function for a *Watson* @f$G@f$ distribution with
 * parameter @f$n@f$.
 */
public static double density (int n, double x) {
      final double MINARG = 0.15;
      final double MAXARG = 1.5;

      if (n < 2)
        throw new IllegalArgumentException ("n < 2");

      if (x <= MINARG || x >= XBIGM)
         return 0.0;

      final double Res;
      if (x > MAXARG)
         Res = 20 * Math.exp (19.0 - 20.0*x) 
            - 15.26 * Math.exp (13.34 - 15.26*x) / Math.sqrt ((double)n);
      else {
         final double EPS = 1.0 / 20.0;
         Res = (cdf(n, x + EPS) - cdf(n, x - EPS)) / (2.0 * EPS);
      }

      if (Res <= 0.0)
         return 0.0;
      return Res;
   }

   // Tables for a spline approximation
   // of the WatsonG distribution
   // Empirical correction in 1/sqrt (n)
   private static double YWA[] = new double[143];
   private static double MWA[] = new double[143];
   private static double CoWA[] = new double[143];

   static {
   /*
    * Initialization for watsonG
    */
      int j;
      YWA[0] = 1.8121832847E-39;
      YWA[1] = 2.0503176304E-32;
      YWA[2] = 4.6139577764E-27;
      YWA[3] = 6.5869745929E-23;
      YWA[4] = 1.2765816107E-19;
      YWA[5] = 5.6251923105E-17;
      YWA[6] = 8.0747150511E-15;
      YWA[7] = 4.8819994144E-13;
      YWA[8] = 1.4996052497E-11;
      YWA[9] = 2.6903519441E-10;
      YWA[10] = 3.1322929018E-9;
      YWA[11] = 2.5659643046E-8;
      YWA[12] = 1.5749759318E-7;
      YWA[13] = 7.6105096466E-7;
      YWA[14] = 3.0113293541E-6;
      YWA[15] = 1.0070166837E-5;
      YWA[16] = 2.9199826692E-5;
      YWA[17] = 7.4970409372E-5;
      YWA[18] = 1.7340586581E-4;
      YWA[19] = 3.6654236297E-4;
      YWA[20] = 7.165864865E-4;
      YWA[21] = 1.3087767385E-3;
      YWA[22] = 2.2522044209E-3;
      YWA[23] = 3.6781862572E-3;
      YWA[24] = 5.7361958631E-3;
      YWA[25] = 8.5877444706E-3;
      YWA[26] = 1.23988738E-2;
      YWA[27] = 1.73320516E-2;
      YWA[28] = 2.35382479E-2;
      YWA[29] = 3.11498548E-2;
      YWA[30] = 4.02749297E-2;
      YWA[31] = 5.09930445E-2;
      YWA[32] = 6.33528333E-2;
      YWA[33] = 7.73711747E-2;
      YWA[34] = 9.30338324E-2;
      YWA[35] = 1.10297306E-1;
      YWA[36] = 1.290916098E-1;
      YWA[37] = 1.493236984E-1;
      YWA[38] = 1.708812741E-1;
      YWA[39] = 1.936367476E-1;
      YWA[40] = 2.174511609E-1;
      YWA[41] = 2.42177928E-1;
      YWA[42] = 2.676662852E-1;
      YWA[43] = 2.937643828E-1;
      YWA[44] = 3.203219784E-1;
      YWA[45] = 3.471927188E-1;
      YWA[46] = 3.742360163E-1;
      YWA[47] = 4.013185392E-1;
      YWA[48] = 4.283153467E-1;
      YWA[49] = 4.551107027E-1;
      YWA[50] = 4.815986082E-1;
      YWA[51] = 5.076830902E-1;
      YWA[52] = 5.332782852E-1;
      YWA[53] = 5.583083531E-1;
      YWA[54] = 5.827072528E-1;
      YWA[55] = 6.064184099E-1;
      YWA[56] = 6.293943006E-1;
      YWA[57] = 6.515959739E-1;
      YWA[58] = 6.729925313E-1;
      YWA[59] = 6.935605784E-1;
      YWA[60] = 7.132836621E-1;
      YWA[61] = 7.321517033E-1;
      YWA[62] = 7.501604333E-1;
      YWA[63] = 7.673108406E-1;
      YWA[64] = 7.836086337E-1;
      YWA[65] = 7.99063723E-1;
      YWA[66] = 8.136897251E-1;
      YWA[67] = 8.275034914E-1;
      YWA[68] = 8.405246632E-1;
      YWA[69] = 8.527752531E-1;
      YWA[70] = 8.642792535E-1;
      YWA[71] = 8.750622738E-1;
      YWA[72] = 8.851512032E-1;
      YWA[73] = 8.945739017E-1;
      YWA[74] = 9.033589176E-1;
      YWA[75] = 9.115352296E-1;
      YWA[76] = 9.19132015E-1;
      YWA[77] = 9.261784413E-1;
      YWA[78] = 9.327034806E-1;
      YWA[79] = 9.387357465E-1;
      YWA[80] = 9.44303351E-1;
      YWA[81] = 9.494337813E-1;
      YWA[82] = 9.541537951E-1;
      YWA[83] = 9.584893325E-1;
      YWA[84] = 9.624654445E-1;
      YWA[85] = 9.661062352E-1;
      YWA[86] = 9.694348183E-1;
      YWA[87] = 9.724732859E-1;
      YWA[88] = 9.752426872E-1;
      YWA[89] = 9.777630186E-1;
      YWA[90] = 9.800532221E-1;
      YWA[91] = 9.821311912E-1;
      YWA[92] = 9.840137844E-1;
      YWA[93] = 9.85716844E-1;
      YWA[94] = 9.872552203E-1;
      YWA[95] = 9.886428002E-1;
      YWA[96] = 9.898925389E-1;
      YWA[97] = 9.910164946E-1;
      YWA[98] = 9.920258656E-1;
      YWA[99] = 9.929310287E-1;
      YWA[100] = 9.937415788E-1;
      YWA[101] = 9.944663692E-1;
      YWA[102] = 9.95113552E-1;
      YWA[103] = 9.956906185E-1;
      YWA[104] = 9.962044387E-1;
      YWA[105] = 9.966613009E-1;
      YWA[106] = 9.970669496E-1;
      YWA[107] = 9.974266225E-1;
      YWA[108] = 9.977450862E-1;
      YWA[109] = 9.980266707E-1;
      YWA[110] = 9.982753021E-1;
      YWA[111] = 9.984945338E-1;
      YWA[112] = 9.98687576E-1;
      YWA[113] = 9.98857324E-1;
      YWA[114] = 9.990063842E-1;
      YWA[115] = 9.991370993E-1;
      YWA[116] = 9.992515708E-1;
      YWA[117] = 9.99351681E-1;
      YWA[118] = 9.994391129E-1;
      YWA[119] = 9.995153688E-1;
      YWA[120] = 9.995817875E-1;
      YWA[121] = 9.996395602E-1;
      YWA[122] = 9.996897446E-1;
      YWA[123] = 9.997332791E-1;
      YWA[124] = 9.997709943E-1;
      YWA[125] = 9.998036243E-1;
      YWA[126] = 9.998318172E-1;
      YWA[127] = 9.998561438E-1;
      YWA[128] = 9.998771066E-1;
      YWA[129] = 9.998951466E-1;
      YWA[130] = 9.999106508E-1;
      YWA[131] = 9.99923958E-1;
      YWA[132] = 9.999353645E-1;
      YWA[133] = 9.999451288E-1;
      YWA[134] = 9.999534765E-1;
      YWA[135] = 9.999606035E-1;
      YWA[136] = 9.999666805E-1;
      YWA[137] = 9.999718553E-1;
      YWA[138] = 9.999762562E-1;
      YWA[139] = 9.999799939E-1;
      YWA[140] = 9.999831643E-1;
      YWA[141] = 9.999858E-1;
      YWA[142] = 9.999883E-1;

      MWA[0] = 0.0;
      MWA[1] = 6.909E-15;
      MWA[2] = 2.763E-14;
      MWA[3] = 1.036E-13;
      MWA[4] = 3.792E-13;
      MWA[5] = 4.773E-12;
      MWA[6] = 4.59E-10;
      MWA[7] = 2.649E-8;
      MWA[8] = 7.353E-7;
      MWA[9] = 1.14E-5;
      MWA[10] = 1.102E-4;
      MWA[11] = 7.276E-4;
      MWA[12] = 3.538E-3;
      MWA[13] = 0.01342;
      MWA[14] = 0.04157;
      MWA[15] = 0.1088;
      MWA[16] = 0.2474;
      MWA[17] = 0.4999;
      MWA[18] = 0.913;
      MWA[19] = 1.53;
      MWA[20] = 2.381;
      MWA[21] = 3.475;
      MWA[22] = 4.795;
      MWA[23] = 6.3;
      MWA[24] = 7.928;
      MWA[25] = 9.602;
      MWA[26] = 11.24;
      MWA[27] = 12.76;
      MWA[28] = 14.1;
      MWA[29] = 15.18;
      MWA[30] = 15.98;
      MWA[31] = 16.47;
      MWA[32] = 16.64;
      MWA[33] = 16.49;
      MWA[34] = 16.05;
      MWA[35] = 15.35;
      MWA[36] = 14.41;
      MWA[37] = 13.28;
      MWA[38] = 12.0;
      MWA[39] = 10.6;
      MWA[40] = 9.13;
      MWA[41] = 7.618;
      MWA[42] = 6.095;
      MWA[43] = 4.588;
      MWA[44] = 3.122;
      MWA[45] = 1.713;
      MWA[46] = 0.3782;
      MWA[47] = -0.8726;
      MWA[48] = -2.031;
      MWA[49] = -3.091;
      MWA[50] = -4.051;
      MWA[51] = -4.91;
      MWA[52] = -5.668;
      MWA[53] = -6.327;
      MWA[54] = -6.893;
      MWA[55] = -7.367;
      MWA[56] = -7.756;
      MWA[57] = -8.064;
      MWA[58] = -8.297;
      MWA[59] = -8.46;
      MWA[60] = -8.56;
      MWA[61] = -8.602;
      MWA[62] = -8.591;
      MWA[63] = -8.533;
      MWA[64] = -8.433;
      MWA[65] = -8.296;
      MWA[66] = -8.127;
      MWA[67] = -7.93;
      MWA[68] = -7.709;
      MWA[69] = -7.469;
      MWA[70] = -7.212;
      MWA[71] = -6.943;
      MWA[72] = -6.663;
      MWA[73] = -6.378;
      MWA[74] = -6.087;
      MWA[75] = -5.795;
      MWA[76] = -5.503;
      MWA[77] = -5.213;
      MWA[78] = -4.927;
      MWA[79] = -4.646;
      MWA[80] = -4.371;
      MWA[81] = -4.103;
      MWA[82] = -3.843;
      MWA[83] = -3.593;
      MWA[84] = -3.352;
      MWA[85] = -3.12;
      MWA[86] = -2.899;
      MWA[87] = -2.689;
      MWA[88] = -2.489;
      MWA[89] = -2.3;
      MWA[90] = -2.121;
      MWA[91] = -1.952;
      MWA[92] = -1.794;
      MWA[93] = -1.645;
      MWA[94] = -1.506;
      MWA[95] = -1.377;
      MWA[96] = -1.256;
      MWA[97] = -1.144;
      MWA[98] = -1.041;
      MWA[99] = -0.9449;
      MWA[100] = -0.8564;
      MWA[101] = -0.775;
      MWA[102] = -0.7001;
      MWA[103] = -0.6315;
      MWA[104] = -0.5687;
      MWA[105] = -0.5113;
      MWA[106] = -0.459;
      MWA[107] = -0.4114;
      MWA[108] = -0.3681;
      MWA[109] = -0.3289;
      MWA[110] = -0.2934;
      MWA[111] = -0.2614;
      MWA[112] = -0.2325;
      MWA[113] = -0.2064;
      MWA[114] = -0.183;
      MWA[115] = -0.1621;
      MWA[116] = -0.1433;
      MWA[117] = -0.1265;
      MWA[118] = -0.1115;
      MWA[119] = -9.813E-2;
      MWA[120] = -8.624E-2;
      MWA[121] = -7.569E-2;
      MWA[122] = -6.632E-2;
      MWA[123] = -5.803E-2;
      MWA[124] = -5.071E-2;
      MWA[125] = -4.424E-2;
      MWA[126] = -3.855E-2;
      MWA[127] = -3.353E-2;
      MWA[128] = -2.914E-2;
      MWA[129] = -2.528E-2;
      MWA[130] = -0.0219;
      MWA[131] = -1.894E-2;
      MWA[132] = -1.637E-2;
      MWA[133] = -1.412E-2;
      MWA[134] = -1.217E-2;
      MWA[135] = -1.046E-2;
      MWA[136] = -8.988E-3;
      MWA[137] = -7.72E-3;
      MWA[138] = -6.567E-3;
      MWA[139] = -5.802E-3;
      MWA[140] = -0.0053;
      MWA[141] = -4.7E-4;
      MWA[142] = -4.3E-4;

      for (j = 5; j <= 11; j++) {
         CoWA[j] = 0.0;
      }
      CoWA[12] = 1.25E-5;
      CoWA[13] = 3.87E-5;
      CoWA[14] = 1.004E-4;
      CoWA[15] = 2.703E-4;
      CoWA[16] = 6.507E-4;
      CoWA[17] = 1.3985E-3;
      CoWA[18] = 2.8353E-3;
      CoWA[19] = 5.1911E-3;
      CoWA[20] = 8.9486E-3;
      CoWA[21] = 1.41773E-2;
      CoWA[22] = 2.16551E-2;
      CoWA[23] = 3.1489E-2;
      CoWA[24] = 4.34123E-2;
      CoWA[25] = 5.78719E-2;
      CoWA[26] = 7.46921E-2;
      CoWA[27] = 9.45265E-2;
      CoWA[28] = 1.165183E-1;
      CoWA[29] = 1.406353E-1;
      CoWA[30] = 1.662849E-1;
      CoWA[31] = 1.929895E-1;
      CoWA[32] = 2.189347E-1;
      CoWA[33] = 2.457772E-1;
      CoWA[34] = 2.704794E-1;
      CoWA[35] = 2.947906E-1;
      CoWA[36] = 3.169854E-1;
      CoWA[37] = 3.377435E-1;
      CoWA[38] = 3.573555E-1;
      CoWA[39] = 3.751205E-1;
      CoWA[40] = 3.906829E-1;
      CoWA[41] = 4.039806E-1;
      CoWA[42] = 4.142483E-1;
      CoWA[43] = 4.22779E-1;
      CoWA[44] = 4.288013E-1;
      CoWA[45] = 4.330353E-1;
      CoWA[46] = 4.34452E-1;
      CoWA[47] = 4.338138E-1;
      CoWA[48] = 4.31504E-1;
      CoWA[49] = 4.272541E-1;
      CoWA[50] = 4.220568E-1;
      CoWA[51] = 4.158229E-1;
      CoWA[52] = 4.083281E-1;
      CoWA[53] = 3.981182E-1;
      CoWA[54] = 3.871678E-1;
      CoWA[55] = 3.755527E-1;
      CoWA[56] = 3.628823E-1;
      CoWA[57] = 3.520135E-1;
      CoWA[58] = 3.400924E-1;
      CoWA[59] = 3.280532E-1;
      CoWA[60] = 3.139477E-1;
      CoWA[61] = 2.997087E-1;
      CoWA[62] = 2.849179E-1;
      CoWA[63] = 2.710475E-1;
      CoWA[64] = 2.576478E-1;
      CoWA[65] = 2.449155E-1;
      CoWA[66] = 2.317447E-1;
      CoWA[67] = 2.193161E-1;
      CoWA[68] = 2.072622E-1;
      CoWA[69] = 1.956955E-1;
      CoWA[70] = 1.846514E-1;
      CoWA[71] = 1.734096E-1;
      CoWA[72] = 1.622678E-1;
      CoWA[73] = 1.520447E-1;
      CoWA[74] = 1.416351E-1;
      CoWA[75] = 1.32136E-1;
      CoWA[76] = 1.231861E-1;
      CoWA[77] = 1.150411E-1;
      CoWA[78] = 1.071536E-1;
      CoWA[79] = 9.9465E-2;
      CoWA[80] = 9.22347E-2;
      CoWA[81] = 8.54394E-2;
      CoWA[82] = 7.87697E-2;
      CoWA[83] = 7.23848E-2;
      CoWA[84] = 6.6587E-2;
      CoWA[85] = 6.15849E-2;
      CoWA[86] = 5.6573E-2;
      CoWA[87] = 5.17893E-2;
      CoWA[88] = 4.70011E-2;
      CoWA[89] = 4.2886E-2;
      CoWA[90] = 3.91224E-2;
      CoWA[91] = 3.53163E-2;
      CoWA[92] = 3.20884E-2;
      CoWA[93] = 2.92264E-2;
      CoWA[94] = 2.66058E-2;
      CoWA[95] = 2.37352E-2;
      CoWA[96] = 2.14669E-2;
      CoWA[97] = 1.94848E-2;
      CoWA[98] = 1.75591E-2;
      CoWA[99] = 1.58232E-2;
      CoWA[100] = 1.40302E-2;
      CoWA[101] = 1.24349E-2;
      CoWA[102] = 1.11856E-2;
      CoWA[103] = 9.9765E-3;
      CoWA[104] = 8.9492E-3;
      CoWA[105] = 8.0063E-3;
      CoWA[106] = 7.1509E-3;
      CoWA[107] = 6.3196E-3;
      CoWA[108] = 5.6856E-3;
      CoWA[109] = 5.0686E-3;
      CoWA[110] = 4.5085E-3;
      CoWA[111] = 3.9895E-3;
      CoWA[112] = 3.4804E-3;
      CoWA[113] = 3.0447E-3;
      CoWA[114] = 2.7012E-3;
      CoWA[115] = 2.2984E-3;
      CoWA[116] = 2.0283E-3;
      CoWA[117] = 1.7399E-3;
      CoWA[118] = 1.5032E-3;
      CoWA[119] = 1.3267E-3;
      CoWA[120] = 1.1531E-3;
      CoWA[121] = 9.92E-4;
      CoWA[122] = 9.211E-4;
      CoWA[123] = 8.296E-4;
      CoWA[124] = 6.991E-4;
      CoWA[125] = 5.84E-4;
      CoWA[126] = 5.12E-4;
      CoWA[127] = 4.314E-4;
      CoWA[128] = 3.593E-4;
      CoWA[129] = 3.014E-4;
      CoWA[130] = 2.401E-4;
      CoWA[131] = 2.004E-4;
      CoWA[132] = 1.614E-4;
      CoWA[133] = 1.257E-4;
      CoWA[134] = 1.112E-4;
      CoWA[135] = 9.22E-5;
      CoWA[136] = 8.77E-5;
      CoWA[137] = 6.22E-5;
      CoWA[138] = 4.93E-5;
      CoWA[139] = 3.92E-5;
      CoWA[140] = 3.15E-5;
      CoWA[141] = 1.03E-5;
      CoWA[142] = 9.6E-6;
   }

/**
 * Computes the Watson @f$G@f$ distribution function @f$F_n(x)@f$, with
 * parameter @f$n@f$. A cubic spline interpolation is used for the asymptotic
 * distribution when @f$n\to\infty@f$, and an empirical correction of order
 * @f$1/\sqrt{n}@f$, obtained empirically from @f$10^7@f$ simulation runs
 * with @f$n = 256@f$ is then added. The absolute error is estimated to be
 * less than 0.01, 0.005, 0.002, 0.0008, 0.0005, 0.0005, 0.0005 for @f$n =
 * 16@f$, 32, 64, 128, 256, 512, 1024, respectively.
 */
public static double cdf (int n, double x) {
     /*
      * Approximation of the cumulative distribution function of the
      * watsonG statistics by the cubic spline function.
      *   Y[.]  - tabular value of the statistic;
      *   M[.]  - tabular value of the first derivative;
      */
      if (n <= 1)
        throw new IllegalArgumentException ("n < 2");

      final double MINARG = 0.15;
      if (x <= MINARG)
         return 0.0;
      if (x >= 10.0)
         return 1.0;

      double R, Res;
      final double MAXARG = 1.5;
      if (x > MAXARG) {
         R = Math.exp (19.0 - 20.0*x);
         Res = 1.0 - R;
         // Empirical Correction in 1/sqrt (n)
         R = Math.exp (13.34 - 15.26*x)/Math.sqrt ((double)n);
         Res += R;
         // The correction in 1/sqrt (n) is not always precise
         if (Res >= 1.0)
            return 1.0;
         else
            return Res;
      }

      final double MINTAB = 0.1;
      final double STEP = 0.01;
      int i, j;
      double Tj;
      double Ti;
      double P;
      double H;

      // Search of the correct slot in the interpolation table
      i = (int)((x - MINTAB)/STEP + 1);
      Ti = MINTAB + i*STEP;
      Tj = Ti - STEP;

      // Approximation within the slot
      j = i - 1;
      H = x - Tj;
      R = Ti - x;
      P = STEP*STEP/6.0;
      Res = ((MWA[j]*R*R*R + MWA[i]*H*H*H)/6.0)/STEP;
      Res += ((YWA[j] - MWA[j]*P)*R + (YWA[i] - MWA[i]*P)*H)/STEP;

      // Empirical correction in 1/sqrt (n)
      Res += (CoWA[i]*H + CoWA[j]*R)/(STEP*Math.sqrt ((double)n));

      if (Res >= 1.0)
         return 1.0;
      return Res;
   }

   /**
    * Computes the complementary distribution function @f$\bar{F}_n(x)@f$
    * with parameter @f$n@f$.
    */
   public static double barF (int n, double x) {
      return 1.0 - cdf(n,x);
   }

   /**
    * Computes @f$x = F_n^{-1}(u)@f$, where @f$F_n@f$ is the *Watson*
    * @f$G@f$ distribution with parameter @f$n@f$.
    */
   public static double inverseF (int n, double u) {
      if (n <= 1)
         throw new IllegalArgumentException ("n < 2");
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u must be in [0,1]");
      if (u == 1.0)
         return Double.POSITIVE_INFINITY;
      if (u == 0.0)
         return 0.0;

      Function f = new Function (n,u);
      return RootFinder.brentDekker (0.0, 10.0, f, 1e-5);
   }

   /**
    * Returns the parameter @f$n@f$ of this object.
    */
   public int getN() {
      return n;
   }

   /**
    * Sets the parameter @f$n@f$ of this object.
    */
   public void setN (int n) {
      if (n <= 1)
         throw new IllegalArgumentException ("n < 2");
      this.n = n;
      supportA = 0.0;
      supportB = 10.0;
   }

   /**
    * Return an array containing the parameter @f$n@f$ of this object.
    */
   public double[] getParams () {
      double[] retour = {n};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : n = " + n;
   }

}