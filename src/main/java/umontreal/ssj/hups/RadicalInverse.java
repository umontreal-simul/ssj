/*
 * Class:        RadicalInverse
 * Description:  Implements radical inverses of integers in an arbitrary basis
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
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
package umontreal.ssj.hups;

/**
 * This class implements basic methods for working with radical inverses of
 * integers in an arbitrary basis @f$b@f$. These methods are used in classes
 * that implement point sets and sequences based on the van der Corput
 * sequence (the Hammersley nets and the Halton sequence, for example).
 *
 * We recall that for a @f$k@f$-digit integer @f$i@f$ whose digital
 * @f$b@f$-ary expansion is
 * @f[
 *   i = a_0 + a_1 b + …+ a_{k-1} b^{k-1},
 * @f]
 * the *radical inverse* in base @f$b@f$ is
 * @f[
 *   \psi_b(i) = a_0 b^{-1} + a_1 b^{-2} + \cdots+ a_{k-1} b^{-k}.
 * @f]
 * The <em>van der Corput sequence in base @f$b@f$</em> is the sequence
 * @f$\psi_b(0), \psi_b(1), \psi_b(2), …@f$
 *
 * Note that @f$\psi_b(i)@f$ cannot always be represented exactly as a
 * floating-point number on the computer (e.g., if @f$b@f$ is not a power of
 * two). For an exact representation, one can use the integer
 * @f[
 *   b^k \psi_b(i) = a_{k-1} + \cdots+ a_1 b^{k-2} + a_0 b^{k-1},
 * @f]
 * which we called the *integer radical inverse* representation. This
 * representation is simply a mirror image of the digits of the usual
 * @f$b@f$-ary representation of @f$i@f$.
 *
 * It is common practice to permute locally the values of the van der Corput
 * sequence. One way of doing this is to apply a permutation to the digits of
 * @f$i@f$ before computing @f$\psi_b(i)@f$. That is, for a permutation
 * @f$\pi@f$ of the digits @f$\{0,…,b-1\}@f$,
 * @f[
 *   \psi_b(i) = \sum_{r=0}^{k-1} a_r b^{-r-1}
 * @f]
 * is replaced by
 * @f[
 *   \sum_{r=0}^{k-1} \pi(a_r) b^{-r-1}.
 * @f]
 * Applying such a permutation only changes the order in which the values of
 * @f$\psi_b(i)@f$ are enumerated. For every integer @f$k@f$, the first
 * @f$b^k@f$ values that are enumerated remain the same (they are the values
 * of @f$\psi_b(i)@f$ for @f$i=0,…,b^k-1@f$), but they are enumerated in a
 * different order. Often, different permutations @f$\pi@f$ will be applied
 * for different coordinates of a point set.
 *
 * The permutation @f$\pi@f$ can be deterministic or random. One
 * (deterministic) possibility implemented here is the Faure permutation
 * @f$\sigma_b@f$ of @f$\{0,…,b-1\}@f$ defined as follows
 * @cite rFAU92a&thinsp;. For @f$b=2@f$, take @f$\sigma= I@f$, the identical
 * permutation. For *even* @f$b=2c > 2@f$, take
 * @f{align}{
 *    \sigma[i] 
 *    & 
 *   =
 *    2\tau[i]\phantom{{} + 1} \qquad i = 0, 1, …, c-1 
 *    \\ 
 *   \sigma[i+c] 
 *    & 
 *   =
 *    2\tau[i] + 1 \qquad i = 0, 1, …, c-1
 * @f}
 * where @f$\tau[i]@f$ is the Faure permutation for base @f$c@f$. For *odd*
 * @f$b=2c+1@f$, take
 * @f{align}{
 *    \sigma[c] 
 *    & 
 *   =
 *    c 
 *    \\ 
 *   \sigma[i] 
 *    & 
 *   =
 *    \tau[i],\phantom{{} + 1} \qquad\mbox{ if } 0 \le\tau[i] < c
 *    \\ 
 *   \sigma[i] 
 *    & 
 *   =
 *    \tau[i] + 1, \qquad\mbox{ if } c \le\tau[i] < 2c
 * @f}
 * for @f$0 \le i < c@f$, and take
 * @f{align}{
 *    \sigma[i] 
 *    & 
 *   =
 *    \tau[i-1],\phantom{{} + 1} \qquad\mbox{ if } 0 \le\tau[i-1] < c
 *    \\ 
 *   \sigma[i] 
 *    & 
 *   =
 *    \tau[i-1]+1, \qquad\mbox{ if } c \le\tau[i-1] < 2c
 * @f}
 * for @f$c < i \le2c@f$, and where @f$\tau[i]@f$ is the Faure permutation
 * for base @f$c@f$. The Faure permutations give very small discrepancies
 * (amongst the best known today) for small bases.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class RadicalInverse {
   private static final int NP = 168;     // First NP primes in table below.
   private static final int PLIM = 1000;  // NP primes < PLIM

   // The first NP prime numbers
   private static final int[] PRIMES = { 
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
    71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 
    149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 
    227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 
    307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 
    389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 
    467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 
    571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 
    653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 
    751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 
    853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 
    947, 953, 967, 971, 977, 983, 991, 997 };

    // Multiplicative factors for permutations as proposed by Faure & Lemieux (2008).
    // The index corresponds to the coordinate.
    private static final int[] FAURE_LEMIEUX_FACTORS = {
      1, 1, 3, 3, 4, 9, 7, 5, 9, 18, 18, 8, 13, 31, 9, 19, 36, 33, 21, 44, 43, 
      61, 60, 56, 26, 71, 32, 77, 26, 95, 92, 47, 29, 61, 57, 69, 115, 63, 92, 
      31, 104, 126, 50, 80, 55, 152, 114, 80, 83, 97, 95, 150, 148, 55, 80, 192, 
      71, 76, 82, 109, 105, 173, 58, 143, 56, 177, 203, 239, 196, 143, 278, 227, 
      87, 274, 264, 84, 226, 163, 231, 177, 95, 116, 165, 131, 156, 105, 188, 
      142, 105, 125, 269, 292, 215, 182, 294, 152, 148, 144, 382, 194, 346, 323, 
      220, 174, 133, 324, 215, 246, 159, 337, 254, 423, 484, 239, 440, 362, 464, 
      376, 398, 174, 149, 418, 306, 282, 434, 196, 458, 313, 512, 450, 161, 315, 
      441, 549, 555, 431, 295, 557, 172, 343, 472, 604, 297, 524, 251, 514, 385, 
      531, 663, 674, 255, 519, 324, 391, 394, 533, 253, 717, 651, 399, 596, 676, 
      425, 261, 404, 691, 604, 274, 627, 777, 269, 217, 599, 447, 581, 640, 666, 
      595, 669, 686, 305, 460, 599, 335, 258, 649, 771, 619, 666, 669, 707, 737, 
      854, 925, 818, 424, 493, 463, 535, 782, 476, 451, 520, 886, 340, 793, 390, 
      381, 274, 500, 581, 345, 363, 1024, 514, 773, 932, 556, 954, 793, 294, 
      863, 393, 827, 527, 1007, 622, 549, 613, 799, 408, 856, 601, 1072, 938, 
      322, 1142, 873, 629, 1071, 1063, 1205, 596, 973, 984, 875, 918, 1133, 
      1223, 933, 1110, 1228, 1017, 701, 480, 678, 1172, 689, 1138, 1022, 682, 
      613, 635, 984, 526, 1311, 459, 1348, 477, 716, 1075, 682, 1245, 401, 774, 
      1026, 499, 1314, 743, 693, 1282, 1003, 1181, 1079, 765, 815, 1350, 1144, 
      1449, 718, 805, 1203, 1173, 737, 562, 579, 701, 1104, 1105, 1379, 827, 
      1256, 759, 540, 1284, 1188, 776, 853, 1140, 445, 1265, 802, 932, 632, 
      1504, 856, 1229, 1619, 774, 1229, 1300, 1563, 1551, 1265, 905, 1333, 493, 
      913, 1397, 1250, 612, 1251, 1765, 1303, 595, 981, 671, 1403, 820, 1404, 
      1661, 973, 1340, 1015, 1649, 855, 1834, 1621, 1704, 893, 1033, 721, 1737, 
      1507, 1851, 1006, 994, 923, 872, 1860
    };

   private static final int NRILIM = 1000;  // For nextRadicalInverse
   private int b;                     // Base
   private double invb;               // 1/b
   private double logb;               // natural log(b)
   private int JMAX;                  // b^JMAX = 2^32
   private int co;                    // Counter for nextRadicalInverse
   private double xx;                 // Current value of x
   private long ix;                   // xx = RadicalInverse (ix)
/*
   // For Struckmeier's algorithm
   private static final int PARTITION_MAX = 54; // Size of partitionL, bi
   private int partM;                 // Effective size of partitionL, bi
   private double[] bi;               // bi[i] = (b + 1)/b^i - 1
   private double[] partitionL;       // L[i] = 1 - 1/b^i
               // Boundaries of Struckmeier partitions Lkp of [0, 1]
*/

   /**
    * Initializes the base of this object to @f$b@f$ and its first value
    * of @f$x@f$ to `x0`.
    *  @param b            Base
    *  @param x0           Initial value of x
    */
   public RadicalInverse (int b, double x0) {
      co = 0;
      this.b = b;
      invb = 1.0 / b;
      logb = Math.log (b);
      JMAX = (int) (32.0 * 0.69314718055994530941 / logb);
      xx = x0;
      ix = computeI (x0);
//      initStruckmeier (b);
   }

   private long computeI (double x) {
      // Compute i such that x = RadicalInverse (i).
      int[] digits = new int[JMAX];              // Digits of x
      int j;
      for (j = 0; (j < JMAX) && (x > 0.5e-15); j++) {
         x *= b;
         digits[j] = (int) x;
         x -= digits[j];
      }
      long i = 0;
      for (j = JMAX - 1; j >= 0; j--) {
         i = i * b + digits[j];
      }
      return i;
   }

/**
 * Provides an elementary method for obtaining the first @f$n@f$ prime
 * numbers larger than 1. Creates and returns an array that contains these
 * numbers. This is useful for determining the prime bases for the different
 * coordinates of the Halton sequence and Hammersley nets.
 *  @param n            number of prime numbers to return
 *  @return an array with the first `n` prime numbers
 */
public static int[] getPrimes (int n) {
      // Allocates an array of size n filled with the first n prime
      // numbers. n must be positive (n > 0). Routine may fail if not enough
      // memory for the array is available. The first prime number is 2. 
      int i;
      boolean moreTests;
      int[] prime = new int[n];

      int n1 = Math.min (NP, n);
      for (i = 0; i < n1; i++)
         prime[i] = PRIMES[i];
      if (NP < n) {
         i = NP;
         for (int candidate = PLIM + 1; i < n; candidate += 2) {
             prime[i] = candidate;
             for (int j = 1; (moreTests = prime[j] <= candidate / prime[j])
                     && ((candidate % prime[j]) > 0); j++);
             if (! moreTests)
                 i++;
         }
      }
      return prime;
   }

   /**
    * Computes the radical inverse of @f$i@f$ in base @f$b@f$. If
    * @f$i=\sum_{r=0}^{k-1} a_r b^r@f$, the method computes and returns
    * @f[
    *   x = \sum_{r=0}^{k-1} a_r b^{-r-1}.
    * @f]
    * @param b            base used for the operation
    *  @param i            the value for which the radical inverse will be
    *                      computed
    *  @return the radical inverse of `i` in base `b`
    */
   public static double radicalInverse (int b, long i) {
      double digit, radical, inverse;
      digit = radical = 1.0 / (double) b;
      for (inverse = 0.0; i > 0; i /= b) {
         inverse += digit * (double) (i % b);
         digit *= radical;
      }
      return inverse;
   }

   /**
    * Computes the radical inverse of @f$x@f$ in base @f$b@f$. If @f$x@f$
    * has more decimals in base @f$b@f$ than
    * @f$\log_b@f$(<tt>Long.MAX_VALUE</tt>), it is truncated to its
    * minimum precision in base @f$b@f$. If @f$x=\sum_{r=0}^{k-1} a_r
    * b^{-r-1}@f$, the method computes and returns
    * @f[
    *   i = \sum_{r=0}^{k-1} a_r b^r.
    * @f]
    * @param b            base used for the operation
    *  @param x            the value for which the radical inverse will be
    *                      computed
    *  @return the radical inverse of `x` in base `b`
    */
   public static int radicalInverseInteger (int b, double x) {
      int digit = 1;
      int inverse = 0;
      int precision = Integer.MAX_VALUE / (2 * b * b);
      while (x > 0 && inverse < precision) {
        int p = digit * b;
        double y = Math.floor(x * p);
        inverse += digit * (int)y;
        x -= y / (double)p;
        digit *= b;
      }
      return inverse;
   }

   /**
    * Computes the radical inverse of @f$x@f$ in base @f$b@f$. If @f$x@f$
    * has more decimals in base @f$b@f$ than
    * @f$\log_b@f$(<tt>Long.MAX_VALUE</tt>), it is truncated to its
    * minimum precision in base @f$b@f$. If @f$x=\sum_{r=0}^{k-1} a_r
    * b^{-r-1}@f$, the method computes and returns
    * @f[
    *   i = \sum_{r=0}^{k-1} a_r b^r.
    * @f]
    * @param b            base used for the operation
    *  @param x            the value for which the radical inverse will be
    *                      computed
    *  @return the radical inverse of `x` in base `b`
    */
   public static long radicalInverseLong (int b, double x) {
      long digit = 1;
      long inverse = 0;
      long precision = Long.MAX_VALUE / (b * b * b);
      while (x > 0 && inverse < precision) {
        long p = digit * b;
        double y = Math.floor(x * p);
        inverse += digit * (long)y;
        x -= y / (double)p;
        digit *= b;
      }
      return inverse;
   }

   /**
    * A fast method that incrementally computes the radical inverse
    * @f$x_{i+1}@f$ in base @f$b@f$ from @f$x_i@f$ = `x` =
    * @f$\psi_b(i)@f$, using addition with *rigthward carry*. The
    * parameter `invb` is equal to @f$1/b@f$. Using long incremental
    * streams (i.e., calling this method several times in a row) cause
    * increasing inaccuracy in @f$x@f$. Thus the user should recompute the
    * radical inverse directly by calling  #radicalInverse every once in a
    * while (i.e. in every few thousand calls).
    *  @param invb         @f$1/b@f$ where @f$b@f$ is the base
    *  @param x            the inverse @f$x_i@f$
    *  @return the radical inverse @f$x_{i+1}@f$
    */
   public static double nextRadicalInverse (double invb, double x) {
      // Calculates the next radical inverse from x in base b.
      // Repeated application causes a loss of accuracy.
      // Note that x can be any number from [0,1).

      final double ALMOST_ONE = 1.0 - 1e-10;
      double nextInverse = x + invb;
      if (nextInverse < ALMOST_ONE)
         return nextInverse;
      else {
         double digit1 = invb;
         double digit2 = invb * invb;
         while (x + digit2 >= ALMOST_ONE) {
            digit1 = digit2;
            digit2 *= invb;
         }
         return x + (digit1 - 1.0) + digit2;
      }
   }

   /**
    * A fast method that incrementally computes the radical inverse
    * @f$x_{i+1}@f$ in base @f$b@f$ from @f$x_i@f$ = @f$\psi_b(i)@f$,
    * using addition with *rigthward carry* as described in
    * @cite vWAN99a&thinsp;. Since using long incremental streams (i.e.,
    * calling this method several times in a row) cause increasing
    * inaccuracy in @f$x@f$, the method recomputes the radical inverse
    * directly from @f$i@f$ by calling  #radicalInverse once in every 1000
    * calls.
    *  @return the radical inverse @f$x_{i+1}@f$
    */
   public double nextRadicalInverse () {
      // Calculates the next radical inverse from xx in base b.
      // Repeated application causes a loss of accuracy.
      // For each NRILIM calls, a direct calculation via radicalInverse
      // is inserted.

      co++;
      if (co >= NRILIM) {
         co = 0;
         ix += NRILIM;
         xx = radicalInverse (b, ix);
         return xx;
      }
      final double ALMOST_ONE = 1.0 - 1e-10;
      double nextInverse = xx + invb;
      if (nextInverse < ALMOST_ONE) {
         xx = nextInverse;
         return xx;
      } else {
         double digit1 = invb;
         double digit2 = invb * invb;
         while (xx + digit2 >= ALMOST_ONE) {
            digit1 = digit2;
            digit2 *= invb;
         }
         xx += (digit1 - 1.0) + digit2;
         return xx;
      }
   }

/**
 * A fast method that incrementally computes the radical inverse
 * @f$x_{i+1}@f$ in base @f$b@f$ from @f$x_i =@f$ ` x`, using the
 * Struckmeier’s algorithm described in @cite rSTR95a&thinsp;. It uses a
 * small precomputed table of values @f$L_k = 1 - b^{-k}@f$. The method
 * returns the next radical inverse @f$x_{i+1} = x_i + (b + 1 - b^k) /
 * b^k@f$, where @f$L_{k-1} \le x < L_k@f$.
 * @remark **Richard:** This method can work only if it is reprogrammed with
 * integers. With floating-point numbers, unavoidable accumulating rounding
 * errors will sooner or later lead to choosing the wrong interval, after
 * which, all subsequent x’s will be completely wrong.
 *  @param b            base
 *  @param x            the inverse @f$x_i@f$
 *  @return the radical inverse @f$x_{i+1}@f$
 */
/*
   private void initStruckmeier (int b) {
      bi = new double[1 + PARTITION_MAX];
      partitionL = new double[1 + PARTITION_MAX];
      logb = Math.log (b);
      partitionL[0] = 0.0;
      bi[0] = 1.0;
      int i = 0;
      while ((i < PARTITION_MAX) && (partitionL[i] < 1.0)) {
         ++i;
         bi[i] = bi[i - 1] / b;
         partitionL[i] = 1.0 - bi[i];
      }
      partM = i - 1;

      for (i = 0; i <= partM + 1; ++i)
         bi[i] = (b + 1) * bi[i] - 1.0;
   }
*/

/*
   public double nextRadicalInverse (double x) {
      int k;
      if (x < partitionL[partM]) {
         k = 1;
         // Find k such:    partitionL[k-1] <= x < partitionL[k]  
         while (x >= partitionL[k])
            ++k;

      } else {           // x >= partitionL [partM]
         k = 1 + (int)(-Math.log(1.0 - x) / logb);
      }
      return x + bi[k];
   }
*/

   /**
    * Given the @f$k@f$ @f$b@f$-ary digits of @f$i@f$ in `bdigits`,
    * returns the @f$k@f$ digits of the integer radical inverse of @f$i@f$
    * in `idigits`. This simply reverses the order of the digits.
    *  @param k            number of digits in arrays
    *  @param bdigits      digits in original order
    *  @param idigits      digits in reverse order
    */
   public static void reverseDigits (int k, int bdigits[], int idigits[]) {
      for (int l = 0; l < k; l++)
         idigits[l] = bdigits[k-l];
   }

   /**
    * Computes the integer radical inverse of @f$i@f$ in base @f$b@f$,
    * equal to @f$b^k \psi_b(i)@f$ if @f$i@f$ has @f$k@f$ @f$b@f$-ary
    * digits.
    *  @param b            base used for the operation
    *  @param i            the value for which the integer radical inverse
    *                      will be computed
    *  @return the integer radical inverse of `i` in base `b`
    */
   public static int integerRadicalInverse (int b, int i) {
      // Simply flips digits of i in base b.
      int inverse;
      for (inverse = 0; i > 0; i /= b)
         inverse = inverse * b + (i % b);
      return inverse;
   }

   /**
    * Given the @f$k@f$ digits of the integer radical inverse of @f$i@f$
    * in `bdigits`, in base @f$b@f$, this method replaces them by the
    * digits of the integer radical inverse of @f$i+1@f$ and returns their
    * number. The array must be large enough to hold this new number of
    * digits.
    *  @param b            base
    *  @param k            initial number of digits in arrays
    *  @param idigits      digits of integer radical inverse
    *  @return new number of digits in arrays
    */
   public static int nextRadicalInverseDigits (int b, int k, int idigits[]) {
      int l;
      for (l = k-1; l >= 0; l--)
         if (idigits[l] == b-1) 
            idigits[l] = 0;
         else {
            idigits[l]++;
            return k;
         }
      if (l == 0) {
         idigits[k] = 1;
         return ++k;
      }
      return 0;
   }

   /**
    * Computes the permutations as proposed in @cite vFAU09a&thinsp;
    * @f$\sigma_b@f$ of the set @f$\{0, …, b - 1\}@f$ and puts it in
    * array `pi`.
    *  @param coordinate   the coordinate
    *  @param pi           an array of size at least `b`, to be filled
    *                      with the permutation
    */
   public static void getFaureLemieuxPermutation (int coordinate, int[] pi) {
      int f = FAURE_LEMIEUX_FACTORS[coordinate];
      int b = PRIMES[coordinate];
      for (int k = 0; k < pi.length; k++)
         pi[k] = f * k % b;
   }

   /**
    * Computes the Faure permutation @cite rFAU92a&thinsp; @f$\sigma_b@f$
    * of the set @f$\{0, …, b - 1\}@f$ and puts it in array `pi`. See the
    * description in the introduction above.
    *  @param b            the base
    *  @param pi           an array of size at least `b`, to be filled
    *                      with the permutation
    */
   public static void getFaurePermutation (int b, int[] pi) {
      // This is a recursive implementation.  
      // Perhaps not the most efficient...
      int i;
      if (b == 2) {
         pi[0] = 0;
         pi[1] = 1;
      }
      else if ((b & 1) != 0) {
         // b is odd.
         b--;
         getFaurePermutation (b, pi);
         for (i = 0; i < b; i++)
            if (pi[i] >= b / 2)
               pi[i]++;
         for (i = b; i > b / 2; i--)
            pi[i] = pi[i - 1];
         pi[b / 2] = b / 2;
      }
      else {
         b /= 2;
         getFaurePermutation (b, pi);
         for (i = 0; i < b; i++) {
            pi[i] *= 2;
            pi[i + b] = pi[i] + 1;
         }
      }
   }

   /**
    * Computes the radical inverse of @f$i@f$ in base @f$b@f$, where the
    * digits are permuted using the permutation @f$\pi@f$. If
    * @f$i=\sum_{r=0}^{k-1} a_r b^r@f$, the method will compute and
    * return
    * @f[
    *   x = \sum_{r=0}^{k-1} \pi[a_r] b^{-r-1}.
    * @f]
    * @param b            base @f$b@f$ used for the operation
    *  @param pi           an array of length at least `b` containing the
    *                      permutation used during the computation
    *  @param i            the value for which the radical inverse will be
    *                      computed
    *  @return the radical inverse of `i` in base `b`
    */
   public static double permutedRadicalInverse (int b, int[] pi, long i) {
      double digit, radical, inverse;
      digit = radical = 1.0 / (double) b;
      for (inverse = 0.0; i > 0; i /= b) {
         inverse += digit * (double) pi[(int)(i % b)];
         digit *= radical;
      }
      return inverse;
   }

}