/*
 * Class:        ArithmeticMod
 * Description:  multiplications of scalars, vectors and matrices modulo m
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
package umontreal.ssj.util;

/**
 * This class provides facilities to compute multiplications of scalars, of
 * vectors and of matrices modulo m. All algorithms are present in three
 * different versions. These allow operations on `double`, `int` and `long`.
 * The `int` and `long` versions work exactly like the `double` ones.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class ArithmeticMod {

   //private constants
   private static final double two17    =  131072.0;
   private static final double two53    =  9007199254740992.0;

   //prevent the creation of the object
   private ArithmeticMod() {};

   /**
    * @name Methods using `double`
    * @{
    */

   /**
    * Computes @f$(a \times s + c) \bmod m@f$. Where `m` must be smaller
    * than @f$2^{35}@f$. Works also if `s` or `c` are negative. The result
    * is always positive (and thus always between 0 and `m` - 1).
    *  @param a            the first factor of the multiplication
    *  @param s            the second factor of the multiplication
    *  @param c            the second term of the addition
    *  @param m            the modulus
    *  @return the result of the multiplication and the addition modulo
    * `m`
    */
   public static double multModM (double a, double s, double c, double m) {
      int a1;
      double v = a * s + c;
      if (v >= two53 || v <= -two53 ) {
         a1 = (int)(a / two17);
         a -= a1 * two17;
         v  = a1 * s;
         a1 = (int)(v / m);
         v -= a1 * m;
         v  = v * two17 + a * s + c;
      }
      a1 = (int)(v / m);
      if ((v -= a1 * m) < 0.0)
         return v += m;
      else
         return v;
   }

   /**
    * Computes the result of @f$\mathtt{A} \times\mathtt{\mathbf{s}}
    * \bmod m@f$ and puts the result in `v`. Where `s` and `v` are both
    * column vectors. This method works even if `s` = `v`.
    *  @param A            the multiplication matrix
    *  @param s            the multiplied vector
    *  @param v            the result of the multiplication
    *  @param m            the modulus
    */
   public static void matVecModM (double A[][], double s[], double v[],
                                  double m) {
      int i;
      double x[] = new double[v.length];
      for (i = 0; i < v.length;  ++i) {
         x[i] = 0.0;
         for(int j = 0; j < s.length; j++)
            x[i] = multModM (A[i][j], s[j], x[i], m);
      }
      for (i = 0; i < v.length;  ++i)
         v[i] = x[i];
   }

   /**
    * Computes @f$\mathtt{A} \times\mathtt{B} \bmod m@f$ and puts the
    * result in `C`. Works even if `A` = `C`, `B` = `C` or `A` = `B` =
    * `C`.
    *  @param A            the first factor of the multiplication
    *  @param B            the second factor of the multiplication
    *  @param C            the result of the multiplication
    *  @param m            the modulus
    */
   public static void matMatModM (double A[][], double B[][], double C[][],
                                  double m) {
      int i, j;
      int r = C.length;    //# of rows of C
      int c = C[0].length; //# of columns of C
      double V[] = new double[r];
      double W[][] = new double[r][c];
      for (i = 0; i < c;  ++i) {
         for (j = 0; j < r;  ++j)
            V[j] = B[j][i];
         matVecModM (A, V, V, m);
         for (j = 0; j < r;  ++j)
            W[j][i] = V[j];
      }
      for (i = 0; i < r;  ++i) {
         for (j = 0; j < c;  ++j)
            C[i][j] = W[i][j];
      }
   }

   /**
    * Computes @f$\mathtt{A}^{2^{\mathtt{e}}} \bmod m@f$ and puts the
    * result in `B`. Works even if `A` = `B`.
    *  @param A            the matrix to raise to a power
    *  @param B            the result of exponentiation
    *  @param m            the modulus
    *  @param e            the @f$\log_2@f$ of the exponent
    */
   public static void matTwoPowModM (double A[][], double B[][], double m,
                                     int e) {
      int i, j;
      /* initialize: B = A */
      if (A != B) {
         for (i = 0; i < A.length; i++) {
            for (j = 0; j < A.length;  ++j)  //A is square
               B[i][j] = A[i][j];
         }
      }
      /* Compute B = A^{2^e} */
      for (i = 0; i < e; i++)
         matMatModM (B, B, B, m);
   }

   /**
    * Computes @f$\mathtt{A}^c \bmod m@f$ and puts the result in `B`.
    * Works even if `A` = `B`.
    *  @param A            the matrix to raise to a power
    *  @param B            the result of the exponentiation
    *  @param m            the modulus
    *  @param c            the exponent
    */
   public static void matPowModM (double A[][], double B[][], double m,
                                  int c) {
      int i, j;
      int n = c;
      int s = A.length;   //we suppose that A is square
      double W[][] = new double[s][s];

      /* initialize: W = A; B = I */
      for (i = 0; i < s; i++) {
         for (j = 0; j < s;  ++j)  {
            W[i][j] = A[i][j];
            B[i][j] = 0.0;
         }
      }
      for (j = 0; j < s;  ++j)
         B[j][j] = 1.0;

      /* Compute B = A^c mod m using the binary decomp. of c */
      while (n > 0) {
         if ((n % 2)==1)
            matMatModM (W, B, B, m);
         matMatModM (W, W, W, m);
         n /= 2;
      }
   }

   /**
    * @}
    */

   /**
    * @name Methods using `int`
    * @{
    */

   /**
    * Computes @f$(a \times s + c) \bmod m@f$. Works also if `s` or `c`
    * are negative. The result is always positive (and thus always between
    * 0 and `m` - 1).
    *  @param a            the first factor of the multiplication
    *  @param s            the second factor of the multiplication
    *  @param c            the second term of the addition
    *  @param m            the modulus
    *  @return the result of the multiplication and the addition modulo
    * `m`
    */
   public static int multModM (int a, int s, int c, int m) {
      int r = (int) (((long)a * s + c) % m);
      return r < 0 ? r + m : r;
   }

   /**
    * Exactly like  #matVecModM(double[][],double[],double[],double) using
    * `double`, but with `int` instead of `double`.
    *  @param A            the multiplication matrix
    *  @param s            the multiplied vector
    *  @param v            the result of the multiplication
    *  @param m            the modulus
    */
   public static void matVecModM (int A[][], int s[], int v[], int m) {
      int i;
      int x[] = new int[v.length];
      for (i = 0; i < v.length;  ++i) {
         x[i] = 0;
         for(int j = 0; j < s.length; j++)
            x[i] = multModM(A[i][j], s[j], x[i], m);
      }
      for (i = 0; i < v.length;  ++i)
         v[i] = x[i];

   }

   /**
    * Exactly like  #matMatModM(double[][],double[][],double[][],double)
    * using `double`, but with `int` instead of `double`.
    *  @param A            the first factor of the multiplication
    *  @param B            the second factor of the multiplication
    *  @param C            the result of the multiplication
    *  @param m            the modulus
    */
   public static void matMatModM (int A[][], int B[][], int C[][], int m) {
      int i, j;
      int r = C.length;    //# of rows of C
      int c = C[0].length; //# of columns of C
      int V[] = new int[r];
      int W[][] = new int[r][c];
      for (i = 0; i < c;  ++i) {
         for (j = 0; j < r;  ++j)
            V[j] = B[j][i];
         matVecModM (A, V, V, m);
         for (j = 0; j < r;  ++j)
            W[j][i] = V[j];
      }
      for (i = 0; i < r;  ++i) {
         for (j = 0; j < c;  ++j)
            C[i][j] = W[i][j];
      }
   }

   /**
    * Exactly like  #matTwoPowModM(double[][],double[][],double,int) using
    * `double`, but with `int` instead of `double`.
    *  @param A            the matrix to raise to a power
    *  @param B            the result of exponentiation
    *  @param m            the modulus
    *  @param e            the @f$\log_2@f$ of the exponent
    */
   public static void matTwoPowModM (int A[][], int B[][], int m, int e) {
      int i, j;
      /* initialize: B = A */
      if (A != B) {
         for (i = 0; i < A.length; i++) {
            for (j = 0; j < A.length;  ++j)  //A is square
               B[i][j] = A[i][j];
         }
      }
      /* Compute B = A^{2^e} */
      for (i = 0; i < e; i++)
         matMatModM (B, B, B, m);
   }

   /**
    * Exactly like  #matPowModM(double[][],double[][],double,int) using
    * `double`, but with `int` instead of `double`.
    *  @param A            the matrix to raise to a power
    *  @param B            the result of the exponentiation
    *  @param m            the modulus
    *  @param c            the exponent
    */
   public static void matPowModM (int A[][], int B[][], int m, int c) {
      int i, j;
      int n = c;
      int s = A.length;   //we suppose that A is square (it must be)
      int W[][] = new int[s][s];

      /* initialize: W = A; B = I */
      for (i = 0; i < s; i++) {
         for (j = 0; j < s;  ++j)  {
            W[i][j] = A[i][j];
            B[i][j] = 0;
         }
      }
      for (j = 0; j < s;  ++j)
         B[j][j] = 1;

      /* Compute B = A^c mod m using the binary decomp. of c */
      while (n > 0) {
         if ((n % 2)==1)
            matMatModM (W, B, B, m);
         matMatModM (W, W, W, m);
         n /= 2;
      }
   }

   /**
    * @}
    */

   /**
    * @name Methods using `long`
    * @{
    */

   /**
    * Computes @f$(a \times s + c) \bmod m@f$. Works also if `s` or `c`
    * are negative. The result is always positive (and thus always between
    * 0 and `m` - 1).
    *  @param a            the first factor of the multiplication
    *  @param s            the second factor of the multiplication
    *  @param c            the second term of the addition
    *  @param m            the modulus
    *  @return the result of the multiplication and the addition modulo
    * `m`
    */
   public static long multModM (long a, long s, long c, long m) {

      /* Suppose que 0 < a < m  et  0 < s < m.   Retourne (a*s + c) % m.
       * Cette procedure est tiree de :
       * L'Ecuyer, P. et Cote, S., A Random Number Package with
       * Splitting Facilities, ACM TOMS, 1991.
       * On coupe les entiers en blocs de d bits. H doit etre egal a 2^d.  */

      final long H = 2147483648L;               // = 2^d  used in MultMod
      long a0, a1, q, qh, rh, k, p;
      if (a < H) {
         a0 = a;
         p = 0;
      } else {
         a1 = a / H;
         a0 = a - H * a1;
         qh = m / H;
         rh = m - H * qh;
         if (a1 >= H) {
            a1 = a1 - H;
            k = s / qh;
            p = H * (s - k * qh) - k * rh;
            if (p < 0)
               p = (p + 1) % m + m - 1;
         } else                      /* p = (A2 * s * h) % m.      */
            p = 0;
         if (a1 != 0) {
            q = m / a1;
            k = s / q;
            p -= k * (m - a1 * q);
            if (p > 0)
               p -= m;
            p += a1 * (s - k * q);
            if (p < 0)
               p = (p + 1) % m + m - 1;
         }                           /* p = ((A2 * h + a1) * s) % m. */
         k = p / qh;
         p = H * (p - k * qh) - k * rh;
         if (p < 0)
            p = (p + 1) % m + m - 1;
      }                               /* p = ((A2 * h + a1) * h * s) % m  */
      if (a0 != 0) {
         q = m / a0;
         k = s / q;
         p -= k * (m - a0 * q);
         if (p > 0)
            p -= m;
         p += a0 * (s - k * q);
         if (p < 0)
            p = (p + 1) % m + m - 1;
      }
      p = (p - m) + c;
      if (p < 0)
         p += m;
      return p;
   }

   /**
    * Exactly like  #matVecModM(double[][],double[],double[],double) using
    * `double`, but with `long` instead of `double`.
    *  @param A            the multiplication matrix
    *  @param s            the multiplied vector
    *  @param v            the result of the multiplication
    *  @param m            the modulus
    */
   public static void matVecModM (long A[][], long s[], long v[], long m) {
      int i;
      long x[] = new long[v.length];
      for (i = 0; i < v.length;  ++i) {
         x[i] = 0;
         for(int j = 0; j < s.length; j++)
            x[i] = multModM(A[i][j], s[j], x[i], m);
      }
      for (i = 0; i < v.length;  ++i)
         v[i] = x[i];

   }

   /**
    * Exactly like  #matMatModM(double[][],double[][],double[][],double)
    * using `double`, but with `long` instead of `double`.
    *  @param A            the first factor of the multiplication
    *  @param B            the second factor of the multiplication
    *  @param C            the result of the multiplication
    *  @param m            the modulus
    */
   public static void matMatModM (long A[][], long B[][], long C[][], long m) {
      int i, j;
      int r = C.length;    //# of rows of C
      int c = C[0].length; //# of columns of C
      long V[] = new long[r];
      long W[][] = new long[r][c];
      for (i = 0; i < c;  ++i) {
         for (j = 0; j < r;  ++j)
            V[j] = B[j][i];
         matVecModM (A, V, V, m);
         for (j = 0; j < r;  ++j)
            W[j][i] = V[j];
      }
      for (i = 0; i < r;  ++i) {
         for (j = 0; j < c;  ++j)
            C[i][j] = W[i][j];
      }
   }

   /**
    * Exactly like  #matTwoPowModM(double[][],double[][],double,int) using
    * `double`, but with `long` instead of `double`.
    *  @param A            the matrix to raise to a power
    *  @param B            the result of exponentiation
    *  @param m            the modulus
    *  @param e            the @f$\log_2@f$ of the exponent
    */
   public static void matTwoPowModM (long A[][], long B[][], long m, int e) {
      int i, j;
      /* initialize: B = A */
      if (A != B) {
         for (i = 0; i < A.length; i++) {
            for (j = 0; j < A.length;  ++j)  //A is square
               B[i][j] = A[i][j];
         }
      }
      /* Compute B = A^{2^e} */
      for (i = 0; i < e; i++)
         matMatModM (B, B, B, m);
   }

   /**
    * Exactly like  #matPowModM(double[][],double[][],double,int) using
    * `double`, but with `long` instead of `double`.
    *  @param A            the matrix to raise to a power
    *  @param B            the result of the exponentiation
    *  @param m            the modulus
    *  @param c            the exponent
    */
   public static void matPowModM (long A[][], long B[][], long m, int c) {
      int i, j;
      int n = c;
      int s = A.length;   //we suppose that A is square (it must be)
      long W[][] = new long[s][s];

      /* initialize: W = A; B = I */
      for (i = 0; i < s; i++) {
         for (j = 0; j < s;  ++j)  {
            W[i][j] = A[i][j];
            B[i][j] = 0;
         }
      }
      for (j = 0; j < s;  ++j)
         B[j][j] = 1;

      /* Compute B = A^c mod m using the binary decomp. of c */
      while (n > 0) {
         if ((n % 2)==1)
            matMatModM (W, B, B, m);
         matMatModM (W, W, W, m);
         n /= 2;
      }
   }

}

/**
 * @}
 */