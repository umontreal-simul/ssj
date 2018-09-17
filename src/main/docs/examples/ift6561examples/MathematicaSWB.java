package ift6561examples;
import umontreal.ssj.rng.RandomStreamBase;


public class MathematicaSWB extends RandomStreamBase {
   private final int r = 8;
   private final int s = 48;
   protected final double norm = 1.0 / (1L << 31); // 1 / 2^31
   private static final int MASK31 = 0x7fffffff; // Masque de 31 bits
   // static final long serialVersionUID;

   private int[] X;
   private int m_c;
   private int m_i;
   private int m_j;

   // =======================================================================

   public MathematicaSWB() {
      X = new int[s];

      // Initialiser l'état du générateur avec des valeurs arbitraires
      X[0] = 12345;
      for (int i = 1; i < s; i++)
         X[i] = (16807 * X[i - 1]) & MASK31;

      m_c = 0;
      m_i = 0;
      m_j = s - r;
   }

   // -------------------------------------------------------------------------

   public void setSeed(int[] seed, int c) {
      for (int i = 0; i < s; i++)
         X[i] = seed[i] & MASK31;
      m_c = c;
      m_i = 0;
      m_j = s - r;
   }

   // =======================================================================

   public double nextValue() {
      int v = X[m_j] - X[m_i] - m_c; // ne peut pas déborder

      if (v < 0) {
         v &= MASK31; // équivalent à v = v + 2^31
         m_c = 1;
      } else {
         m_c = 0;
      }

      X[m_i] = v;

      m_j++;
      if (m_j == s)
         m_j = 0;

      m_i++;
      if (m_i == s)
         m_i = 0;

      return v * norm;
   }

   // =======================================================================

   public double nextDouble() {
      // Retourne un nombre réel dans [0, 1)
      double u = nextValue();
      return u * norm + nextValue();
   }

   // =======================================================================

   public int nextInt(int a, int b) {
      // Retourne un entier dans [a, b]
      return a + (int) (nextDouble() * (b - a + 1.0));
   }

   // =======================================================================

   public String toString() {
      return "MathematicaSWB";
   }

   public void resetStartStream() {
      throw new UnsupportedOperationException ();
   }

   public void resetStartSubstream() {
      throw new UnsupportedOperationException ();
   }

   public void resetNextSubstream() {
      throw new UnsupportedOperationException ();
   }

   //=======================================================================

   public static void main(String[] args) {
      double x;
      MathematicaSWB gen = new MathematicaSWB();

      for (int i = 0; i < 20; i++) {
         x = gen.nextDouble();
         System.out.printf("%12.7f%n", x);
      }
   }

}
