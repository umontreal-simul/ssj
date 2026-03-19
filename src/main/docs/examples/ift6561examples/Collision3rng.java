package ift6561examples;

import umontreal.ssj.rng.*;

// We throw m points in k=d^3 locations in a 3-dim cube and count the number 
// C of collisions. The points are produced by a RNG, with specified lags.
// They have the form \bu_i = (u_ji, u_{ji+i1}, u_{ji+i2}) 
// where $j = i2+1$ and i2 > i1 > 0 are fixed integers.

public class Collision3rng extends Collision {
   int d; // Number of divisions in each dimension.
   int[] point; // To store integer coordinates of point \bu_i
   int i1, i2;

   public Collision3rng(int d, int m, int i1, int i2) {
      super(d * d * d, m, 1);
      this.d = d;
      this.i1 = i1;
      this.i2 = i2;
      point = new int[i2 + 1];
   }

   // Generates the number of collisions, for fixed m.
   public int simulate(RandomStream stream) {
      C = 0;
      int address;
      for (int i = 0; i < k; i++)
         used[i] = false;
      for (int j = 0; j < m; j++) {
         stream.nextArrayOfInt(0, d - 1, point, 0, i2 + 1);
         address = (point[0] * d + point[i1]) * d + point[i2];
         if (used[address] == true)
            C++;
         else
            used[address] = true;
      }
      return C;
   }

   public String toString() {
      String s = "\nCollisions in a hashing system\n" + "k = " + k + " locations \n" + "m = " + m + " items \n"
            + "Theorical mean = lambda = " + lambda + "\n" + " i1 = " + i1 + "  i2 = " + i2 + "\n";
      return s;
   }

   // Performs n indep. runs using stream
   public void simulateRuns(int n, RandomStream stream) {
      System.out.println("\n" + stream.toString());
      for (int i = 0; i < n; i++)
         System.out.printf("Nombre de collisions =  %d%n", simulate(stream));
   }

   public static void main(String[] args) {
      int d = 100;
      int m = 10000;

      Collision3rng col = new Collision3rng(d, m, 1, 2);
      System.out.println(col.toString());
      col.simulateRuns(10, new MRG32k3a());
      col.simulateRuns(10, new LFSR113());
      col.simulateRuns(10, new MathematicaSWB());

      col = new Collision3rng(d, m, 20, 24);
      System.out.println(col.toString());
      col.simulateRuns(10, new MRG32k3a());
      col.simulateRuns(10, new LFSR113());
      col.simulateRuns(10, new MathematicaSWB());
   }
}
