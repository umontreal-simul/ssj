package ift6561examples;
import java.io.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.hups.*;
import umontreal.ssj.util.*;


/**
 * Reads a set of points, applied a digital random shift to them, then prints
 * the shifted points in a file.
 */

public class DigitalShift
{
   private String formatBase (int base, double x) {
      if (x >= 1)
         return PrintfFormat.formatBase(base, (long)x);
      else {
         StringBuffer sb = new StringBuffer("0.");
         long y;
         for (int j = 0; j < 20; ++j) {
            x *= 2;
            y = (long) x;
            if (y == 1)
               sb.append("1");
            else
               sb.append("0");
            x -= y;
         }
         return sb.toString();
      }
   }


   private String formatPoints (PointSet set, int n, int d, String mess) {
      return formatPoints (set, n, d, 10, mess);
   }

   private String formatPoints (PointSet set, int n, int d, int base, String mess) {
      if (set.getNumPoints() < n)
         n = set.getNumPoints();
      if (set.getDimension() < d)
         d = set.getDimension();
      StringBuffer sb = new StringBuffer(mess);
      sb.append (":" + PrintfFormat.LINE_SEPARATOR);
      PointSetIterator itr = set.iterator();
      for (int i = 0; i < n; i++) {
         for (int j = 0; j < d; j++) {
            sb.append ("  ");
            if (base == 10)
               sb.append (itr.nextCoordinate());
            else
               sb.append (formatBase (base, itr.nextCoordinate()));
         }
         sb.append (PrintfFormat.LINE_SEPARATOR);
         itr.resetToNextPoint();
      }
      return sb.toString();
   }


   public DigitalShift(int n) {
      DigitalSequenceBase2 p0 = new SobolSequence(n, 1);
      DigitalNetBase2 p = p0.toNetShiftCj();
      int dim = 2;
      System.out.println(formatPoints(p, n, dim, "Sobol points"));
      System.out.println();

      RandomStream stream = new MRG32k3a();
      p.addRandomShift(0, 1, stream);
      System.out.println(formatPoints(p, n, dim, "Sobol shifted"));
      System.out.println();

      System.out.println(formatPoints(p, n, dim, 2, "Sobol shifted, base = 2"));
      System.out.println();
   }


   public static void main(String[] args) throws IOException {
      new DigitalShift(16);

   }
}
