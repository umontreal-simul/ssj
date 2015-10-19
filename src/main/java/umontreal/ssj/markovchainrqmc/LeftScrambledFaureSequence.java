package umontreal.ssj.markovchainrqmc;

import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.hups.*;
import java.util.*;

/**
 * DEPRECATED: Should use a  @ref umontreal.ssj.hups.FaureSequence with a
 * @ref umontreal.ssj.hups.LMScrambleShift.
 *
 * A Faure sequence randomized by a left matrix scramble followed by a
 * digital random shift.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
@Deprecated
public class LeftScrambledFaureSequence extends FaureSequence {

   /**
    * Same as <tt>FaureSequence(b, k, r, w, dim)</tt>, except that its
    * #randomize method will do a left matrix scramble followed by a
    * random digital shift.
    */
   public LeftScrambledFaureSequence (int b, int k, int r, int w, int dim) {
       super (b, k, r, w, dim);
   }


   public void randomize (RandomStream noise) {
      leftMatrixScramble (noise);
      addRandomShift (noise);
   }

   public void randomize (int d1, int d2, RandomStream noise)  {
      leftMatrixScramble (noise);
      addRandomShift (d1, d2, noise);
   }

   public String toString() {
      StringBuffer sb = new StringBuffer ("LeftScrambledFaureNet:" +
          PrintfFormat.NEWLINE);
      sb.append (super.toString());
      return sb.toString();
   }
}