package umontreal.ssj.markovchainrqmc;

import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.hups.*;
import java.util.*;

/**
 * DEPRECATED: Should use a  @ref umontreal.ssj.hups.SobolSequence with a
 * @ref umontreal.ssj.hups.LMScrambleShift.
 *
 * A Sobol sequence randomized by a left matrix scramble followed by a
 * digital random shift.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
@Deprecated
public class LeftScrambledSobolSequence extends SobolSequence {

   /**
    * Same as <tt>SobolSequence(k, w, dim)</tt>, except that its
    * #randomize method will do a left matrix scramble followed by a
    * random digital shift.
    */
   public LeftScrambledSobolSequence (int k, int w, int dim) {
       super (k, w, dim);
   }


   public void addRandomShift ()  {
      leftMatrixScramble (shiftStream);
      addRandomShift (0, dimShift, shiftStream);
   }

   public String toString() {
      StringBuffer sb = new StringBuffer ("LeftScrambledSobolNet:" +
          PrintfFormat.NEWLINE);
      sb.append (super.toString());
      return sb.toString();
   }
}