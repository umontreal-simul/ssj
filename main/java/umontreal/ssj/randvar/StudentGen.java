/*
 * Class:        StudentGen
 * Description:  Student-t random variate generators 
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
package umontreal.ssj.randvar;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;

/**
 * This class implements methods for generating random variates from the
 * *Student* distribution with @f$n>0@f$ degrees of freedom. Its density
 * function is
 * @anchor REF_randvar_StudentGen_eq_fstudent
 * @f[
 *   f (x) = \frac{\Gamma\left((n + 1)/2 \right)}{\Gamma(n/2) \sqrt{\pi n}} \left[1 + \frac{x^2}{n}\right]^{-(n+1)/2} \qquad\qquad\mbox{for } -\infty< x < \infty, \tag{fstudent}
 * @f]
 * where @f$\Gamma(x)@f$ is the gamma function defined in (
 * {@link REF_randvar_GammaGen_eq_Gamma Gamma} ).
 *
 * The `nextDouble` method simply calls `inverseF` on the distribution.
 *
 * The following table gives the CPU time needed to generate @f$10^7@f$
 * Student random variates using the different implementations available in
 * SSJ. The second test (Q) was made with the inverse in
 * @ref umontreal.ssj.probdist.StudentDistQuick, while the first test was
 * made with the inverse in  @ref umontreal.ssj.probdist.StudentDist. These
 * tests were made on a machine with processor AMD Athlon 4000, running Red
 * Hat Linux, with clock speed at 2400 MHz.
 *
 * <center>
 *
 * <table class="SSJ-table SSJ-has-hlines">
 * <tr class="bt">
 *   <td class="l bl br">Generator</td>
 *   <td class="c bl br">time in seconds</td>
 * </tr><tr class="bt">
 *   <td class="l bl br">`StudentGen`</td>
 *   <td class="c bl br">22.4</td>
 * </tr><tr>
 *   <td class="l bl br">`StudentGen(Q)`</td>
 *   <td class="c bl br">&ensp;6.5</td>
 * </tr><tr>
 *   <td class="l bl br">`StudentPolarGen`</td>
 *   <td class="c bl br">&ensp;1.4</td>
 * </tr>
 * </table>
 *
 * </center>
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class StudentGen extends RandomVariateGen {
   protected int n = -1;

   /**
    * Creates a Student random variate generator with @f$n@f$ degrees of
    * freedom, using stream `s`.
    */
   public StudentGen (RandomStream s, int n) {
      super (s, new StudentDist(n));
      setN (n);
   }

   /**
    * Creates a new generator for the Student distribution `dist` and
    * stream `s`.
    */
   public StudentGen (RandomStream s, StudentDist dist) {
      super (s, dist);
      if (dist != null)
         setN (dist.getN ());
   }

   /**
    * Generates a new variate from the Student distribution with @f$n =
    * @f$&nbsp;`n` degrees of freedom, using stream `s`.
    */
   public static double nextDouble (RandomStream s, int n) {
      return StudentDist.inverseF (n, s.nextDouble());
    }

   /**
    * Returns the value of @f$n@f$ for this object.
    */
   public int getN() {
      return n;
   }
   protected void setN (int nu) {
      if (nu <= 0)
         throw new IllegalArgumentException ("n <= 0");
      this.n = nu;
   }
}