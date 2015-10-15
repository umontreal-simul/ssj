/*
 * Class:        UniformIntGen
 * Description:  random variate generator for the uniform distribution over integers
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since

 * SSJ is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, either version 3 of the License, or
 * any later version.

 * SSJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * A copy of the GNU General Public License is available at
   <a href="http://www.gnu.org/licenses">GPL licence site</a>.
 */
package umontreal.ssj.randvar;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;

/**
 * This class implements a random variate generator for the *uniform*
 * distribution over integers, over the interval @f$[i,j]@f$. Its mass
 * function is
 * @anchor REF_randvar_UniformIntGen_eq_fmassuniformint
 * @f[
 *   p(x) = \frac{1}{j - i + 1} \qquad\mbox{ for } x = i, i + 1, …, j \tag{fmassuniformint}
 * @f]
 * and 0 elsewhere.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_discrete
 */
public class UniformIntGen extends RandomVariateGenInt {
   protected int left;     // the left limit of the interval
   protected int right;    // the right limit of the interval

   /**
    * Creates a uniform random variate generator over the integers in the
    * closed interval @f$[i, j]@f$, using stream `s`.
    */
   public UniformIntGen (RandomStream s, int i, int j) {
      super (s, new UniformIntDist (i, j));
      setParams (i, j);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public UniformIntGen (RandomStream s, UniformIntDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getI(), dist.getJ());
   }

   /**
    * Generates a new *uniform* random variate over the interval
    * @f$[i,j]@f$, using stream `s`, by inversion.
    */
   public static int nextInt (RandomStream s, int i, int j) {
      return UniformIntDist.inverseF (i, j, s.nextDouble());
   }

   /**
    * Returns the parameter @f$i@f$.
    */
   public int getI() {
      return left;
   }

   /**
    * Returns the parameter @f$j@f$.
    */
   public int getJ() {
      return right;
   }
   protected  void setParams (int i, int j) {
      if (j < i)
        throw new IllegalArgumentException ("j < i");
      left = i;
      right = j;
   }
}