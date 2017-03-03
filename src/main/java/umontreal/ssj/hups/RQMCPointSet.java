/*
 * Class:        RQMCPointSet
 * Description:  randomized quasi-Monte Carlo simulations
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
 * This class is used for *randomized quasi-Monte Carlo* (RQMC) simulations
 * @cite vLEC00b, @cite vLEC02a, @cite vOWE97a, @cite vOWE97b&thinsp;. The
 * idea is to randomize a point set so that:
 *
 * <ul><li>
 * it retains its high uniformity when taken as a set and
 * </li>
 * <li>
 * each individual point is a random vector with the uniform distribution
 * over @f$(0, 1)^s@f$.
 * </li>
 * </ul>
 *
 *  A RQMC point set is one that satisfies these two conditions. One simple
 * randomization that satisfies these conditions for an arbirary point set
 * @f$P_n@f$ is a random shift modulo 1 @cite vCRA76a, @cite vLEC00b,
 * @cite vSLO94a&thinsp;: Generate a single point @f$\mathbf{U}@f$ uniformly
 * over @f$(0, 1)^s@f$ and add it to each point of @f$P_n@f$, modulo 1,
 * coordinate-wise. Another one is a random digital shift in base @f$b@f$
 * @cite vLEC99a, @cite vLEC02a, @cite mMAT99a&thinsp;: generate again
 * @f$\mathbf{U}@f$ uniformly over @f$(0, 1)^s@f$, expand each of its
 * coordinates in base @f$b@f$, and add the digits, modulo @f$b@f$, to the
 * corresponding digits of each point of @f$P_n@f$. <div
 * class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class RQMCPointSet {
   private PointSet set;
   private PointSetRandomization rand;

   /**
    * Constructor with the point set `set` and the randomization `rand`.
    *  @param set          the point set
    *  @param rand         the randomization
    */
   public RQMCPointSet (PointSet set, PointSetRandomization rand) {
      this.rand = rand;
      this.set = set;
   }

   /**
    * Randomizes the point set. The randomization and the point set are
    * those of this object.
    */
   public void randomize() {
       rand.randomize(set);
   }

   /**
    * Returns a new point set iterator for the point set associated to
    * this object.
    *  @return point set iterator for the point set
    */
   public PointSetIterator iterator() {
      return set.iterator();
   }

   /**
    * Returns the point set associated to this object.
    *  @return the point set associated to this object
    */
   public PointSet getPointSet() {
      return set;
   }

   /**
    * Returns the randomization associated to this object.
    *  @return the randomization associated to this object
    */
   public PointSetRandomization getRandomization() {
      return rand;
   }

   /**
    * Returns the number of points in the associated point set.
    *  @return the number of points in point set
    */
   public int getNumPoints() {
      return set.getNumPoints();
   }

}