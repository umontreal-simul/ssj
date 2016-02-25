/*
 * Class:        RandomStream
 * Description:  basic structures to handle multiple streams of uniform
                 (pseudo)-random numbers and tools to move around within
                 and across these streams
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Pierre L'Ecuyer
 * @since        2000
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
package umontreal.ssj.rng;

/**
 * This interface defines the basic structures to handle multiple streams of
 * uniform (pseudo)random numbers and convenient tools to move around within
 * and across these streams. The actual random number generators (RNGs) are
 * provided in classes that implement this `RandomStream` interface. Each
 * stream of random numbers is an object of the class that implements this
 * interface, and can be viewed as a virtual random number generator.
 *
 * For each type of base RNG (i.e., each implementation of the `RandomStream`
 * interface), the full period of the generator is cut into adjacent
 * *streams* (or segments) of length @f$Z@f$, and each of these streams is
 * partitioned into @f$V@f$ *substreams* of length @f$W@f$, where @f$Z =
 * VW@f$. The values of @f$V@f$ and @f$W@f$ depend on the specific RNG, but
 * are usually larger than @f$2^{50}@f$. Thus, the distance @f$Z@f$ between
 * the starting points of two successive streams provided by an RNG usually
 * exceeds @f$2^{100}@f$. The initial seed of the RNG is the starting point
 * of the first stream. It has a default value for each type of RNG, but this
 * initial value can be changed by calling `setPackageSeed` for the
 * corresponding class. Each time a new `RandomStream` is created, its
 * starting point (initial seed) is computed automatically, @f$Z@f$ steps
 * ahead of the starting point of the previously created stream of the same
 * type, and its current state is set equal to this starting point.
 *
 * For each stream, one can advance by one step and generate one value, or go
 * ahead to the beginning of the next substream within this stream, or go
 * back to the beginning of the current substream, or to the beginning of the
 * stream, or jump ahead or back by an arbitrary number of steps. Denote by
 * @f$C_g@f$ the current state of a stream @f$g@f$, @f$I_g@f$ its initial
 * state, @f$B_g@f$ the state at the beginning of the current substream, and
 * @f$N_g@f$ the state at the beginning of the next substream. The following
 * diagram shows an example of a stream whose state is at the 6th value of
 * the third substream, i.e., @f$2W+5@f$ steps ahead of its initial state
 * @f$I_g@f$ and 5 steps ahead of its state @f$B_g@f$. The form of the state
 * of a stream depends on its type. For example, the state of a stream of
 * class  @ref MRG32k3a is a vector of six 32-bit integers represented
 * internally as floating-point numbers (in <tt>double</tt>).
 *
 * @image html rng_randomstream_01.svg
 * <!--
 * LaTeX code used to generate the picture:
 *
 *  \def\tick #1{\vrule height 0pt depth #1pt} \def\enskip {\hskip .5em\relax
 * } \def\ld {\hbox to 0.24 in{\vtop {\kern 3.0pt\hbox {\dotfill }}}} \def\ts
 * {\enskip \tick 4} \def\suba {\hbox to 1.2in {\vtop {\hbox to 1.2in{\hbox
 * to .66in{\hrulefill }\hbox to .30in{\dotfill }\hrulefill } \hbox to
 * 1.2in{\tick 9\ts \ts \ts \ts \ts \ts \ts \hfill \ts \enskip }}}} \def\subb
 * {\hbox to 1.2in {\vtop {\hbox to 1.2in{\hbox to .84in {\hrulefill } \hbox
 * to .24in {\dotfill }\hbox to .12in {\hrulefill }} \hbox to 1.2in{\tick
 * 9\ts \ts \ts \ts \ts \ts \ts \ts \hfill \ts \enskip }}}}
 *
 * \[  \vbox{\offinterlineskip \hbox to 5.5in{\hskip173.448pt\hbox to
 * 1.2in{\hfil \hskip11.77pt$C_ g$\hfil }\hfill } \vskip0.1pt\hbox to
 * 5.5in{\hskip173.448pt\hbox to 1.2in{\hfil \hskip11.77pt$\Downarrow $\hfil
 * } \hfill } \vtop {\offinterlineskip \hskip0.0pt\hbox to 5.5 true in {\hbox
 * to 1.2in {\vtop {\hbox to 1.2in{\hbox to .66in{\hrulefill }\hbox to
 * .30in{\dotfill }\hrulefill } \hbox to 1.2in{\vrule height 0pt depth
 * 9pt\hskip5.5pt\relax \vrule height 0pt depth 4pt\hskip5.5pt\relax \vrule
 * height 0pt depth 4pt\hskip5.5pt\relax \vrule height 0pt depth
 * 4pt\hskip5.5pt\relax \vrule height 0pt depth 4pt\hskip5.5pt\relax \vrule
 * height 0pt depth 4pt\hskip5.5pt\relax \vrule height 0pt depth
 * 4pt\hskip5.5pt\relax \vrule height 0pt depth 4pt\hfill \hskip5.5pt\relax
 * \vrule height 0pt depth 4pt\hskip5.5pt\relax }}}\hbox to 1.2in {\vtop
 * {\hbox to 1.2in{\hbox to .66in{\hrulefill }\hbox to .30in{\dotfill
 * }\hrulefill } \hbox to 1.2in{\vrule height 0pt depth 9pt\hskip5.5pt\relax
 * \vrule height 0pt depth 4pt\hskip5.5pt\relax \vrule height 0pt depth
 * 4pt\hskip5.5pt\relax \vrule height 0pt depth 4pt\hskip5.5pt\relax \vrule
 * height 0pt depth 4pt\hskip5.5pt\relax \vrule height 0pt depth
 * 4pt\hskip5.5pt\relax \vrule height 0pt depth 4pt\hskip5.5pt\relax \vrule
 * height 0pt depth 4pt\hfill \hskip5.5pt\relax \vrule height 0pt depth
 * 4pt\hskip5.5pt\relax }}}\hbox to 1.2in {\vtop {\hbox to 1.2in{\hbox to
 * .84in {\hrulefill } \hbox to .24in {\dotfill }\hbox to .12in {\hrulefill
 * }} \hbox to 1.2in{\vrule height 0pt depth 9pt\hskip5.5pt\relax \vrule
 * height 0pt depth 4pt\hskip5.5pt\relax \vrule height 0pt depth
 * 4pt\hskip5.5pt\relax \vrule height 0pt depth 4pt\hskip5.5pt\relax \vrule
 * height 0pt depth 4pt\hskip5.5pt\relax \vrule height 0pt depth
 * 4pt\hskip5.5pt\relax \vrule height 0pt depth 4pt\hskip5.5pt\relax \vrule
 * height 0pt depth 4pt\hskip5.5pt\relax \vrule height 0pt depth 4pt\hfill
 * \hskip5.5pt\relax \vrule height 0pt depth 4pt\hskip5.5pt\relax }}}\hbox to
 * 1.2in {\vtop {\hbox to 1.2in{\hbox to .66in{\hrulefill }\hbox to
 * .30in{\dotfill }\hrulefill } \hbox to 1.2in{\vrule height 0pt depth
 * 9pt\hskip5.5pt\relax \vrule height 0pt depth 4pt\hskip5.5pt\relax \vrule
 * height 0pt depth 4pt\hskip5.5pt\relax \vrule height 0pt depth
 * 4pt\hskip5.5pt\relax \vrule height 0pt depth 4pt\hskip5.5pt\relax \vrule
 * height 0pt depth 4pt\hskip5.5pt\relax \vrule height 0pt depth
 * 4pt\hskip5.5pt\relax \vrule height 0pt depth 4pt\hfill \hskip5.5pt\relax
 * \vrule height 0pt depth 4pt\hskip5.5pt\relax }}}\vrule height 0pt depth
 * 9pt \hbox to .12in {\hrulefill }\hbox to .24in{\dotfill
 * }}}\vskip0.1pt\hskip-46.9755pt\hbox to 1.2 in{\hfil $I_ g$\hfil
 * }\hskip90.3375pt\hbox to 1.2 in{\hfil $B_ g$\hfil }\hbox to 1.2 in{\hfil
 * $N_ g$\hfil }\hfill \vskip0.1pt}  \]
 *
 * -->
 *
 * The methods for manipulating the streams and generating random numbers are
 * implemented differently for each type of RNG. The methods whose formal
 * parameter types do not depend on the RNG type are specified in the
 * interface `RandomStream`. The others (e.g., for setting the seeds) are
 * given only in the classes that implement the specific RNG types.
 *
 * See @cite sLAW00a, @cite rLEC91a, @cite rLEC02a&thinsp; for examples of
 * situations where the multiple streams offered here are useful.
 *
 * Methods for generating random variates from non-uniform distributions are
 * provided in the  @ref umontreal.ssj.randvar package.
 *
 * <div class="SSJ-bigskip"></div>
 */
public interface RandomStream {

   /**
    * Reinitializes the stream to its initial state @f$I_g@f$: @f$C_g@f$
    * and @f$B_g@f$ are set to @f$I_g@f$.
    */
   public void resetStartStream();

   /**
    * Reinitializes the stream to the beginning of its current substream:
    * @f$C_g@f$ is set to @f$B_g@f$.
    */
   public void resetStartSubstream();

   /**
    * Reinitializes the stream to the beginning of its next substream:
    * @f$N_g@f$ is computed, and @f$C_g@f$ and @f$B_g@f$ are set to
    * @f$N_g@f$.
    */
   public void resetNextSubstream();

   /**
    * Returns a string containing the current state of this stream.
    *  @return the state of the generator formated as a string
    */
   public String toString();

   /**
    * Returns a (pseudo)random number from the uniform distribution over
    * the interval @f$(0,1)@f$, using this stream, after advancing its
    * state by one step. The generators programmed in SSJ never return the
    * values 0 or 1.
    *  @return the next generated uniform
    */
   public double nextDouble();

   /**
    * Generates `n` (pseudo)random numbers from the uniform distribution
    * and stores them into the array `u` starting at index `start`.
    *  @param u            array that will contain the generated uniforms
    *  @param start        starting index, in the array `u`, to write
    *                      uniforms from
    *  @param n            number of uniforms to generate
    */
   public void nextArrayOfDouble (double[] u, int start, int n);

   /**
    * Returns a (pseudo)random number from the discrete uniform
    * distribution over the integers @f$\{i,i+1,…,j\}@f$, using this
    * stream. (Calls `nextDouble` once.)
    *  @param i            smallest integer that can be generated
    *  @param j            greatest integer that can be generated
    *  @return the generated integer
    */
   public int nextInt (int i, int j);

   /**
    * Generates `n` (pseudo)random numbers from the discrete uniform
    * distribution over the integers @f$\{i,i+1,…,j\}@f$, using this
    * stream and stores the result in the array `u` starting at index
    * `start`. (Calls `nextInt` `n` times.)
    *  @param i            smallest integer that can be generated
    *  @param j            greatest integer that can be generated
    *  @param u            array that will contain the generated values
    *  @param start        starting index, in the array `u`, to write
    *                      integers from
    *  @param n            number of values being generated
    */
   public void nextArrayOfInt (int i, int j, int[] u, int start, int n);
 
}