/*
 * Class:        Chrono
 * Description:  computes the CPU time for the current thread only
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Éric Buist
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
 * The  @ref Chrono class extends the  @ref umontreal.ssj.util.AbstractChrono
 * class and computes the CPU time for the current thread only. This is the
 * simplest way to use chronos. Classes `AbstractChrono`,
 * @ref umontreal.ssj.util.SystemTimeChrono,
 * @ref umontreal.ssj.util.GlobalCPUTimeChrono and
 * @ref umontreal.ssj.util.ThreadCPUTimeChrono provide different chronos
 * implementations. See these classes to learn more about SSJ chronos, if
 * problems appear with class `Chrono`.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class Chrono extends AbstractChrono {
   private ThreadCPUTimeChrono chrono = new ThreadCPUTimeChrono();

   protected void getTime (long[] tab) {
         chrono.getTime(tab);
   }

   /**
    * Constructs a `Chrono` object and initializes it to zero.
    */
   public Chrono() {
      chrono.init();
      init();
   }

   /**
    * Creates a `Chrono` instance adapted for a program using a single
    * thread. Under Java 1.5, this method returns an instance of
    * @ref ChronoSingleThread which can measure CPU time for one thread.
    * Under Java versions prior to 1.5, this returns an instance of this
    * class. This method must not be used to create a timer for a
    * multi-threaded program, because the obtained CPU times will differ
    * depending on the used Java version.
    *  @return the constructed timer.
    */
   public static Chrono createForSingleThread () {
         return new Chrono();
   }

}