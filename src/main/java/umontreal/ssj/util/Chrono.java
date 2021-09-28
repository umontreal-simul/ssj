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
 * class and computes the CPU time for the current thread only.  It is equivalent to
 * @ref ChronoSingleThread.
 * The class @ref umontreal.ssj.util.ChronoWall provides another option.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class Chrono extends AbstractChrono {
   private ChronoSingleThread chrono = new ChronoSingleThread();

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
    * Creates a `Chrono` instance adapted for a program using a single thread.
    * It is equivalent to @ref ChronoSingleThread.
    * This class should not be used to create a timer for a
    * multi-threaded program, because the obtained CPU times will differ
    * depending on the used Java version.
    *  @return the constructed timer.
    */
   public static Chrono createForSingleThread () {
         return new Chrono();
   }

}