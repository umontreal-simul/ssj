/*
 * Class:        SystemTimeChrono
 * Description:  
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
 * Extends the  @ref AbstractChrono class to compute the total system time
 * using Java’s builtin `System.nanoTime`. The system can be used as a rough
 * approximation of the CPU time taken by a program if no other tasks are
 * executed on the host while the program is running.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class SystemTimeChrono extends AbstractChrono {

   protected void getTime (long[] tab) {
      long rawTime = System.nanoTime();
      final long DIV = 1000000000L;
      long seconds = rawTime/DIV;
      long micros = (rawTime % DIV)/1000L;
      tab[0] = seconds;
      tab[1] = micros;
   }

   /**
    * Constructs a new chrono object and initializes it to zero.
    */
   public SystemTimeChrono() {
      super();
      init();
   }

}