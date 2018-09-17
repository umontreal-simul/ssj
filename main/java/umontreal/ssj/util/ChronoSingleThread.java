/*
 * Class:        ChronoSingleThread
 * Description:  deprecated
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
package umontreal.ssj.util;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;

/**
 * **This class is deprecated** but kept for compatibility with older
 * versions of SSJ.  @ref Chrono should be used instead of
 * @ref ChronoSingleThread. The  @ref ChronoSingleThread class extends the
 * @ref AbstractChrono class and computes the CPU time for the current thread
 * only. This is the simplest way to use chronos. Classes
 * @ref AbstractChrono,  @ref SystemTimeChrono,  @ref GlobalCPUTimeChrono and
 * @ref ThreadCPUTimeChrono provide different chronos implementations (see
 * these classes to learn more about SSJ chronos).
 *
 * <div class="SSJ-bigskip"></div>
 */
@Deprecated
public class ChronoSingleThread extends AbstractChrono {

   private ThreadCPUTimeChrono chrono = new ThreadCPUTimeChrono();

   protected void getTime (long[] tab) {
         chrono.getTime(tab);
   }

   /**
    * Constructs a `ChronoSingleThread` object and initializes it to zero.
    */
   public ChronoSingleThread() {
      chrono.init();
      init();
   }

}