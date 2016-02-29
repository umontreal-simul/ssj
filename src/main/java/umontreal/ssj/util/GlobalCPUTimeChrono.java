/*
 * Class:        GlobalCPUTimeChrono
 * Description:  Compute the global CPU time used by the Java Virtual Machine
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

/**
 * Extends the  @ref AbstractChrono class to compute the global CPU time used
 * by the Java Virtual Machine. This includes CPU time taken by any thread,
 * including the garbage collector, class loader, etc.
 *
 * Part of this class is implemented in the C language and the implementation
 * is unfortunately operating system-dependent. The C functions for the
 * current class have been compiled on a 32-bit machine running Linux. For a
 * *platform-independent* CPU timer (valid only with Java–1.5 or later), one
 * should use the class  @ref ThreadCPUTimeChrono which is programmed
 * directly in Java.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class GlobalCPUTimeChrono extends AbstractChrono {
   
   private static boolean ssjUtilLoaded = false;
   private static native void Heure (long[] tab);

   protected void getTime (long[] tab) {
      if (!ssjUtilLoaded) {
         // We cannot load the library in a static block
         // since the subclass ChronoSingleThread does not
         // need the native method.
         try {
            NativeUtils.loadLibraryFromJar("/jni/" + System.mapLibraryName("ssjutil"));
         }
         catch (java.io.IOException e) {
            throw new UnsatisfiedLinkError(e.getMessage());
         }

         ssjUtilLoaded = true;
      }
      
      Heure (tab);
   }

   /**
    * Constructs a `Chrono` object and initializes it to zero.
    */
   public GlobalCPUTimeChrono() {
      init();
   }

}
