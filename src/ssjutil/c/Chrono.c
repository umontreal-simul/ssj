/* Interface C --> Java pour GlobalCPUTimeChrono.java */

#include <stdio.h>
#include <stdlib.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include <time.h>
#ifdef _linux
#include <unistd.h>
#include <sys/times.h>
#endif
#include "umontreal_ssj_util_GlobalCPUTimeChrono.h"

#ifdef _WIN32
static HANDLE currentProcess;

JNIEXPORT jint JNICALL
JNI_OnLoad (JavaVM* vm, void* reserved) {
  currentProcess = GetCurrentProcess();
  return JNI_VERSION_1_2;
}

/*
 * A helper function for converting FILETIME to a LONGLONG [safe from memory
 * alignment point of view].
 */
static ULONGLONG
fileTimeToInt64 (const FILETIME * time)
{
    ULARGE_INTEGER _time;

    _time.LowPart = time->dwLowDateTime;
    _time.HighPart = time->dwHighDateTime;

    return _time.QuadPart;
}
#endif

JNIEXPORT void JNICALL 
Java_umontreal_ssj_util_GlobalCPUTimeChrono_Heure (JNIEnv *env, jclass class, jlongArray array){

#ifdef _linux
  struct tms us;
#endif
#ifdef _WIN32
  FILETIME creationTime, exitTime, kernelTime, userTime;
#endif
  jlong *jarray = (*env)->GetLongArrayElements(env, array, 0);

#if defined(_linux)
  long TIC = sysconf(_SC_CLK_TCK);

  times(&us);

  jarray[1] = us.tms_utime + us.tms_stime;
  jarray[0] = jarray[1] / TIC;
  jarray[1] = (jarray[1] % TIC) * 1000000 / TIC;
#elif defined(_WIN32)
  /* Strongly inspired from
   * http://www.javaworld.com/javaworld/javaqa/2002-11/01-qa-1108-cpu.html */
  GetProcessTimes (currentProcess, &creationTime, &exitTime,
		   &kernelTime, &userTime);
  ULONGLONG rawTime = (ULONGLONG)(fileTimeToInt64 (&kernelTime) +
                                  fileTimeToInt64 (&userTime));
  /* We have to divide by 10000 to get milliseconds out of
   * the computed time. So we divide by 10000*1000 to get seconds. */
  jarray[0] = (unsigned long)(rawTime / 10000000);
  /* One raw time unit corresponds to 10 microseconds.
   */
  jarray[1] = (unsigned long)((rawTime % 10000000) / 10);
#else
  /* This one is bad but portable across Unixes.
   * The clock function wraps after approximately 72 minutes. */
  jarray[1] = clock();
  jarray[0] = jarray[1] / CLOCKS_PER_SEC;
  jarray[1] = (jarray[1] % CLOCKS_PER_SEC) * 1000000 / CLOCKS_PER_SEC;
#endif

  (*env)->ReleaseLongArrayElements(env, array, jarray, 0);
}
