#include <string.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include <unuran.h>
#include "umontreal_ssj_randvar_RandUnuran.h"

/* Native interface for UNURAN with Java
 *
 * UNURAN has been set up to use a custom generator represented
 * by a data structure containing a function pointer and a parameter
 * pointer. This setup is different from the two supported configuration.
 * When using UNUR_URNG_POINTER, only a function pointer is used.
 * When using UNUR_URNG_PRNG, a PRNG generator is used and an extra
 * library (prng) is needed to compile UNURAN and RandUnuran.
 * When a uniform must be generated, the function pointed to
 * by the pointer getrand is called using the parameter pointer params.
 * The function has to return a double, which will be the generated
 * uniform. This structure has been defined in unuran.h so UNUR_URNG
 * is a typedef of it. This structure, called UNUR_URNG_FUNCTION_POINTER
 * is general and can solve other cases than the JNI.
 *
 * For our special case, we would like to be able to call a Java
 * method of a RandomStream Java object to get the uniform numbers.
 * The UNUR_URNG will act as a wrapper around the Java method.
 * So we need a function to issue a JNI call. However, such a call
 * requires some information: the JNI environment, the reference
 * to the object and the ID of the method.
 * We are using a struct GenParams as the parameter
 * block. This block contains all the needed information to allow
 * the function unif, which getrand will be assigned to, to call the JNI.
 * We have one GenParams structure for each UNUR_URNG generator
 * used and two UNUR_URNG's are used for each RandUnuran Java object,
 * one for the main stream and one for the auxiliary stream.
 *
 * This uniform generator structure allows us to change the function
 * dynamically. We can therefore perform optimization by executing
 * as many RandomStream method calls as possible on the Java side
 * instead of the native side. When generating a single variate,
 * we cache one uniform and use the cons function to get it.
 * When generating arrays of variate, the function arr is used.
 * It gets elements from an array and can refill it with
 * new random number if necessary.
 *
 * Additionally, we need another structure called NativeParams.
 * One such structure is assigned to each RandUnuran Java object.
 * It contains native-side specific parameters such as pointers
 * to the UNUR_GEN structure and the two UNUR_URNG's.
 *
 *
 * The life of a RandUnuran object
 *
 * When constructed, a RandUnuran object will always call
 * the native method init. The Java call init will be
 * translated to a C call to Java_umontreal_ssj_randvar_RandUnuran_init with
 * the proper arguments. This method will allocate the
 * NativeParams structure and assign its pointer to an integer
 * instance field named nativeParams. This will allow us to attach the
 * NativeParams structure to the Java object and get it back
 * from calls to calls. Then UNURAN is called in order
 * to create a generator object which will be stored in
 * the NativeParams C structure. Two UNUR_URNG structures are then
 * allocated to store information about main and auxiliary RandomStream
 * objects. The newly-allocated uniform random number generators will
 * be attached to the UNURAN generator by UNURAN function calls.
 *
 * When the user asks for a random number, a Java method uses
 * the mainStream RandomStream object to get a first uniform from the
 * main stream. This uniform u is passed to the native side by a call
 * to getRandDisc or getRandCont and stored
 * inside the GenParams structure stored in the UNUR_URNG structure
 * for the main stream. The getrand pointer of the main generator is also
 * changed to cons. Environment and object reference are also
 * stored into the GenParams structure. UNURAN can then be called.
 * During the computation process performed by UNURAN, at least
 * one uniform random number will be needed. This will call the cons
 * function with a void pointer. This pointer will be casted to
 * the GenParams pointer assigned to the main stream. The cons
 * function will set the function pointer of the UNUR_URNG structure
 * assigned to the main stream to unif. If UNURAN requests
 * extra numbers from the main stream during the same variate generation
 * process, unif is called.
 * It will issue a JNI call that will return to Java side in order
 * to use the nextDouble method of the RandomStream object.
 *
 * When no references point to the RandUnuran object, the Java garbage
 * collector will free it. Before releasing the memory, the finalize
 * method is called, which will call the native close method. This method
 * frees all memory allocated by the native part of the interface,
 * including the UNURAN generator object.
 */

/* cached data for efficiency
 * These ids are needed to issue JNI calls for getting fields
 * and calling methods.
 */
/* fields for RandUnuran objects */
static jfieldID fidNativeParams = 0;
static jfieldID fidMainStream = 0;
static jfieldID fidAuxStream = 0;

/* methods for RandomStream objects */
static jmethodID midNextDouble = 0;
static jmethodID midNextArrayOfDouble;

struct UrngWithParams {
   UNUR_URNG* urng;
   struct GenParams* params;
};

/* We must keep a default uniform generator for the
 * generators' setup phases (call to unur_str2gen).
 * This generator does not creates any Java RandomStream objects. */
static struct UrngWithParams default_urng = { NULL, NULL };

struct GenParams {
   /* Parameters block attached to one uniform random number generator
    * This will be passed as a void pointer to the function pointed to by
    * the urng. */

   /* JNI information */
   JNIEnv* env; /* the JNI environment */
   jobject rsObj; /* the RandomStream object reference */

   /* caching of one uniform */
   double u; /* first random number obtained on the Java side */
   UNUR_URNG* urng; /* needed to set back to unif after cons is called */

   /* array of uniforms */
   jdoubleArray junifArray; /* Reference to the Java array */
   jdouble* unifArray; /* array of uniforms */
   int n; /* number of generated values */
   int nextIndex; /* index of the next uniform */

   double (*unif)(struct GenParams* p);
};

static double unif_wrapper (void* p) {
   /* Function used to generate uniforms using a Java RandomStream. */

   struct GenParams* params = (struct GenParams*)p;
   return params->unif (params);
}

static double unif (struct GenParams* params) {
   /* Function used to generate uniforms using a Java RandomStream. */

   return (*(params->env))->CallDoubleMethod (params->env, params->rsObj,
					      midNextDouble);
}

static double cons (struct GenParams* params) {
   /* We can prevent some native Java calls by precomputing one
    * uniform before entering the native side. If UNURAN asks for
    * additional random numbers, this will call the unif function.
    * This function avoids an if inside the unif function. */

   params->unif = unif;
   return params->u;
}

static double arr (void* p) {
   /* Gets a uniform from a previously generated
    * array of uniforms. If no more uniforms are available,
    * the array is refilled. */

   struct GenParams* params = (struct GenParams*)p;
   if (params->nextIndex >= params->n) {
      /* a refill is needed */
      (*(params->env))->ReleaseDoubleArrayElements (params->env,
						    params->junifArray,
						    params->unifArray, 0);
      /* equivalent to
       * mainStream.nextArrayOfDouble (unifArray, 0, n)
       * or
       * auxStream.nextArrayOfDouble (...)
       */
      (*(params->env))->CallVoidMethod (params->env, params->rsObj,
					midNextArrayOfDouble,
					params->junifArray, 0, params->n);
      params->unifArray =
	 (*(params->env))->GetDoubleArrayElements (params->env,
                                                   params->junifArray, 0);
      params->nextIndex = 0;
   }

   /* gets a uniform */
   return params->unifArray[params->nextIndex++];
}

/* Per-object parameters used in C */
struct NativeParams {
   /* We store it in a single structure in order to limit
    * the necessary calls to the JNI. Additionnally, some
    * information cannot be directly stored in Java fields, like
    * pointers.
    * The objMainStream and objAuxStream fields are replicated inside
    * the parameter block of the urng and urng_aux structures but
    * it is more simple to free the references inside the
    * NativeParams structure. */

   UNUR_GEN* unurgen; /* non-uniform generator */
   jobject objMainStream; /* uniform generator RandomStream */
   jobject objAuxStream; /* uniform auxiliary generator RandomStream */
   struct UrngWithParams urng; /* uniform generator structure */
   struct UrngWithParams urng_aux; /* uniform auxiliary generator structure */

   /* The dimension of the distribution is used for array size
    * checking in case of multivariate distributions. */
   int dim;
};

static void urng_jni_destroy (struct UrngWithParams* gen) {
   /* This function frees the custom JNI generator. */

   if (!gen)
      return;
   if (gen->params) {
      free (gen->params);
      gen->params = NULL;
   }
   if (gen->urng) {
      unur_urng_free (gen->urng);
      gen->urng = NULL;
   }
}

static void urng_jni_create (struct UrngWithParams* gen) {
   /* Creates a new JNI generator structure. */

   gen->params = malloc (sizeof (struct GenParams));
   if (!gen->params) {
      urng_jni_destroy (gen);
      return;
   }
   gen->params->unif = unif;
   gen->urng = unur_urng_new (unif_wrapper, gen->params);
   if (!gen->urng) {
      urng_jni_destroy (gen);
      return;
   }
}

static int initClass (JNIEnv* env) {
   /* we perform some class initialization tasks.
    * If this function returns 0, the initialization failed. */

   jclass cls;

   /* Create and set the default uniform random number generator. */
   if (!default_urng.urng) {
      urng_jni_create (&default_urng);
      if (!default_urng.urng) {
	 cls = (*env)->FindClass (env, "umontreal/ssj/randvar/UnuranException");
	 if (!cls)
	    return 0;
	 (*env)->ThrowNew (env, cls, "cannot create default uniform rng");
	 return 0;
      }
      unur_set_default_urng (default_urng.urng);
      /* After this setup, UNURAN won't work properly
       * until we set a valid parameter block for the unif
       * function. This will be done in the init method only. */
   }

   cls = (*env)->FindClass (env, "umontreal/ssj/randvar/RandUnuran");
   if (!cls)
      return 0;
   if (!fidNativeParams &&
       !(fidNativeParams = (*env)->GetFieldID (env, cls,
					       "nativeParams", "J")))
      return 0;
   if (!fidMainStream &&
       !(fidMainStream = (*env)->GetFieldID
	 (env, cls, "mainStream", "Lumontreal/ssj/rng/RandomStream;")))
      return 0;
   if (!fidAuxStream &&
       !(fidAuxStream = (*env)->GetFieldID
	 (env, cls, "auxStream", "Lumontreal/ssj/rng/RandomStream;")))
      return 0;

   cls = (*env)->FindClass (env, "umontreal/ssj/rng/RandomStream");
   if (!cls)
      return 0;
   if (!midNextDouble &&
       !(midNextDouble = (*env)->GetMethodID (env, cls, "nextDouble", "()D")))
      return 0;
   if (!midNextArrayOfDouble &&
       !(midNextArrayOfDouble = (*env)->GetMethodID (env, cls,
						     "nextArrayOfDouble",
						  "([DII)V")))
      return 0;

   return 1;
}

JNIEXPORT void JNICALL Java_umontreal_ssj_randvar_RandUnuran_init
  (JNIEnv *env, jobject obj, jstring genStr)
{
   jclass expcls;
   jobject ob;
   const char* genStrUTF;
   struct NativeParams* np;
   struct GenParams *gp;

   if (!initClass (env))
      return;

   /* Allocate the object parameter block. */
   np = (struct NativeParams*)malloc (sizeof (struct NativeParams));
   if (!np) {
      expcls = (*env)->FindClass (env, "umontreal/ssj/randvar/UnuranException");
      if (!expcls)
	 return;
      (*env)->ThrowNew (env, expcls, "cannot create UNURAN generator");
      return;
   }
   /* This will allow us to call close upon errors */
   memset (np, 0, sizeof (struct NativeParams));
   /*
   np->unurgen = 0;
   np->objMainStream = 0;
   np->objAuxStream = 0;
   np->urng = 0;
   np->urng_aux = 0; */

   /* save the parameters with the Java object */
   (*env)->SetLongField (env, obj, fidNativeParams, (jlong)np);

   /* Set the parameter block for the default generator.
    * This will allow setup-time random number generation
    * by the Java RandomStream. The setup has to be made
    * all the times because the env and obj vary from calls to calls. */
   gp = default_urng.params;
   gp->env = env;
   gp->rsObj = (*env)->GetObjectField (env, obj, fidMainStream);

   /* create the UNURAN generator */
   genStrUTF = (*env)->GetStringUTFChars (env, genStr, 0);
   if (!genStrUTF) {
      Java_umontreal_ssj_randvar_RandUnuran_close (env, obj);
      return;
   }
   np->unurgen = unur_str2gen (genStrUTF);
   (*env)->ReleaseStringUTFChars (env, genStr, genStrUTF);
   if (!np->unurgen) {
      /* The creation of the generator has failed.
       * Throws an exception that will report the error,
       * including the UNURAN error message. */
      const char* unurerr = unur_get_strerror (unur_errno);
      const char* errmsg = "cannot create UNURAN generator: ";
      char* errstr;
      expcls = (*env)->FindClass (env, "umontreal/ssj/randvar/UnuranException");
      if (!expcls) {
	 Java_umontreal_ssj_randvar_RandUnuran_close (env, obj);
	 return;
      }
      errstr = (char*)malloc (strlen(errmsg) + strlen(unurerr) + 1);
      if (!errstr) {
	 Java_umontreal_ssj_randvar_RandUnuran_close (env, obj);
	 (*env)->ThrowNew (env, expcls, errmsg);
	 return;
      }
      strcpy (errstr, errmsg);
      strcat (errstr, unurerr);
      Java_umontreal_ssj_randvar_RandUnuran_close (env, obj);
      (*env)->ThrowNew (env, expcls, errstr);
      free (errstr);
      return;
   }
   np->dim = unur_get_dimension (np->unurgen);

   /* now we must setup the uniform generator */
   urng_jni_create (&np->urng);
   if (!np->urng.urng) {
      Java_umontreal_ssj_randvar_RandUnuran_close (env, obj);
      expcls = (*env)->FindClass (env, "umontreal/ssj/randvar/UnuranException");
      if (!expcls)
	 return;
      (*env)->ThrowNew (env, expcls, "cannot allocate uniform "
                                     "random number generator");
      return;
   }
   gp = np->urng.params;
   gp->urng = np->urng.urng;
   unur_chg_urng (np->unurgen, np->urng.urng);

   /* now we must setup the auxiliary uniform generator */
   urng_jni_create (&np->urng_aux);
   if (!np->urng_aux.urng) {
      Java_umontreal_ssj_randvar_RandUnuran_close (env, obj);
      expcls = (*env)->FindClass (env, "umontreal/ssj/randvar/UnuranException");
      if (!expcls)
	 return;
      (*env)->ThrowNew (env, expcls, "cannot allocate uniform auxiliary "
                                     "random number generator.");
      return;
   }
   unur_chg_urng_aux (np->unurgen, np->urng_aux.urng);

   /* Get reference to the RandomStream mainStream.
    * The references to the RandomStream objects have
    * to be global because they are kept inside the NativeParams
    * structure which is passed to other instance methods. */
   ob = (*env)->GetObjectField (env, obj, fidMainStream);
   if (!ob) {
      Java_umontreal_ssj_randvar_RandUnuran_close (env, obj);
      return;
   }
   np->objMainStream = (*env)->NewGlobalRef (env, ob);
   if (!np->objMainStream) {
      Java_umontreal_ssj_randvar_RandUnuran_close (env, obj);
      return;
   }
   /* store a copy of the global reference inside the
    * random number generator parameter block. */
   gp = np->urng.params;
   gp->rsObj = np->objMainStream;

   /* Get reference to the RandomStream auxStream. */
   ob = (*env)->GetObjectField (env, obj, fidAuxStream);
   if (!ob) {
      Java_umontreal_ssj_randvar_RandUnuran_close (env, obj);
      return;
   }
   np->objAuxStream = (*env)->NewGlobalRef (env, ob);
   if (!np->objAuxStream) {
      Java_umontreal_ssj_randvar_RandUnuran_close (env, obj);
      return;
   }
   gp = np->urng_aux.params;
   gp->rsObj = np->objAuxStream;
}

JNIEXPORT void JNICALL Java_umontreal_ssj_randvar_RandUnuran_close
  (JNIEnv *env, jobject obj)
{
   struct NativeParams* np;
   jthrowable ex;

   /* Since this method is also called upon errors inside init,
    * we must ensure that no exception is pending before
    * calling JNI. */
   ex = (*env)->ExceptionOccurred (env);
   (*env)->ExceptionClear (env);

   np = (struct NativeParams*)((*env)->GetLongField
			       (env, obj, fidNativeParams));
   /* We check whether the pointer is null to prevent segmentation
    * fault that would result in a VM crash.
    * Without this checking, a VM crash would occur every time
    * the user calls this method manually because the finalizer
    * would call it and get an invalid pointer. */
   if (!np) {
      if (ex)
	 (*env)->Throw (env, ex);
      return;
   }

   /* We free memory previously allocated */
   urng_jni_destroy (&np->urng_aux);
   urng_jni_destroy (&np->urng);
   unur_free (np->unurgen);

   if (np->objAuxStream)
      (*env)->DeleteGlobalRef (env, np->objAuxStream);
   if (np->objMainStream)
      (*env)->DeleteGlobalRef (env, np->objMainStream);
   free (np);

   /* to avoid receiving an invalid pointer in a subsequent
    * call to close. */
   (*env)->SetLongField (env, obj, fidNativeParams, 0);

   /* If an exception was pending at the beginning of this method,
    * we can now throw it back. */
   if (ex)
      (*env)->Throw (env, ex);
}

JNIEXPORT jint JNICALL Java_umontreal_ssj_randvar_RandUnuran_getRandDisc
  (JNIEnv *env, jobject obj, jdouble u, jlong _np)
{
   struct NativeParams* np;
   struct GenParams *gp;

   np = (struct NativeParams*)_np;

   gp = np->urng.params;
   gp->env = env;

   /* We put the pregenerated uniform into the u field
    * and change the generation function. */
   gp->u = u;
   np->urng.params->unif = cons;

   gp = np->urng_aux.params;
   gp->env = env;
   np->urng_aux.params->unif = unif;

   return unur_sample_discr (np->unurgen);
}

JNIEXPORT jdouble JNICALL Java_umontreal_ssj_randvar_RandUnuran_getRandCont
  (JNIEnv *env, jobject obj, jdouble u, jlong _np)
{
   struct NativeParams* np;
   struct GenParams *gp;

   np = (struct NativeParams*)_np;

   gp = np->urng.params;

   gp->env = env;

   /* We put the pregenerated uniform into the u field
    * and change the generation function. */
   gp->u = u;
   np->urng.params->unif = cons;

   gp = np->urng_aux.params;
   gp->env = env;
   np->urng_aux.params->unif = unif;

   return unur_sample_cont (np->unurgen);
}

JNIEXPORT void JNICALL Java_umontreal_ssj_randvar_RandUnuran_getRandVec
  (JNIEnv *env, jobject obj, jdouble u, jlong _np, jdoubleArray jvec)
{
   struct NativeParams* np;
   struct GenParams *gp;
   jdouble* vec;

   np = (struct NativeParams*)_np;

   gp = np->urng.params;

   /* We put the pregenerated uniform into the u field
    * and change the generation function. */
   gp->u = u;
   np->urng.params->unif = cons;

   gp->env = env;

   gp = np->urng_aux.params;
   gp->env = env;
   np->urng_aux.params->unif = unif;

   /* size checking is required because a VM crash
    * would occur if the array is too short. */
   if ((*env)->GetArrayLength (env, jvec) < np->dim) {
      jclass expcls = (*env)->FindClass (env, "umontreal/ssj/randvar/UnuranException");
      if (!expcls)
	 return;
      (*env)->ThrowNew (env, expcls, "array too short");
      return;
   }
   vec = (*env)->GetDoubleArrayElements (env, jvec, 0);
   if (!vec)
      return;
   unur_sample_vec (np->unurgen, vec);
   (*env)->ReleaseDoubleArrayElements (env, jvec, vec, 0);
}

JNIEXPORT void JNICALL Java_umontreal_ssj_randvar_RandUnuran_getRandDiscArray
  (JNIEnv *env, jobject obj, jlong _np,
   jdoubleArray ju, jdoubleArray juaux, jintArray jv,
   jint start, jint n)
{
   struct NativeParams* np;
   struct GenParams *gp;
   jdouble* u;
   jdouble* uaux;
   jint* v;
   int i;

   np = (struct NativeParams*)_np;

   gp = np->urng.params;
   gp->env = env;
   gp->junifArray = ju;
   u = gp->unifArray = (*env)->GetDoubleArrayElements (env, ju, 0);
   gp->n = gp->nextIndex = n;    /* schedule an array refill */
   np->urng.params->unif = arr;

   gp = np->urng_aux.params;
   if (juaux) {
      gp->env = env;
      gp->junifArray = juaux;
      uaux = gp->unifArray = (*env)->GetDoubleArrayElements (env, juaux, 0);
      gp->n = gp->nextIndex = n;
      np->urng.params->unif = arr;
   }
   else {
      /* If the main and auxiliary streams are the same,
       * we use the same array for both, so the parameter
       * blocks must be equal. */
      np->urng_aux.params = np->urng.params;
      uaux = 0;
   }

   v = (*env)->GetIntArrayElements (env, jv, 0);
   for (i = 0; i < n; i++)
      v[start+i] = unur_sample_discr (np->unurgen);
   (*env)->ReleaseIntArrayElements (env, jv, v, 0);
   if (juaux) {
      (*env)->ReleaseDoubleArrayElements (env, juaux, uaux, 0);
   }
   else {
      np->urng_aux.params = gp;
   }
   (*env)->ReleaseDoubleArrayElements (env, ju, u, 0);
}

JNIEXPORT void JNICALL Java_umontreal_ssj_randvar_RandUnuran_getRandContArray
  (JNIEnv *env, jobject obj, jlong _np,
   jdoubleArray ju, jdoubleArray juaux, jdoubleArray jv,
   jint start, jint n)
{
   struct NativeParams* np;
   struct GenParams *gp;
   jdouble* u;
   jdouble* uaux;
   jdouble* v;
   int i;

   np = (struct NativeParams*)_np;

   gp = np->urng.params;
   gp->env = env;
   gp->junifArray = ju;
   u = gp->unifArray = (*env)->GetDoubleArrayElements (env, ju, 0);
   gp->n = gp->nextIndex = n;
   np->urng.params->unif = arr;

   gp = np->urng_aux.params;
   if (juaux) {
      gp->env = env;
      gp->junifArray = juaux;
      uaux = gp->unifArray = (*env)->GetDoubleArrayElements (env, juaux, 0);
      gp->n = gp->nextIndex = n;
      np->urng.params->unif = arr;
   }
   else {
      np->urng_aux.params = np->urng.params;
      uaux = 0;
   }

   v = (*env)->GetDoubleArrayElements (env, jv, 0);
   for (i = 0; i < n; i++)
      v[start+i] = unur_sample_cont (np->unurgen);
   (*env)->ReleaseDoubleArrayElements (env, jv, v, 0);
   if (juaux) {
      (*env)->ReleaseDoubleArrayElements (env, juaux, uaux, 0);
   }
   else {
      np->urng_aux.params = gp;
   }
   (*env)->ReleaseDoubleArrayElements (env, ju, u, 0);
}

JNIEXPORT jboolean JNICALL Java_umontreal_ssj_randvar_RandUnuran_isDiscrete
  (JNIEnv *env, jobject obj)
{
   const UNUR_DISTR* distr;
   struct NativeParams* np = (struct NativeParams*)(*env)->GetLongField
      (env, obj, fidNativeParams);

   if (!np)
      return JNI_FALSE;

   distr = unur_get_distr (np->unurgen);
   return unur_distr_is_discr (distr) ? JNI_TRUE : JNI_FALSE;
}

JNIEXPORT jboolean JNICALL Java_umontreal_ssj_randvar_RandUnuran_isContinuous
  (JNIEnv *env, jobject obj)
{
   const UNUR_DISTR* distr;
   struct NativeParams* np = (struct NativeParams*)(*env)->GetLongField
      (env, obj, fidNativeParams);

   if (!np)
      return JNI_FALSE;

   distr = unur_get_distr (np->unurgen);
   return unur_distr_is_cont (distr) ? JNI_TRUE : JNI_FALSE;
}

JNIEXPORT jboolean JNICALL Java_umontreal_ssj_randvar_RandUnuran_isContinuousMultivariate
  (JNIEnv *env, jobject obj)
{
   const UNUR_DISTR* distr;
   struct NativeParams* np = (struct NativeParams*)(*env)->GetLongField
      (env, obj, fidNativeParams);

   if (!np)
      return JNI_FALSE;

   distr = unur_get_distr (np->unurgen);
   return unur_distr_is_cvec (distr) ? JNI_TRUE : JNI_FALSE;
}

JNIEXPORT jboolean JNICALL Java_umontreal_ssj_randvar_RandUnuran_isEmpirical
  (JNIEnv *env, jobject obj)
{
   const UNUR_DISTR* distr;
   struct NativeParams* np = (struct NativeParams*)(*env)->GetLongField
      (env, obj, fidNativeParams);

   if (!np)
      return JNI_FALSE;

   distr = unur_get_distr (np->unurgen);
   return unur_distr_is_cemp (distr) ? JNI_TRUE : JNI_FALSE;
}

JNIEXPORT jboolean JNICALL Java_umontreal_ssj_randvar_RandUnuran_isEmpiricalMultivariate
  (JNIEnv *env, jobject obj)
{
   const UNUR_DISTR* distr;
   struct NativeParams* np = (struct NativeParams*)(*env)->GetLongField
      (env, obj, fidNativeParams);

   if (!np)
      return JNI_FALSE;

   distr = unur_get_distr (np->unurgen);
   return unur_distr_is_cvemp (distr) ? JNI_TRUE : JNI_FALSE;
}
