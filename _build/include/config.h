/* include/config.h.  Generated from config.h.in by configure.  */
/* include/config.h.in.  Generated from configure.ac by autoheader.  */

/* VECTOR ALIGNED pragma */
#define ALIGNED _Pragma("vector aligned")

/* Define to turn on debugging checks */
/* #undef DEBUG */

/* Define when using the profiler tool */
/* #undef ENABLE_PROF */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to 1 if your Fortran compiler doesn't accept -c and -o together. */
/* #undef FC_NO_MINUS_C_MINUS_O */

/* Define to 1 if you have the <assert.h> header file. */
#define HAVE_ASSERT_H 1

/* Define if you have a BLAS library. */
#define HAVE_BLAS 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <errno.h> header file. */
#define HAVE_ERRNO_H 1

/* Activate floating-point exceptions */
/* #undef HAVE_FPE */

/* Activate HPC routines */
/* #undef HAVE_HPC */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if you have LAPACK library. */
#define HAVE_LAPACK 1

/* Define to 1 if you have the `efence' library (-lefence). */
/* #undef HAVE_LIBEFENCE */

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the `pthread' library (-lpthread). */
#define HAVE_LIBPTHREAD 1

/* Define to 1 if you have the `qmckldgemm' library (-lqmckldgemm). */
/* #undef HAVE_LIBQMCKLDGEMM */

/* Define to 1 if you have the `trexio' library (-ltrexio). */
/* #undef HAVE_LIBTREXIO */

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H 1

/* Define to 1 if you have the <math.h> header file. */
#define HAVE_MATH_H 1

/* Define to 1 if you have the `memset' function. */
#define HAVE_MEMSET 1

/* Define if OpenMP is enabled */
#define HAVE_OPENMP 1

/* Define to 1 if you have the posix_memalign function. */
#define HAVE_POSIX_MEMALIGN 1

/* Define to 1 if you have the <pthread.h> header file. */
#define HAVE_PTHREAD_H 1

/* Define if QMCKL DGEMM support is available */
/* #undef HAVE_QMCKLDGEMM */

/* Define to 1 if you have the <qmckl_dgemm.h> header file. */
/* #undef HAVE_QMCKL_DGEMM_H */

/* Define to 1 if you have the <stdbool.h> header file. */
#define HAVE_STDBOOL_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `strerror' function. */
#define HAVE_STRERROR 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define if TREXIO support is available */
/* #undef HAVE_TREXIO */

/* Define to 1 if you have the <trexio.h> header file. */
/* #undef HAVE_TREXIO_H */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if the system has the type `_Bool'. */
#define HAVE__BOOL 1

/* IVDEP pragma */
#define IVDEP _Pragma("ivdep")

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Define to use debugging malloc/free */
/* #undef MALLOC_TRACE */

/* Define to 1 if your C compiler doesn't accept -c and -o together. */
/* #undef NO_MINUS_C_MINUS_O */

/* Name of package */
#define PACKAGE "qmckl"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "https://github.com/TREX-CoE/qmckl/issues"

/* Define to the full name of this package. */
#define PACKAGE_NAME "qmckl"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "qmckl 1.0.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "qmckl"

/* Define to the home page for this package. */
#define PACKAGE_URL "https://trex-coe.github.io/qmckl/index.html"

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.0.0"

/* major version */
#define QMCKL_VERSION_MAJOR 1

/* minor version */
#define QMCKL_VERSION_MINOR 0

/* patch version */
#define QMCKL_VERSION_PATCH 0

/* Length of SIMD vectors */
#define SIMD_LENGTH 2

/* Define to 1 if all of the C90 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "1.0.0"

/* Define for Solaris 2.5.1 so the uint32_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT32_T */

/* Define for Solaris 2.5.1 so the uint64_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT64_T */

/* Define to the type of a signed integer type of width exactly 32 bits if
   such a type exists and the standard includes do not define it. */
/* #undef int32_t */

/* Define to the type of a signed integer type of width exactly 64 bits if
   such a type exists and the standard includes do not define it. */
/* #undef int64_t */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to the type of an unsigned integer type of width exactly 32 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint32_t */

/* Define to the type of an unsigned integer type of width exactly 64 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint64_t */
