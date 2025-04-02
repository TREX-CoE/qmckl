#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
# QMCkl - Quantum Monte Carlo kernel library
#
# BSD 3-Clause License
#
# Copyright (c) 2020, TREX Center of Excellence
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

AC_PREREQ([2.69])

AC_INIT([qmckl],[1.0.0],[https://github.com/TREX-CoE/qmckl/issues],[],[https://trex-coe.github.io/qmckl/index.html])
AC_CONFIG_AUX_DIR([tools])
AM_INIT_AUTOMAKE([subdir-objects color-tests parallel-tests silent-rules 1.11])

# Activate developer mode when the .maintainer_mode file is found
# Otherwise, it is the source distribution and the developer mode should not be activated.
AS_IF([test -f ${srcdir}/.maintainer_mode],
        [enable_maintainer_mode="yes"],
        [enable_maintainer_mode="no"]
)

# Intel Fortran compilers
AC_ARG_WITH([ifort],
            [AS_HELP_STRING([--with-ifort],
                            [Use Intel ifort Fortran compiler])],
            [with_ifort=$withval],
            [with_ifort=no])

AS_IF([test "x$with_ifort" = "xyes"], [
  FC=ifort
  FCFLAGS="-march=native -ip -O3 -ftz -finline -g -qmkl=sequential" ]
  ax_blas_ok="yes"
  ax_lapack_ok="yes"
  BLAS_LIBS=""
  LAPACK_LIBS=""
  ])

AC_ARG_WITH([ifx],
            [AS_HELP_STRING([--with-ifx],
                            [Use Intel ifx Fortran compiler])],
            [with_ifx=$withval],
            [with_ifx=no])

AS_IF([test "x$with_ifx" = "xyes"], [
  FC=ifx
  FCFLAGS="-march=native -O3 -ftz -g -qmkl=sequential"
  ax_blas_ok="yes"
  ax_lapack_ok="yes"
  BLAS_LIBS=""
  LAPACK_LIBS=""
  ])


AC_ARG_WITH([icx],
            [AS_HELP_STRING([--with-icx],
                            [Use Intel icx C compiler])],
            [with_icx=$withval],
            [with_icx=no])

AS_IF([test "x$with_icx" = "xyes"], [
  CC=icx
  CFLAGS="-march=native -O3 -ftz -finline -g -qmkl=sequential"
  ax_blas_ok="yes"
  ax_lapack_ok="yes"
  BLAS_LIBS=""
  LAPACK_LIBS=""
  ])

AC_ARG_WITH([icc],
            [AS_HELP_STRING([--with-icc],
                            [Use Intel icc C compiler])],
            [with_icc=$withval],
            [with_icc=no])

AS_IF([test "x$with_icc" = "xyes"], [
  CC=icc
  CFLAGS="-march=native -ip -O3 -ftz -finline -g -mkl=sequential"
  ax_blas_ok="yes"
  ax_lapack_ok="yes"
  BLAS_LIBS=""
  LAPACK_LIBS=""
  ])

# Default Automake variables
F77=$FC
AM_PROG_AR
AM_MAINTAINER_MODE()
LT_INIT
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_HEADERS([include/config.h])
AC_CONFIG_MACRO_DIR([m4])

VERSION_MAJOR=`echo ${PACKAGE_VERSION} | cut -d. -f1`
VERSION_MINOR=`echo ${PACKAGE_VERSION} | cut -d. -f2`
VERSION_PATCH=`echo ${PACKAGE_VERSION} | cut -d. -f3 | cut -d- -f1`
AC_DEFINE_UNQUOTED([QMCKL_VERSION_MAJOR], [$VERSION_MAJOR], [major version])
AC_DEFINE_UNQUOTED([QMCKL_VERSION_MINOR], [$VERSION_MINOR], [minor version])
AC_DEFINE_UNQUOTED([QMCKL_VERSION_PATCH], [$VERSION_PATCH], [patch version])
AC_SUBST([VERSION_MAJOR])
AC_SUBST([VERSION_MINOR])
AC_SUBST([VERSION_PATCH])

AM_SILENT_RULES(yes)

AC_LANG(C)

# Checks for programs.
AC_PROG_CC
AC_PROG_F77

# Make sure the c compiler supports C99
m4_version_prereq([2.70],[], [AC_PROG_CC_C99])
AS_IF([test "x$ac_cv_prog_cc_c99" = "xno"], [AC_MSG_ERROR([The compiler does not support C99])])
AC_PROG_CC_C_O
AM_PROG_CC_C_O
AC_PROG_FC
AC_PROG_FC_C_O
AC_FC_PP_DEFINE
AC_FC_SRCEXT([F90])
AC_FC_FREEFORM
AC_PROG_INSTALL
AC_PROG_LN_S
PKG_PROG_PKG_CONFIG([])

# Checks for libraries.
PKG_LIBS=""
PKG_CFLAGS=""
PKG_LIBS2=""
PKG_CFLAGS2=""

AC_FC_LIBRARY_LDFLAGS
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_HEADERS([assert.h errno.h math.h pthread.h stdbool.h stdint.h stdio.h stdlib.h string.h])

## pthread
AC_CHECK_LIB([pthread], [pthread_create])

# OpenMP
AC_ARG_WITH(openmp, [AS_HELP_STRING([--without-openmp],[deactivate OpenMP])], with_omp=$withval, with_omp=yes)
AS_IF([test "x$with_omp" = xyes], [
   AC_DEFINE([HAVE_OPENMP], [1], [Define to use OpenMP threading.])
   AX_OPENMP([],
    [AC_MSG_ERROR([Could not find OpenMP flags; configure with --without-openmp])])
   CFLAGS="${CFLAGS} ${OPENMP_CFLAGS}"
   FCFLAGS="${CFLAGS} ${OPENMP_FCFLAGS}"
])

# CHAMELEON
#AC_ARG_WITH(chameleon,
#  AS_HELP_STRING([--without-chameleon],
#                 [Do not use Chameleon. Default: auto-detect]),
#                 [case "$with_chameleon" in
#                    no)
#                        : ;;
#                    yes)
#                        PKG_CHECK_MODULES([LIBCHAMELEON],[chameleon >= 1.0.0],
#                        [PKG_CFLAGS="$PKG_CFLAGS $LIBCHAMELEON_CFLAGS"
#                        PKG_LIBS2="$PKG_LIBS2 $LIBCHAMELEON_LIBS"]
#                                        ,[
#
#                        ## something went wrong.
#                        ## try to find the package without pkg-config
#
#                        ## check that the library is actually new enough.
#                        ## by testing for a 1.0.0+ function which we use
#                        AC_CHECK_LIB(chameleon,CHAMELEON_finalize,[LIBCHAMELEON_LIBS="-lchameleon"])
#                        ])
#                        AH_TEMPLATE(HAVE_CHAMELEON,[Chameleon support is available])
#                    ;;
#                    *)
#                        if test ! -d "$withval" ; then
#                        AC_MSG_ERROR([--with-chameleon argument does not point to a directory])
#                        fi
#                        LIBCHAMELEON_LIBS="-L$with_chameleon/lib -lchameleon -lchameleon_starpu -lhqr -lcoreblas"
#                        LIBCHAMELEON_CFLAGS="-I$with_chameleon/include $CFLAGS"
#                        LIBCHAMELEON_CPPFLAGS="-I$with_chameleon/include $CPPFLAGS"
#                    esac
#                    ])
#
#if test "x$LIBCHAMELEON_LIBS" != "x" ; then
#  LIBS="$LIBS $LIBCHAMELEON_LIBS"
#  CFLAGS="$CFLAGS $LIBCHAMELEON_CFLAGS"
#  CPPFLAGS="$CPPFLAGS $LIBCHAMELEON_CPPFLAGS"
#  AC_CHECK_HEADERS([chameleon.h], [], [AC_MSG_ERROR("chamelon.h not found")])
#  AC_DEFINE_UNQUOTED([HAVE_CHAMELEON],1,[CHAMELEON support is available])
#fi
#
#AC_MSG_NOTICE([CHAMELEON library support: ${with_CHAMELEON:=auto} ${LIBCHAMELEON_PATH} ${LIBCHAMELEON_LIBS}])


# TREXIO
AC_ARG_WITH([trexio],
            [AS_HELP_STRING([--without-trexio],
                            [disable support for TREXIO])],
            [with_trexio=$withval],
            [with_trexio=yes])

AS_CASE([$with_trexio],
  [no], [],
  [yes],
  [AC_CHECK_LIB([trexio], [trexio_open], [],
                [PKG_CHECK_MODULES([TREXIO], [trexio]) ])
   CFLAGS="$CFLAGS $TREXIO_CFLAGS"
   LIBS="$LIBS $TREXIO_LIBS"
   AC_CHECK_HEADERS([trexio.h], [],
                    [AC_MSG_ERROR("trexio.h not found")])
   AC_DEFINE([HAVE_TREXIO], [1], [Define if TREXIO support is available])
   ARGS="${ARGS} trexio"],
  [AS_IF([test ! -d "$withval"],
         [AC_MSG_ERROR([--with-trexio=$withval: argument is not a directory])])
   TREXIO_LIBS="-L$with_trexio/lib -ltrexio"
   TREXIO_CFLAGS="-I$with_trexio/include"
   CFLAGS="$CFLAGS $TREXIO_CFLAGS"
   LIBS="$LIBS $TREXIO_LIBS"
   AC_CHECK_HEADERS([trexio.h], [],
                    [AC_MSG_ERROR("trexio.h not found")])
   AC_DEFINE([HAVE_TREXIO], [1], [Define if TREXIO support is available])
   ARGS="${ARGS} trexio"
  ])

PKG_CFLAGS2="$PKG_CFLAGS2 $TREXIO_CFLAGS"
PKG_LIBS2="$PKG_LIBS2 $TREXIO_LIBS"


# QMCKLDGEMM
AC_ARG_WITH([qmckldgemm],
            [AS_HELP_STRING([--with-qmckldgemm],
                            [enable support for QMCKL DGEMM])],
            [with_qmckldgemm=$withval],
            [with_qmckldgemm=no])

AS_CASE([$with_qmckldgemm],
  [no], [],
  [yes],
  [AC_CHECK_LIB([qmckldgemm], [qmckl_packed_matrix_create], [],
                [PKG_CHECK_MODULES([QMCKLDGEMM], [qmckldgemm]) ])
   CFLAGS="$CFLAGS $QMCKLDGEMM_CFLAGS"
   LIBS="$LIBS $QMCKLDGEMM_LIBS"
   AC_CHECK_HEADERS([qmckl_dgemm.h], [],
                    [AC_MSG_ERROR("qmckl_dgemm.h not found")])
   AC_DEFINE([HAVE_QMCKLDGEMM], [1], [Define if QMCKL DGEMM support is available])
   ARGS="${ARGS} qmckldgemm"],
  [AS_IF([test ! -d "$withval"],
         [AC_MSG_ERROR([--with-qmckldgemm=$withval: argument is not a directory])])
   QMCKLDGEMM_LIBS="-L$with_qmckldgemm/lib -lqmckldgemm"
   QMCKLDGEMM_CFLAGS="-I$with_qmckldgemm/include"
   CFLAGS="$CFLAGS $QMCKLDGEMM_CFLAGS"
   LIBS="$LIBS $QMCKLDGEMM_LIBS"
   echo $CFLAGS
   AC_CHECK_HEADERS([qmckl_dgemm.h], [],
                    [AC_MSG_ERROR("qmckl_dgemm.h not found")])
   AC_DEFINE([HAVE_QMCKLDGEMM], [1], [Define if QMCKL DGEMM support is available])
   ARGS="${ARGS} qmckldgemm"
  ])

PKG_CFLAGS2="$PKG_CFLAGS2 $QMCKLDGEMM_CFLAGS"
PKG_LIBS2="$PKG_LIBS2 $QMCKLDGEMM_LIBS"

## BLAS
AX_BLAS([], [AC_MSG_ERROR([BLAS was not found.])])

## LAPACK
AX_LAPACK([], [AC_MSG_ERROR([LAPACK was not found.])])
AS_IF([test "$BLAS_LIBS" = "$LAPACK_LIBS"], [BLAS_LIBS=""])


# Specific options required with some compilers
AS_CASE([$FC],
 [*ifort*], [FCFLAGS="$FCFLAGS -nofor-main"],
 [*nvfortran*], [FCFLAGS="$FCFLAGS -fPIC -Mnomain"],
 [*pgf90*], [FCFLAGS="$FCFLAGS -fPIC -Mnomain"],
 [*pgf77*], [FCFLAGS="$FCFLAGS -fPIC -Mnomain"],
 [])

AS_CASE([$CC],
  [*gcc*], [CFLAGS="$CFLAGS -fPIC"],
  [*nvc*], [CFLAGS="$CFLAGS -fPIC"],
  [*pgcc*], [CFLAGS="$CFLAGS -fPIC"],
  [])


# Options
AC_ARG_ENABLE([hpc],
              [AS_HELP_STRING([--enable-hpc],
                              [Use HPC-optimized functions])],
              [enable_hpc=$enableval],
              [enable_hpc=no])

AS_IF([test "x$enable_hpc" = "xyes"],
      [AC_DEFINE([HAVE_HPC], [1], [Activate HPC routines])])


AC_ARG_ENABLE([fpe],
              [AS_HELP_STRING([--enable-fpe],
                              [Enable floating-point exceptions])],
              [enable_fpe=$enableval],
              [enable_fpe=no])

AS_IF([test "x$enable_fpe" = "xyes"],
      [AC_DEFINE([HAVE_FPE], [1], [Activate floating-point exceptions])])



AC_ARG_ENABLE([doc],
            [AS_HELP_STRING([--disable-doc],
                            [Disable documentation])],
            [enable_doc=$enableval],
            [enable_doc=yes])

AM_CONDITIONAL([HAVE_DOC],[test "x$enable_doc" = "xyes"])


# Enable Verificarlo tests
AC_ARG_ENABLE([vfc_ci],
              [AS_HELP_STRING([--enable-vfc_ci],
                              [Build the library with vfc_ci support])],
              [enable_vfc_ci=$enableval],
              [enable_vfc_ci=no])

AS_IF([test "x$enable_vfc_ci" = "xyes"],
       [FCFLAGS="-D VFC_CI $FCFLAGS"
	CFLAGS="-D VFC_CI $CFLAGS"
        LIBS="$LIBS -lvfc_probes"])

AM_CONDITIONAL([VFC_CI], [test "x$vfc_ci" = "xtrue"])

AS_IF([test "$FC" = "verificarlo-f"], [
  AC_MSG_NOTICE(verificarlo-f detected)
  # Arguments order is important here
  FCFLAGS="-Mpreprocess $FCFLAGS"
])

AC_ARG_ENABLE(malloc-trace, [AS_HELP_STRING([--enable-malloc-trace],[use debug malloc/free])], ok=$enableval, ok=no)
AS_IF([test "x$ok" = "xyes"], [
	AC_DEFINE(MALLOC_TRACE,"malloc_trace.dat",[Define to use debugging malloc/free])
        ARGS="${ARGS} malloc-trace"
])

AC_ARG_ENABLE(prof, [AS_HELP_STRING([--enable-prof],[compile for profiling])], ok=$enableval, ok=no)
AS_IF([test "x$ok" = "xyes"], [
	CFLAGS="${CFLAGS} -pg"
	AC_DEFINE(ENABLE_PROF,1,[Define when using the profiler tool])
        ARGS="${ARGS} prof"
])

AC_ARG_WITH(efence, [AS_HELP_STRING([--with-efence],[use ElectricFence library])], ok=$withval, ok=no)
AS_IF([test "x$ok" = "xyes"], [
	AC_CHECK_LIB([efence], [malloc])
        ARGS="${ARGS} efence"
])


AC_ARG_ENABLE(sanitizer, [AS_HELP_STRING([--enable-sanitizer],[enable sanitizer debug flags])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  CFLAGS="${CFLAGS} -fsanitize=address -fsanitize=undefined -fsanitize=leak -fsanitize=pointer-compare -fsanitize=pointer-subtract -fsanitize=bounds -fsanitize=bounds-strict"
  FCFLAGS="${FCFLAGS} -fsanitize=address -fsanitize=undefined -fsanitize=leak -fsanitize=pointer-compare -fsanitize=pointer-subtract -fsanitize=bounds -fsanitize=bounds-strict"
fi



##

AC_ARG_ENABLE(debug, [AS_HELP_STRING([--enable-debug],[compile for debugging])], ok=$enableval, ok=no)
AS_IF([test "x$ok" = "xyes"], [
        AS_IF([test "x$GCC" = "xyes"], [
                CPPFLAGS="-Wdate-time -D_FORTIFY_SOURCE=2"
	        CFLAGS="$CFLAGS -g \
-Wall -W -Wbad-function-cast -Wcast-qual -Warray-bounds -Wdisabled-optimization \
-fno-omit-frame-pointer -fstack-protector-strong -Wformat -Werror=format-security \
-Wpointer-arith -Wcast-align -Wpedantic -Wextra -Walloc-zero -Werror \
"
		LDFLAGS="$LDFLAGS"
        ])
        AS_IF([test "x$GFC" = "xyes"], [
	        FCFLAGS="$FCFLAGS \
-g -fcheck=all -Waliasing -Wampersand -Wconversion \
-Wsurprising -ffpe-trap=zero,overflow,underflow \
-Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation \
-Wreal-q-constant -Wuninitialized  -fbacktrace -finit-real=nan"
        ])
        AS_IF([test "x$FC" = "xifort"], [
	        FCFLAGS="$FCFLAGS \
-traceback -check all -debug all -fpe-all=0  -implicitnone"
        ])
	AC_DEFINE(DEBUG,1,[Define to turn on debugging checks])
        ARGS="${ARGS} debug"
])

# Checks for header files.

## qmckl
AC_CHECK_HEADERS([inttypes.h malloc.h stdint.h stdlib.h string.h strings.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.

## qmckl
AC_CHECK_HEADER_STDBOOL
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T


# Checks for library functions.

## qmckl
# AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset strerror])

# Development mode
  QMCKL_DEVEL=
AS_IF([test "x$enable_maintainer_mode" = "xyes"],
      [QMCKL_DEVEL=" -- Developer mode"])

AM_CONDITIONAL([QMCKL_DEVEL],[test "x$QMCKL_DEVEL" != x])

AC_ARG_ENABLE([python],
            [AS_HELP_STRING([--enable-python],
                            [Enable building the Python binding])],
            [enable_python=$enableval],
            [enable_python=no])

AM_CONDITIONAL([HAVE_PYTHON],[test "x$enable_python" = "xyes"])

AS_IF([test "x${QMCKL_DEVEL}" != "x"], [
  AC_PROG_AWK
  AM_PATH_PYTHON

  AC_CHECK_PROGS([EMACS],[emacs26 emacs],[no])
  AS_IF([test "x${EMACS}" = "xno"],
    [AC_MSG_ERROR([
    --------------------------------------
    Error: Emacs is required for org-mode.
    --------------------------------------
    ])])

  AC_CHECK_PROGS([HAS_CPPCHECK],[cppcheck],[no])
  AS_IF([test "x${HAS_CPPCHECK}" != "xno"], [HAS_CPPCHECK=1])

  AS_IF([test "x$enable_python" = "xyes"], [
    AX_PKG_SWIG(4.0.0, [], AC_MSG_WARN([SWIG is required to build Python API.]) )
  ] )
])

# Checking SIMD
AC_MSG_CHECKING([for SIMD length])
SIMD_LENGTH=1
AC_RUN_IFELSE(
    [AC_LANG_PROGRAM([],[
int simd=1;
#if defined(__AVX512F__)
   simd=8;
#elif defined(__AVX2__)
   simd=4;
#elif defined(__AVX__)
   simd=4;
#elif defined(__SSE2__)
   simd=2;
#elif defined(__ARM_NEON__)
   simd=2;
#endif
  return simd;
])], [SIMD_LENGTH=1],
[ AS_CASE([$?],
  [1], [SIMD_LENGTH=1],
  [2], [SIMD_LENGTH=2],
  [4], [SIMD_LENGTH=4],
  [8], [SIMD_LENGTH=8],
  [16], [SIMD_LENGTH=16],
  [SIMD_LENGTH=1])],
  [SIMD_LENGTH=1]
)
AC_MSG_RESULT([$SIMD_LENGTH])
AC_DEFINE_UNQUOTED([SIMD_LENGTH], [$SIMD_LENGTH], [Length of SIMD vectors])

# Checking IVDEP
ivdep=""
AC_MSG_CHECKING([for ivdep pragma])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
    #include <stdio.h>
]], [[
    int main() {
        #pragma ivdep
        for (int i = 0; i < 10; ++i) {
            printf("Testing: %d\n", i);
        }
        return 0;
    }
]])],
    [ivdep='_Pragma("ivdep")'], [
])

AS_IF([test "x$ivdep" = "x"], [
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      #include <stdio.h>
  ]], [[
      int main() {
          #pragma clang loop vectorize(enable)
          for (int i = 0; i < 10; ++i) {
              printf("Testing: %d\n", i);
          }
          return 0;
      }
  ]])],
      [ivdep='_Pragma("clang loop vectorize(enable)")'], [
  ])
])

AS_IF([test "x$ivdep" = "x"], [
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      #include <stdio.h>
  ]], [[
      int main() {
          #pragma GCC ivdep
          for (int i = 0; i < 10; ++i) {
              printf("Testing: %d\n", i);
          }
          return 0;
      }
  ]])],
      [ivdep='_Pragma("GCC ivdep")'], [
  ])
])

AC_DEFINE_UNQUOTED([IVDEP], [$ivdep], [IVDEP pragma])
AS_IF([test "x$ivdep" = "x"], [
  ivdep="no"
])
AC_MSG_RESULT([$ivdep])


# Checking ALIGNED

AC_MSG_CHECKING([for posix_memalign])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <stdlib.h>
#include <stdio.h>
]], [[
    int main() {
      void *ptr;
      int ret = posix_memalign(&ptr, 64, 1024);
      if (ret != 0) {
        return EXIT_FAILURE;
      }
      free(ptr);
      return 0;
    }
]])],
    [have_posix_memalign=yes], [have_posix_memalign=no
])
AS_IF([test "x$have_posix_memalign" = "xyes"], [
  AC_DEFINE([HAVE_POSIX_MEMALIGN], [1], [Define to 1 if you have the posix_memalign function.])
])
AC_MSG_RESULT([$have_posix_memalign])

aligned=""
AC_MSG_CHECKING([for vector aligned pragma])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
]], [[
    int main() {
        double __attribute__((aligned(8))) a[10] ;
        #pragma vector aligned
        for (int i = 0; i < 10; ++i) {
            a[i] = (double) i;
        }
        return 0;
    }
]])],
    [aligned='_Pragma("vector aligned")'], [
])

AS_IF([test "x$have_posix_memalign" = "xno"], [
  aligned=""
])

AC_DEFINE_UNQUOTED([ALIGNED], [$aligned], [VECTOR ALIGNED pragma])
AS_IF([test "x$aligned" = "x"], [
  aligned="no"
])
AC_MSG_RESULT([$aligned])



#PKG-CONFIG
#mkl-dynamic-lp64-seq

PKG_LIBS="$PKG_LIBS $LIBS"
LIBS="$BLAS_LIBS $LAPACK_LIBS $PKG_LIBS"
CFLAGS="$CFLAGS $PKG_CFLAGS"
PKG_LIBS="$PKG_LIBS $PKG_LIBS2"
PKG_CFLAGS="$PKG_CFLAGS $PKG_CFLAGS2"
AC_SUBST([PKG_LIBS])
AC_SUBST([PKG_CFLAGS])
AC_SUBST([HAS_CPPCHECK])


AC_CONFIG_FILES([Makefile
                 pkgconfig/qmckl.pc
                ])


AC_OUTPUT

echo \
"-------------------------------------------------

${PACKAGE_NAME} Version ${PACKAGE_VERSION} ${QMCKL_DEVEL}

Prefix: '${prefix}'.

CC..............: ${CC}
CFLAGS..........: ${CFLAGS}
CPPFLAGS........: ${CPPFLAGS}
FC..............: ${FC}
FCFLAGS.........: ${FCFLAGS}
HPC version.....: ${HAVE_HPC}
LDFLAGS:........: ${LDFLAGS}
LIBS............: ${LIBS}
OPENMP..........: ${with_omp}
PYTHON..........: ${enable_python}
QMCKL_DGEMM.....: ${with_qmckldgemm}
TREXIO..........: ${with_trexio}

Package features:
  ${ARGS}

Now type 'make @<:@<target>@:>@'
where the optional <target> is:
  all           - build the library
  check         - run tests
  install       - install ${PACKAGE_NAME}
--------------------------------------------------"
