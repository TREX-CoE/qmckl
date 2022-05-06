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

AC_INIT([qmckl],[0.1.1],[https://github.com/TREX-CoE/qmckl/issues],[],[https://trex-coe.github.io/qmckl/index.html])
AC_CONFIG_AUX_DIR([tools])
AM_INIT_AUTOMAKE([subdir-objects color-tests parallel-tests silent-rules 1.11])

# Activate developer mode when the source is the git repository.
# Otherwise, it is the source distribution and the developer mode should not be activated.
AS_IF([test -d ${srcdir}/.git],
        [enable_maintainer_mode="yes"],
        [enable_maintainer_mode="no"]
)

AC_ARG_WITH(ifort, [AS_HELP_STRING([--with-ifort],[Use Intel Fortran compiler])], with_ifort=$withval, with_ifort=no)
AS_IF([test "$with_ifort" == "yes"], [
  FC=ifort
  FCFLAGS="-march=native -ip -Ofast -ftz -finline -g -mkl=sequential" ])

AC_ARG_WITH(icc, [AS_HELP_STRING([--with-icc],[Use Intel C compiler])], with_icc=$withval, with_icc=no)
AS_IF([test "$with_icc" == "yes"], [
  CC=icc
  CFLAGS="-march=native -ip -Ofast -ftz -finline -g -mkl=sequential" ])

AS_IF([test "$with_icc"."$with_ifort" == "yes.yes"], [
  ax_blas_ok="yes"
  ax_lapack_ok="yes"
  BLAS_LIBS=""
  LAPACK_LIBS=""])

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

#AM_SILENT_RULES(yes)
#AC_SUBST(SHARED_VERSION_INFO)
#AC_ENABLE_SHARED([no]) dnl shared libs cause too many headaches to be default

AC_LANG(C)

# Checks for programs.
AC_PROG_CC
AC_PROG_F77

# Make sure the c compiler supports C99
m4_version_prereq([2.70],[], [AC_PROG_CC_C99])
AS_IF([test "$ac_cv_prog_cc_c99" = "no"], [AC_MSG_ERROR([The compiler does not support C99])])
AC_PROG_CC_C_O
AM_PROG_CC_C_O
AC_PROG_FC
AC_PROG_FC_C_O
AC_FC_PP_DEFINE
AC_FC_SRCEXT([F90])
AC_FC_FREEFORM
LT_INIT
AC_PROG_INSTALL
AC_PROG_LN_S
PKG_PROG_PKG_CONFIG([])
PKG_LIBS=""
PKG_CFLAGS=""

# Checks for libraries.

AC_FC_LIBRARY_LDFLAGS

AC_CHECK_LIB([m], [sqrt])
AC_CHECK_HEADERS([assert.h errno.h math.h pthread.h stdbool.h stdint.h stdio.h stdlib.h string.h])

## pthread
AC_CHECK_LIB([pthread], [pthread_create])

# OpenMP
AC_ARG_WITH(openmp, [AS_HELP_STRING([--with-openmp],[activate OpenMP])], with_omp=$withval, with_omp=no)
if test "x$with_omp" = xyes; then
   AC_DEFINE([HAVE_OPENMP], [1], [Define to use OpenMP threading.])
   AX_OPENMP([],
    [AC_MSG_ERROR([Could not find OpenMP flags; configure with --without-openmp])])
   CFLAGS="${CFLAGS} ${OPENMP_CFLAGS}"
   FCFLAGS="${CFLAGS} ${OPENMP_FCFLAGS}"
fi

# CHAMELEON
AC_ARG_WITH(chameleon,
  AS_HELP_STRING([--without-chameleon],
                 [Do not use Chameleon. Default: auto-detect]), [
case "$with_chameleon" in
  no)
    : ;;
  yes)
    PKG_CHECK_MODULES([LIBCHAMELEON],[chameleon >= 1.0.0],
      [PKG_CFLAGS="$PKG_CFLAGS $LIBCHAMELEON_CFLAGS"
      PKG_LIBS="$PKG_LIBS $LIBCHAMELEON_LIBS"]
                      ,[

      ## something went wrong.
      ## try to find the package without pkg-config

      ## check that the library is actually new enough.
      ## by testing for a 1.0.0+ function which we use
      AC_CHECK_LIB(chameleon,CHAMELEON_finalize,[LIBCHAMELEON_LIBS="-lchameleon"])
    ])
    AH_TEMPLATE(HAVE_CHAMELEON,[Chameleon support is available])
  ;;
  *)
    if test ! -d "$withval" ; then
      AC_MSG_ERROR([--with-chameleon path does not point to a directory])
    fi
    LIBCHAMELEON_LIBS="-L$with_chameleon/lib -lchameleon -lchameleon_starpu -lhqr -lcoreblas"
    LIBCHAMELEON_CFLAGS="-I$with_chameleon/include $CFLAGS"
    LIBCHAMELEON_CPPFLAGS="-I$with_chameleon/include $CPPFLAGS"
  esac
])

if test "x$LIBCHAMELEON_LIBS" != "x" ; then
  LIBS="$LIBS $LIBCHAMELEON_LIBS"
  CFLAGS="$CFLAGS $LIBCHAMELEON_CFLAGS"
  CPPFLAGS="$CPPFLAGS $LIBCHAMELEON_CPPFLAGS"
  AC_CHECK_HEADERS([chameleon.h], [], [AC_MSG_ERROR("chamelon.h not found")])
  AC_DEFINE_UNQUOTED([HAVE_CHAMELEON],1,[CHAMELEON support is available])
fi

#AS_IF([test "x$with_chameleon" != "xno"], [
#])

AC_MSG_NOTICE([CHAMELEON library support: ${with_CHAMELEON:=auto} ${LIBCHAMELEON_PATH} ${LIBCHAMELEON_LIBS}])


# TREXIO

AC_ARG_WITH(trexio, [AS_HELP_STRING([--without-trexio],[disable support for TREXIO])],
            with_trexio=$withval, with_trexio=yes)

AS_IF([test "x$with_trexio" != xno], [
  AC_DEFINE([HAVE_TREXIO], [1], [Define if your have libtrexio])
  ARGS="${ARGS} trexio"
  AC_CHECK_LIB([trexio], [trexio_open],
       [],
       [AS_IF([test "x$with_trexio" != xcheck],
       [PKG_CHECK_MODULES([TREXIO], [trexio]) ])
  ])
])

PKG_CFLAGS="$PKG_CFLAGS $TREXIO_CFLAGS"
PKG_LIBS="$PKG_LIBS $TREXIO_LIBS"

## BLAS
AX_BLAS([], [AC_MSG_ERROR([BLAS was not found.])])

## LAPACK
AX_LAPACK([], [AC_MSG_ERROR([LAPACK was not found.])])
AS_IF([test "$BLAS_LIBS" == "$LAPACK_LIBS"], [BLAS_LIBS=""])


# Specific options required with some compilers
case $FC in

 *ifort*)
      FCFLAGS="$FCFLAGS -nofor-main"
      ;;

 *nvfortran*)
      FCFLAGS="$FCFLAGS -fPIC -Mnomain"
      ;;

esac

case $CC in

  *gcc*)
        CFLAGS="$CFLAGS -fPIC"
        ;;
  *nvc*)
        CFLAGS="$CFLAGS -fPIC"
        ;;
esac


# Options.
AC_ARG_ENABLE(hpc, [AS_HELP_STRING([--enable-hpc],[Use HPC-optimized functions])], HAVE_HPC=$enableval, HAVE_HPC=no)
AS_IF([test "$HAVE_HPC" = "yes"], [
   AC_DEFINE([HAVE_HPC], [1], [If defined, activate HPC routines])
])

# Enable Verificarlo tests
AC_ARG_ENABLE([vfc_ci],
[  --enable-vfc_ci    Build the library with vfc_ci support],
[case "${enableval}" in
  yes) vfc_ci=true &&  FCFLAGS="-D VFC_CI $FCFLAGS" && CFLAGS="-D VFC_CI $CFLAGS";;
  no)  vfc_ci=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable_vfc_ci]) ;;
esac],[vfc_ci=false])
AM_CONDITIONAL([VFC_CI], [test x$vfc_ci = xtrue])

if test "$FC" = "verificarlo-f"; then
  AC_MSG_NOTICE(verificarlo-f detected)
  # Arguments order is important here
  FCFLAGS="-Mpreprocess $FCFLAGS"
fi

## Enable GPU offloading

# GPU offloading
AC_ARG_ENABLE(gpu, [AS_HELP_STRING([--enable-gpu],[openmp|openacc : Use GPU-offloaded functions])], enable_gpu=$enableval, enable_gpu=no)
AS_IF([test "$enable_gpu" = "yes"], [enable_gpu="openmp"])

# OpenMP offloading
HAVE_OPENMP_OFFLOAD="no"
AS_IF([test "$enable_gpu" = "openmp"], [
  AC_DEFINE([HAVE_OPENMP_OFFLOAD], [1], [If defined, activate OpenMP-offloaded routines])
  HAVE_OPENMP_OFFLOAD="yes"
  case $CC in

    *gcc*)
          CFLAGS="$CFLAGS -fopenmp"
          ;;
    *nvc*)
          CFLAGS="$CFLAGS -mp=gpu"
          ;;
  esac

  case $FC in

    *gfortran*)
          FCFLAGS="$FCFLAGS -fopenmp"
          ;;
    *nvfortran*)
          FCFLAGS="$FCFLAGS -mp=gpu"
          ;;
  esac]
)

# OpenMP offloading
HAVE_OPENACC_OFFLOAD="no"
AS_IF([test "$enable_gpu" = "openacc"], [
  AC_DEFINE([HAVE_OPENACC_OFFLOAD], [1], [If defined, activate OpenACC-offloaded routines])
  HAVE_OPENACC_OFFLOAD="yes"
  case $CC in

    *gcc*)
          CFLAGS="$CFLAGS -fopenacc"
          ;;
    *nvc*)
          CFLAGS="$CFLAGS -acc=gpu"
          ;;
  esac

  case $FC in

    *gfortran*)
          FCFLAGS="$FCFLAGS -fopenacc"
          ;;
    *nvfortran*)
          FCFLAGS="$FCFLAGS -acc=gpu"
          ;;
  esac

])

# cuBLAS offloading
AC_ARG_WITH(cublas, [AS_HELP_STRING([--with-cublas],[Use cuBLAS-offloaded functions])], HAVE_CUBLAS_OFFLOAD=$withval, HAVE_CUBLAS_OFFLOAD=no)
AS_IF([test "$HAVE_CUBLAS_OFFLOAD" = "yes"], [
  AC_DEFINE([HAVE_CUBLAS_OFFLOAD], [1], [If defined, activate cuBLAS-offloaded routines])
  HAVE_OPENACC_OFFLOAD="yes"
  case $CC in

    *gcc*)
          CFLAGS="$CFLAGS -fopenmp"
          LDFLAGS="-lcublas"
          ;;
    *nvc*)
          CFLAGS="$CFLAGS -mp=gpu -cudalib=cublas"
          ;;
  esac

  case $FC in

    *gfortran*)
          FCFLAGS="$FCFLAGS -fopenmp"
          ;;
    *nvfortran*)
          FCFLAGS="$FCFLAGS -mp=gpu -cudalib=cublas"
          ;;
  esac
])

AC_ARG_ENABLE(malloc-trace, [AS_HELP_STRING([--enable-malloc-trace],[use debug malloc/free])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
	AC_DEFINE(MALLOC_TRACE,"malloc_trace.dat",[Define to use debugging malloc/free])
        ARGS="${ARGS} malloc-trace"
fi

AC_ARG_ENABLE(prof, [AS_HELP_STRING([--enable-prof],[compile for profiling])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
	CFLAGS="${CFLAGS} -pg"
	AC_DEFINE(ENABLE_PROF,1,[Define when using the profiler tool])
        ARGS="${ARGS} prof"
fi

AC_ARG_WITH(efence, [AS_HELP_STRING([--with-efence],[use ElectricFence library])], ok=$withval, ok=no)
if test "$ok" = "yes"; then
	AC_CHECK_LIB([efence], [malloc])
        ARGS="${ARGS} efence"
fi



##

AC_ARG_ENABLE(debug, [AS_HELP_STRING([--enable-debug],[compile for debugging])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
        if test "$GCC" = "yes"; then
	        CFLAGS="$CFLAGS \
-g -Wall -W -Wbad-function-cast -Wcast-qual -Warray-bounds -Wdisabled-optimization \
-Wpointer-arith -Wcast-align -Wpedantic -Wextra -Walloc-zero -Werror \
"
        fi
        if test "$GFC" = "yes"; then
	        FCFLAGS="$FCFLAGS \
-g -fcheck=all -Waliasing -Wampersand -Wconversion \
-Wsurprising -ffpe-trap=zero,overflow,underflow \
-Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation \
-Wreal-q-constant -Wuninitialized  -fbacktrace -finit-real=nan"
        fi
        if test "$FC" = "ifort"; then
	        FCFLAGS="$FCFLAGS \
-traceback -check all -debug all -fpe-all=0  -implicitnone"
        fi
	AC_DEFINE(DEBUG,1,[Define to turn on debugging checks])
        ARGS="${ARGS} debug"
fi

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
if test "x$enable_maintainer_mode" == "xyes" ; then
  QMCKL_DEVEL=" -- Developer mode"
else
  QMCKL_DEVEL=""
fi

AM_CONDITIONAL([QMCKL_DEVEL],[test "x$QMCKL_DEVEL" != x])

if test "x${QMCKL_DEVEL}" != "x"; then

  AC_PROG_AWK
  AM_PATH_PYTHON

  AC_CHECK_PROGS([EMACS],[emacs26 emacs],[no])
  if test x${EMACS} == xno ; then
    AC_MSG_ERROR([
    --------------------------------------
    Error: Emacs is required for org-mode.
    --------------------------------------
    ])
  fi

  AC_CHECK_PROGS([HAS_CPPCHECK],[cppcheck],[no])
  if test x${HAS_CPPCHECK} != xno ; then
        HAS_CPPCHECK=1
  fi

fi


#PKG-CONFIG
#mkl-dynamic-lp64-seq

PKG_LIBS="$PKG_LIBS $LIBS"
LIBS="$BLAS_LIBS $LAPACK_LIBS $PKG_LIBS"
CFLAGS="$CFLAGS $PKG_CFLAGS"
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
CPPFLAGS........: ${CPPFLAGS}
CFLAGS..........: ${CFLAGS}
FC..............: ${FC}
FCFLAGS.........: ${FCFLAGS}
LDFLAGS:........: ${LDFLAGS}
LIBS............: ${LIBS}
USE CHAMELEON...: ${with_chameleon}
HPC version.....: ${HAVE_HPC}
OpenMP offload..: ${HAVE_OPENMP_OFFLOAD}
OpenACC offload.: ${HAVE_OPENACC_OFFLOAD}
cuBLAS offload..: ${HAVE_CUBLAS_OFFLOAD}

Package features:
  ${ARGS}

Now type 'make @<:@<target>@:>@'
where the optional <target> is:
  all           - build the library
  check         - run tests
  install       - install ${PACKAGE_NAME}
--------------------------------------------------"
