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

AC_INIT([qmckl],[0.3.1],[https://github.com/TREX-CoE/qmckl/issues],[],[https://trex-coe.github.io/qmckl/index.html])
AC_CONFIG_AUX_DIR([tools])
AM_INIT_AUTOMAKE([subdir-objects color-tests parallel-tests silent-rules 1.11])

# Activate developer mode when the source is the git repository.
# Otherwise, it is the source distribution and the developer mode should not be activated.
AS_IF([test -d ${srcdir}/.git],
        [enable_maintainer_mode="yes"],
        [enable_maintainer_mode="no"]
)

# Intel Fortran compiler
AC_ARG_WITH([ifort],
            [AS_HELP_STRING([--with-ifort],
                            [Use Intel Fortran compiler])],
            [with_ifort=$withval],
            [with_ifort=no])

AS_IF([test "$with_ifort" = "yes"], [
  FC=ifort
  FCFLAGS="-march=native -ip -Ofast -ftz -finline -g -mkl=sequential" ])

# Intel C compiler
AC_ARG_WITH([icc],
            [AS_HELP_STRING([--with-icc],
                            [Use Intel C compiler])],
            [with_icc=$withval],
            [with_icc=no])

AS_IF([test "$with_icc" = "yes"], [
  CC=icc
  CFLAGS="-march=native -ip -Ofast -ftz -finline -g -mkl=sequential" ])

AS_IF([test "$with_icc"."$with_ifort" = "yes.yes"], [
  ax_blas_ok="yes"
  ax_lapack_ok="yes"
  BLAS_LIBS=""
  LAPACK_LIBS=""])

# Default Automake variables
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
AS_IF([test "$ac_cv_prog_cc_c99" = "no"], [AC_MSG_ERROR([The compiler does not support C99])])
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
AC_ARG_WITH(openmp, [AS_HELP_STRING([--with-openmp],[activate OpenMP])], with_omp=$withval, with_omp=no)
if test "x$with_omp" = xyes; then
   AC_DEFINE([HAVE_OPENMP], [1], [Define to use OpenMP threading.])
   AX_OPENMP([],
    [AC_MSG_ERROR([Could not find OpenMP flags; configure with --without-openmp])])
   CFLAGS="${CFLAGS} ${OPENMP_CFLAGS}"
   FCFLAGS="${CFLAGS} ${OPENMP_FCFLAGS}"
fi

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
 [])

AS_CASE([$CC],
  [*gcc*], [CFLAGS="$CFLAGS -fPIC"],
  [*nvc*], [CFLAGS="$CFLAGS -fPIC"],
  [])


# Options
AC_ARG_ENABLE([hpc],
              [AS_HELP_STRING([--enable-hpc],
                              [Use HPC-optimized functions])],
              [enable_hpc=$enableval],
              [enable_hpc=no])

AS_IF([test "x$enable_hpc" = "yes"],
      [AC_DEFINE([HAVE_HPC], [1], [Activate HPC routines])])

AC_ARG_ENABLE([doc],
            [AS_HELP_STRING([--disable-doc],
                            [Disable documentation])],
            [enable_doc=$withval],
            [enable_doc=yes])

AM_CONDITIONAL([HAVE_DOC],[test "x$enable_doc" = "no"])

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
  AS_CASE([$CC],
    [*gcc*], [CFLAGS="$CFLAGS -fopenmp"],
    [*nvc*], [CFLAGS="$CFLAGS -mp=gpu"],
    [])
  AS_CASE([$FC],
    [*gfortran*],  [FCFLAGS="$FCFLAGS -fopenmp"],
    [*nvfortran*], [FCFLAGS="$FCFLAGS -mp=gpu"],
    [])
  ]
)

# OpenMP offloading
HAVE_OPENACC_OFFLOAD="no"
AS_IF([test "$enable_gpu" = "openacc"], [
  AC_DEFINE([HAVE_OPENACC_OFFLOAD], [1], [If defined, activate OpenACC-offloaded routines])
  HAVE_OPENACC_OFFLOAD="yes"
  AS_CASE([$CC],
    [*gcc*], [CFLAGS="$CFLAGS -fopenacc"],
    [*nvc*], [CFLAGS="$CFLAGS -acc=gpu"],
    [])
  AS_CASE([$FC],
    [*gfortran*],  [FCFLAGS="$FCFLAGS -fopenacc"],
    [*nvfortran*], [FCFLAGS="$FCFLAGS -acc=gpu"],
    [])
  ]

])

# cuBLAS offloading
AC_ARG_WITH(cublas, [AS_HELP_STRING([--with-cublas],[Use cuBLAS-offloaded functions])], HAVE_CUBLAS_OFFLOAD=$withval, HAVE_CUBLAS_OFFLOAD=no)
AS_IF([test "$HAVE_CUBLAS_OFFLOAD" = "yes"], [
  AC_DEFINE([HAVE_CUBLAS_OFFLOAD], [1], [If defined, activate cuBLAS-offloaded routines])
  HAVE_OPENACC_OFFLOAD="yes"
  AS_CASE([$CC],
    [*gcc*], [CFLAGS="$CFLAGS -fopenmp"
              LDFLAGS="-lcublas"],
    [*nvc*], [CFLAGS="$CFLAGS -mp=gpu -cudalib=cublas"],
    [])
  AS_CASE([$FC],
    [*gfortran*],  [FCFLAGS="$FCFLAGS -fopenmp"],
    [*nvfortran*], [FCFLAGS="$FCFLAGS -mp=gpu -cudalib=cublas"],
    [])
  ]
])

AC_ARG_ENABLE(malloc-trace, [AS_HELP_STRING([--enable-malloc-trace],[use debug malloc/free])], ok=$enableval, ok=no)
AS_IF([test "$ok" = "yes"], [
	AC_DEFINE(MALLOC_TRACE,"malloc_trace.dat",[Define to use debugging malloc/free])
        ARGS="${ARGS} malloc-trace"
])

AC_ARG_ENABLE(prof, [AS_HELP_STRING([--enable-prof],[compile for profiling])], ok=$enableval, ok=no)
AS_IF([test "$ok" = "yes"], [
	CFLAGS="${CFLAGS} -pg"
	AC_DEFINE(ENABLE_PROF,1,[Define when using the profiler tool])
        ARGS="${ARGS} prof"
])

AC_ARG_WITH(efence, [AS_HELP_STRING([--with-efence],[use ElectricFence library])], ok=$withval, ok=no)
AS_IF([test "$ok" = "yes"], [
	AC_CHECK_LIB([efence], [malloc])
        ARGS="${ARGS} efence"
])



##

AC_ARG_ENABLE(debug, [AS_HELP_STRING([--enable-debug],[compile for debugging])], ok=$enableval, ok=no)
AS_IF([test "$ok" = "yes"], [
        AS_IF([test "$GCC" = "yes"], [
	        CFLAGS="$CFLAGS \
-g -Wall -W -Wbad-function-cast -Wcast-qual -Warray-bounds -Wdisabled-optimization \
-fsanitize=address -fno-omit-frame-pointer \
-Wpointer-arith -Wcast-align -Wpedantic -Wextra -Walloc-zero -Werror \
"
		LDFLAGS="$LDFLAGS -fsanitize=address"
        ])
        AS_IF([test "$GFC" = "yes"], [
	        FCFLAGS="$FCFLAGS \
-g -fcheck=all -Waliasing -Wampersand -Wconversion \
-Wsurprising -ffpe-trap=zero,overflow,underflow \
-Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation \
-Wreal-q-constant -Wuninitialized  -fbacktrace -finit-real=nan"
        ])
        AS_IF([test "$FC" = "ifort"], [
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

  AX_PKG_SWIG(4.0.0, [], AC_MSG_WARN([SWIG is required to build Python API.]) )
])


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
CPPFLAGS........: ${CPPFLAGS}
CFLAGS..........: ${CFLAGS}
FC..............: ${FC}
FCFLAGS.........: ${FCFLAGS}
LDFLAGS:........: ${LDFLAGS}
LIBS............: ${LIBS}
TREXIO..........: ${with_trexio}
QMCKL_DGEMM.....: ${with_qmckldgemm}
HPC version.....: ${HAVE_HPC}

Package features:
  ${ARGS}

Now type 'make @<:@<target>@:>@'
where the optional <target> is:
  all           - build the library
  check         - run tests
  install       - install ${PACKAGE_NAME}
--------------------------------------------------"
