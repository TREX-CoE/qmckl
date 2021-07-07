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

AC_INIT([qmckl],[0.1.1],
        [https://github.com/TREX-CoE/qmckl/issues], [],
        [https://trex-coe.github.io/qmckl/index.html])
AC_CONFIG_AUX_DIR(tools)
AM_INIT_AUTOMAKE([subdir-objects color-tests parallel-tests silent-rules 1.11])
AM_MAINTAINER_MODE()
LT_INIT
AC_CONFIG_SRCDIR([Makefile.in])
AC_CONFIG_HEADERS([include/config.h])
AC_CONFIG_MACRO_DIR([m4])

VERSION_MAJOR=`echo ${PACKAGE_VERSION} | cut -d. -f1`
VERSION_MINOR=`echo ${PACKAGE_VERSION} | cut -d. -f2`
VERSION_PATCH=`echo ${PACKAGE_VERSION} | cut -d. -f3 | cut -d- -f1`
AC_DEFINE_UNQUOTED(QMCKL_VERSION_MAJOR, [$VERSION_MAJOR], [major version])
AC_DEFINE_UNQUOTED(QMCKL_VERSION_MINOR, [$VERSION_MINOR], [minor version])
AC_DEFINE_UNQUOTED(QMCKL_VERSION_PATCH, [$VERSION_PATCH], [patch version])
AC_SUBST([VERSION_MAJOR])
AC_SUBST([VERSION_MINOR])
AC_SUBST([VERSION_PATCH])

#AM_SILENT_RULES(yes)
#AC_SUBST(SHARED_VERSION_INFO)
#AM_ENABLE_SHARED(no) dnl shared libs cause too many headaches to be default

AC_LANG(C)

# Checks for programs.
AC_PROG_CC
# Make sure the c compiler supports C99
m4_version_prereq([2.70],[], [AC_PROG_CC_C99])
AS_IF([test "$ac_cv_prog_cc_c99" = "no"], [AC_MSG_ERROR([The compiler does not support C99])])
AC_PROG_CC_C_O
AC_PROG_FC
AC_PROG_FC_C_O
AC_FC_SRCEXT([f90])
AC_FC_FREEFORM
AC_PROG_LIBTOOL
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
#AC_ARG_WITH(openmp, [AC_HELP_STRING([--with-openmp],[enable OpenMP])], with_omp=$withval, with_omp=no)
#if test "x$with_omp" = xyes; then
#   AC_DEFINE([USE_OPENMP], [1], [Define to use OpenMP threading.])
#   AX_OPENMP([],
#    [AC_MSG_ERROR([Could not find OpenMP flags; configure with --without-openmp])])
#   CFLAGS="${CFLAGS} ${OPENMP_CFLAGS}"
#fi


## BLAS
#AX_BLAS([], [AC_MSG_ERROR([BLAS was not found.])])

## LAPACK
#AX_LAPACK([], [AC_MSG_ERROR([LAPACK was not found.])])

# Options.

AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],[compile for debugging])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
        if test "$GCC" = "yes"; then
	        CFLAGS="$CFLAGS \
-Wall -W -Wbad-function-cast -Wcast-qual \
-Wpointer-arith -Wcast-align -Wpedantic -Wextra -fmax-errors=3"
        fi
        if test "$GFC" = "yes"; then
	        FCFLAGS="$FCFLAGS \
-fcheck=all -Waliasing -Wampersand -Wconversion \
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

AC_ARG_ENABLE(malloc-trace, [AC_HELP_STRING([--enable-malloc-trace],[use debug malloc/free])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
	AC_DEFINE(MALLOC_TRACE,"malloc_trace.dat",[Define to use debugging malloc/free])
        ARGS="${ARGS} malloc-trace"
fi

AC_ARG_ENABLE(prof, [AC_HELP_STRING([--enable-prof],[compile for profiling])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
	CFLAGS="${CFLAGS} -pg"
	AC_DEFINE(ENABLE_PROF,1,[Define when using the profiler tool])
        ARGS="${ARGS} prof"
fi

AC_ARG_WITH(efence, [AC_HELP_STRING([--with-efence],[use ElectricFence library])], ok=$withval, ok=no)
if test "$ok" = "yes"; then
	AC_CHECK_LIB(efence, malloc)
        ARGS="${ARGS} efence"
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
AM_CONDITIONAL([QMCKL_DEVEL],[test "x$QMCKL_DEVEL" != x])

if test "x${QMCKL_DEVEL}" != "x"; then

  QMCKL_DEVEL=" -- Developer mode"
  AC_PROG_AWK
  AM_PATH_PYTHON
  ${PYTHON} ${srcdir}/tools/build_makefile.py

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

# Enable Verificarlo tests
AC_ARG_ENABLE([vfc_ci],
[  --enable-vfc_ci    Build the library with vfc_ci support],
[case "${enableval}" in
  yes) vfc_ci=true &&  FCFLAGS="-D VFC_CI $FCFLAGS" && CFLAGS="-D VFC_CI $CFLAGS";;
  no)  vfc_ci=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable_vfc_ci]) ;;
esac],[vfc_ci=false])
AM_CONDITIONAL([VFC_CI], [test x$vfc_ci = xtrue])

# Enable Fortran preprocessor
if test "$FC" = "gfortran"; then
  AC_MSG_NOTICE(gfortran detected)
  # Arguments order is important here
  FCFLAGS="-cpp $FCFLAGS"
fi

if test "$FC" = "verificarlo-f"; then
  AC_MSG_NOTICE(verificarlo-f detected)
  # Arguments order is important here
  FCFLAGS="-Mpreprocess $FCFLAGS"
fi

#PKG-CONFIG
#mkl-dynamic-lp64-seq

LIBS="${LAPACK_LIBS} ${BLAS_LIBS} ${LIBS}"
PKG_LIBS="${PKG_LIBS} ${LIBS}"
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

CC..........: ${CC}
CPPFLAGS....: ${CPPFLAGS}
CFLAGS......: ${CFLAGS}
FC..........: ${FC}
FCLAGS......: ${FCFLAGS}
LDFLAGS:....: ${LDFLAGS}
LIBS........: ${LIBS}

Package features:
  ${ARGS}

Now type 'make @<:@<target>@:>@'
where the optional <target> is:
  all           - build the library
  check         - run tests
  install       - install ${PACKAGE_NAME}
--------------------------------------------------"
