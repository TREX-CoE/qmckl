#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
# configure script for maintainers.

AC_PREREQ([2.69])

AC_INIT([QMCkl],[0.1.1],
        [https://github.com/TREX-CoE/qmckl/issues], [],
        [https://trex-coe.github.io/qmckl/index.html])

AC_CONFIG_SRCDIR([src/README.org])
AC_CONFIG_HEADERS([include/config.h])
AC_CONFIG_MACRO_DIR([m4])
#AM_MAINTAINER_MODE

VERSION_MAJOR=`echo ${PACKAGE_VERSION} | cut -d. -f1`
VERSION_MINOR=`echo ${PACKAGE_VERSION} | cut -d. -f2`
VERSION_PATCH=`echo ${PACKAGE_VERSION} | cut -d. -f3 | cut -d- -f1`
AC_DEFINE_UNQUOTED(QMCKL_VERSION_MAJOR, [$VERSION_MAJOR], [major version])
AC_DEFINE_UNQUOTED(QMCKL_VERSION_MINOR, [$VERSION_MINOR], [minor version])
AC_DEFINE_UNQUOTED(QMCKL_VERSION_PATCH, [$VERSION_PATCH], [patch version])

#AM_INIT_AUTOMAKE([foreign color-tests parallel-tests silent-rules 1.11])
#AM_SILENT_RULES(yes)
#AC_SUBST(SHARED_VERSION_INFO)
#AM_ENABLE_SHARED(no) dnl shared libs cause too many headaches to be default
#AC_PROG_LIBTOOL

AC_LANG(C)

# Checks for programs.
AC_PROG_CC
AC_PROG_CC_C_O
AC_PROG_FC
AC_PROG_FC_C_O
AC_FC_SRCEXT([f90])
AC_FC_FREEFORM
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_CC_C99
AC_CHECK_PROGS([EMACS],[emacs26 emacs],[no])
AC_CHECK_PROGS([HAS_CPPCHECK],[cppcheck],[no])
PKG_PROG_PKG_CONFIG()
PKG_LIBS=""
PKG_CFLAGS=""

if test x${EMACS} == xno ; then
  AC_MSG_ERROR([
  --------------------------------------
  Error: Emacs is required for org-mode.
  --------------------------------------
  ])
fi

if test x${HAS_CPPCHECK} == xno ; then
        HAS_CPPCHECK=0
else
        HAS_CPPCHECK=1
fi
AC_SUBST([HAS_CPPCHECK])

# fPIC option
AC_MSG_CHECKING(whether fPIC compiler option is accepted in C)
SAVED_CFLAGS="${CFLAGS}"
CFLAGS="${CFLAGS} -fPIC -Werror"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [return 0;])],
    [AC_MSG_RESULT(yes)
     CFLAGS="${SAVED_CFLAGS} -fPIC"],
    [AC_MSG_RESULT(no)
     CFLAGS="${SAVED_CFLAGS}"])

AC_MSG_CHECKING(whether fPIC compiler option is accepted in Fortran)
SAVED_FCFLAGS="${FCFLAGS}"
FCFLAGS="${FCFLAGS} -fPIC -Werror"
AC_LANG_PUSH(Fortran)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [])],
    [AC_MSG_RESULT(yes)
     FCFLAGS="${SAVED_FCFLAGS} -fPIC"],
    [AC_MSG_RESULT(no)
     FCFLAGS="${SAVED_FCFLAGS}"])
AC_LANG_POP()

# Checks for libraries.

AC_FC_LIBRARY_LDFLAGS

AC_CHECK_LIB([m], [sqrt])
AC_CHECK_HEADERS([assert.h errno.h malloc.h math.h pthread.h stdbool.h stdint.h stdio.h stdlib.h string.h])

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

AC_ARG_WITH(munit, [AC_HELP_STRING([--without-munit],[use munit for testing])], ok=yes, ok=no)
if test "$ok" = "no"; then
AC_CHECK_FILE([munit/munit.h], [],
  AC_MSG_ERROR([
  -------------------------------------------------
  Error: Git submodule for munit is not downloaded.
  $ git submodule update --init
  -------------------------------------------------
]) )
else
  ARGS="${ARGS} munit"
  AC_DEFINE(HAVE_MUNIT,1,[Define when using munit])
fi


# Checks for header files.

## qmckl
AC_CHECK_HEADERS([inttypes.h malloc.h stdint.h stdlib.h string.h strings.h unistd.h])

## munit
AC_CHECK_HEADERS([fcntl.h langinfo.h limits.h locale.h mach/mach.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.

## qmckl
AC_CHECK_HEADER_STDBOOL
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

## munit
AC_TYPE_PID_T
AC_C_INLINE
AC_TYPE_SSIZE_T

# Checks for library functions.

## qmckl
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset strerror])

## munit
AC_FUNC_REALLOC
AC_FUNC_STRERROR_R
AC_FUNC_FORK
AC_CHECK_FUNCS([clock_gettime dup2 gettimeofday nl_langinfo setlocale strcasecmp strtoul])



#PKG-CONFIG
#mkl-dynamic-lp64-seq

LIBS="${LAPACK_LIBS} ${BLAS_LIBS} ${LIBS}"
PKG_LIBS="${PKG_LIBS} ${LIBS}"
AC_SUBST([PKG_LIBS])
AC_SUBST([PKG_CFLAGS])

AC_CONFIG_FILES([Makefile
                qmckl.pc
                configure.ac.dist
		src/Makefile.generated
                src/Makefile])


AC_OUTPUT

echo \
"-------------------------------------------------

${PACKAGE_NAME} Version ${PACKAGE_VERSION}

Prefix: '${prefix}'.

CC:           ${CC}
CPPFLAGS:     ${CPPFLAGS}
CFLAGS:       ${CFLAGS}
FC:           ${FC}
FCLAGS:       ${FCFLAGS}
LDFLAGS:      ${LDFLAGS}
LIBS:         ${LIBS}

Package features:
  ${ARGS}

Now type 'make @<:@<target>@:>@'
where the optional <target> is:
  all           - build the library
  check         - run tests
  install       - install ${PACKAGE_NAME}
--------------------------------------------------"

