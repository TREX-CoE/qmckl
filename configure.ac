VERSION=[1.0]
AC_SUBST([VERSION])

AC_REVISION([c2d02f7 Added DESTDIR in Makefile])
AC_INIT([QMCkl],[1.0],
        [https://github.com/TREX-CoE/qmckl/issues], [],
        [https://trex-coe.github.io/qmckl/index.html])

AC_CONFIG_SRCDIR([src/README.org])

AC_LANG_PUSH([C])
AC_PROG_CC([gcc icc clang nvcc xlc])

AC_CHECK_HEADERS([assert.h errno.h malloc.h math.h pthread.h stdbool.h stdint.h stdio.h stdlib.h string.h])

AC_PROG_FC([ifort gfortran flang nvfortran xlf90],[Fortran])
AC_PROG_FC_C_O
AC_FC_SRCEXT([f90])
AC_FC_FREEFORM

AC_CHECK_HEADER([munit/munit.h], [echo found], [echo not found] )

AC_CHECK_LIB([pthread], [pthread_create])

AC_SEARCH_LIBS([dgemm], [blas mkl],
                [],
                AC_MSG_ERROR([Unable to find a BLAS library]) )

AC_CONFIG_FILES(Makefile)



AC_CHECK_PROGS([HAS_EMACS],[emacs],[])

AC_OUTPUT
