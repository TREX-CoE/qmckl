#!/bin/bash

# This scripts is meant to be used by Verificarlo CI to automatically install
# the library dependencies and build it with Verificarlo with vfc_probes support
# enabled.

./autogen.sh
# Verificarlo 2.x has one wrapper per language: C is instrumented by
# verificarlo-c and Fortran by verificarlo-f. Using verificarlo-f for C makes
# configure fail its C compiler check.
QMCKL_DEVEL=1 ./configure --prefix=$PWD/_install \
--enable-silent-rules --enable-maintainer-mode --enable-vfc_ci --host=x86_64 \
CC="verificarlo-c" FC="verificarlo-f"

make all
# Here we build the test suite, but expect it to fail because it is run without
# specifying VFC_BACKENDS. However, the generated executables will be reused
# individually by the CI.
#
# The C-only test drivers are linked by the C driver, but libqmckl is built as
# a static library (--host=x86_64 makes autoconf disable shared libraries), so
# its Fortran objects leave the flang runtime symbols unresolved at that step.
# Linking the drivers with the Fortran driver resolves them.
make check CCLD="verificarlo-f"

exit 0
