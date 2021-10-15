#!/bin/bash

# This scripts is meant to be used by Verificarlo CI to automatically install
# the library dependencies and build it with Verificarlo with vfc_probes support
# enabled.

./autogen.sh
QMCKL_DEVEL=1 ./configure --prefix=$PWD/_install \
--enable-silent-rules --enable-maintainer-mode --enable-vfc_ci --host=x86_64 \
CC="verificarlo-f" FC="verificarlo-f"

make all
# Here we build the test suite, but expect it to fail because it is run without
# specifying VFC_BACKENDS. However, the generated executables will be reused
# individually by the CI.
make check

exit 0
