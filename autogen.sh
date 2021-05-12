#!/bin/bash

./tools/build_makefile.sh
autoreconf -i
echo "To configure in maintainer mode, use:
$ QMCKL_DEVEL=1 ./configure --enable-silent-rules --enable-maintainer-mode
"
