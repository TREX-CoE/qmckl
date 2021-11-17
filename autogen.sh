#!/bin/bash

export srcdir="."
${PYTHON} ${srcdir}/tools/build_makefile.py
autoreconf -i -Wall --no-recursive
